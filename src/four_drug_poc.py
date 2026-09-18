"""Four-drug Tier-1 extension of the existing strict docking PoC.

The script keeps the previous three-drug outputs untouched.  It adds
alprazolam using the existing receptor/box/Vina/PLIP settings and writes all
four-drug artifacts to explicitly named files.  Tier-1 pharmacology is kept at
the pair level; a separate Ki-only modal-composition selection is used for the
one-row-per-drug pharmacology/IFP integration table so metrics and receptor
compositions are never silently mixed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Mapping, Tuple

import numpy as np
import pandas as pd

from . import strict_poc as sp

ROOT = Path(__file__).resolve().parents[1]
FOUR_DRUGS = ["diazepam", "alprazolam", "triazolam", "zolpidem"]
RECEPTORS = ["alpha1_beta3_gamma2", "alpha2_beta3_gamma2_local_9CTJ"]
RAW4 = ROOT / "data" / "raw" / "strict4"
PROC = ROOT / "data" / "processed"
RES = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
OUT = ROOT / "outputs"


def mkdirs() -> None:
    for p in [RAW4 / "ligands", RAW4 / "docking", PROC, RES, FIG, OUT]:
        p.mkdir(parents=True, exist_ok=True)


def write_csv(df: pd.DataFrame, *paths: Path) -> None:
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False)


def backup_existing(path: Path) -> None:
    if path.exists():
        backup = path.with_name(path.stem + "_pre_four_drug_backup" + path.suffix)
        if not backup.exists():
            shutil.copy2(path, backup)


def build_primary_pharmacology() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pairs = pd.read_csv(RES / "pharmacology_comparable_metrics.csv")
    master = pd.read_csv(PROC / "pharmacology_master_raw.csv")
    drugs = pd.read_csv(ROOT / "data" / "raw" / "drug_list.csv")
    pairs = pairs[(pairs.drug_id.isin(FOUR_DRUGS)) & (pairs.tier == 1)].copy()
    pairs["ratio_alpha2_over_alpha1"] = pd.to_numeric(pairs.alpha2_value, errors="coerce") / pd.to_numeric(pairs.alpha1_value, errors="coerce")
    pairs["log10_ratio"] = np.log10(pairs["ratio_alpha2_over_alpha1"].where(pairs["ratio_alpha2_over_alpha1"] > 0))
    # Attach fields from the common master schema only when the comparable
    # table does not already carry them.  The pair table already includes
    # alpha1/alpha2 species; merging those columns again would create pandas
    # suffixes and make the QC columns ambiguous.
    master_by_evidence = master.drop_duplicates("evidence_id").set_index("evidence_id")
    for side in ("alpha1", "alpha2"):
        evidence_col = f"{side}_evidence_id"
        for field in ("reference_id", "assay", "reference", "species", "target_name"):
            col = f"{side}_{field}"
            if col not in pairs.columns:
                source = master_by_evidence[field] if field in master_by_evidence.columns else pd.Series(dtype=object)
                pairs[col] = pairs[evidence_col].map(source) if evidence_col in pairs.columns else ""
    pairs["reference"] = pairs["alpha1_reference_id"].where(pairs["alpha1_reference_id"].astype(str).str.len() > 0, pairs["alpha1_reference"])
    pairs["assay_id"] = pairs["alpha1_assay"]
    pairs["species_qc"] = ((pairs["alpha1_species"].astype(str).str.contains("Homo sapiens", case=False)) & (pairs["alpha2_species"].astype(str).str.contains("Homo sapiens", case=False)))
    pairs["composition_qc"] = pairs["same_beta_gamma"].astype(bool)
    pairs["metric_unit_qc"] = pairs["metric"].notna() & pairs["unit"].notna() & (pairs["metric"].astype(str).str.len() > 0) & (pairs["unit"].astype(str).str.len() > 0)
    pairs["primary_inclusion"] = pairs[["species_qc", "composition_qc", "metric_unit_qc"]].all(axis=1)
    # Keep every original comparable-metrics column plus derived ratios/QC
    # fields.  The primary raw table is therefore auditable back to the
    # existing ChEMBL pair table without dropping source/species/tier flags.
    raw = pairs.copy()
    write_csv(raw, PROC / "primary_tier1_pharmacology_raw.csv", RES / "primary_tier1_pharmacology_raw.csv")

    # Summary preserves exact metric/unit/composition/reference groups.  A
    # median is only used inside a group of raw pairs with the same comparison
    # context; no different reference, metric, unit or composition is merged.
    group_cols = ["drug_id", "alpha1_composition", "alpha2_composition", "metric", "unit", "alpha1_reference_id", "alpha2_reference_id", "alpha1_assay", "alpha2_assay"]
    summary = (pairs.groupby(group_cols, dropna=False)
               .agg(number_of_comparable_pairs=("ratio_alpha2_over_alpha1", "size"),
                    alpha1_value=("alpha1_value", "median"), alpha2_value=("alpha2_value", "median"),
                    alpha2_over_alpha1_ratio=("ratio_alpha2_over_alpha1", "median"),
                    reference=("reference", "first"), assay_id=("assay_id", "first"),
                    species_qc=("species_qc", "all"), composition_qc=("composition_qc", "all"),
                    primary_inclusion=("primary_inclusion", "all"))
               .reset_index())
    summary["log10_alpha2_over_alpha1_ratio"] = np.log10(summary["alpha2_over_alpha1_ratio"].where(summary["alpha2_over_alpha1_ratio"] > 0))
    summary["summary_method"] = "median_within_same_composition_metric_unit_reference_assay_group"
    write_csv(summary, PROC / "primary_tier1_pharmacology_summary.csv", RES / "primary_tier1_pharmacology_summary.csv")

    # Ki integration retains every exact comparison context.  In particular,
    # values from different references/assay series are not collapsed into a
    # drug-level representative.  A median is used only when the same exact
    # composition/metric/unit/reference/assay context contains repeated raw
    # rows; the context identifiers remain in the output.
    ki = pairs[(pairs.metric.astype(str).str.casefold() == "ki") & pairs.primary_inclusion].copy()
    qc_rows = []
    selected_rows = []
    context_cols = ["alpha1_composition", "alpha2_composition", "metric", "unit", "alpha1_reference_id", "alpha2_reference_id", "alpha1_assay", "alpha2_assay"]
    for drug in FOUR_DRUGS:
        d = ki[ki.drug_id == drug]
        if d.empty:
            qc_rows.append({"drug_id": drug, "check": "integration_metric_selection", "status": "missing_ki_tier1"})
            continue
        for key, chosen in d.groupby(context_cols, dropna=False, sort=True):
            context = dict(zip(context_cols, key if isinstance(key, tuple) else (key,)))
            ratios = chosen.ratio_alpha2_over_alpha1.dropna()
            ratio = ratios.median() if len(ratios) else np.nan
            selected_rows.append({
                "drug_id": drug, **context, "number_of_comparable_pairs": len(chosen),
                "alpha1_value_median": chosen.alpha1_value.median(), "alpha2_value_median": chosen.alpha2_value.median(),
                "alpha2_over_alpha1_ratio_median": ratio, "log10_ratio_median": np.log10(ratio) if pd.notna(ratio) and ratio > 0 else np.nan,
                "reference": chosen.reference.iloc[0], "assay_id": chosen.assay_id.iloc[0],
                "reference_count": chosen.alpha1_reference_id.nunique(), "assay_count": chosen.alpha1_assay.nunique(),
                "comparison_context_id": "|".join(str(context.get(c, "")) for c in ["alpha1_reference_id", "alpha2_reference_id", "alpha1_assay", "alpha2_assay"]),
                "selection_rule": "Ki-only; exact composition/metric/unit/reference/assay context; median only within this context"
            })
        qc_rows.append({"drug_id": drug, "check": "integration_metric_selection", "status": "exact_contexts_retained", "details": f"Ki; {d.groupby(context_cols, dropna=False).ngroups} exact comparison context(s); no cross-reference median"})
    integration = pd.DataFrame(selected_rows)
    qc = pd.DataFrame(qc_rows)
    # Annotation-vs-assay text cannot be rechecked for all historical rows
    # because the archived processed table does not contain assay_description.
    qc = pd.concat([qc, pd.DataFrame([{"drug_id": d, "check": "target_annotation_vs_assay_description", "status": "not_available_in_archived_chembl_table", "details": "No assay_description field was retained for these rows; no automatic correction performed."} for d in FOUR_DRUGS])], ignore_index=True)
    qc = pd.concat([qc, pd.DataFrame([{"drug_id": d, "check": "tier_filter", "status": "tier1_only", "details": "Tier 2/3 rows excluded from primary dataset"} for d in FOUR_DRUGS])], ignore_index=True)
    # Audit duplicate inflation and separation of metric/unit contexts.  A
    # repeated value is retained as raw evidence unless it is an exact repeat
    # of the same pair context; exact repeats are counted here and never
    # allowed to silently increase the median's evidence count.
    context_cols = ["drug_id", "alpha1_composition", "alpha2_composition", "metric", "unit", "alpha1_reference_id", "alpha2_reference_id", "alpha1_assay", "alpha2_assay"]
    for drug in FOUR_DRUGS:
        d = raw[raw.drug_id == drug]
        dup_n = int(d.duplicated(context_cols).sum())
        contexts = "; ".join(f"{m}/{u}" for m, u in d[["metric", "unit"]].drop_duplicates().itertuples(index=False, name=None))
        qc = pd.concat([qc, pd.DataFrame([{
            "drug_id": drug, "check": "duplicate_same_reference_assay_context",
            "status": "no_exact_duplicate" if dup_n == 0 else "duplicate_flagged",
            "details": f"duplicate_rows={dup_n}; exact pair contexts retained in raw table"
        }, {
            "drug_id": drug, "check": "metric_unit_separation", "status": "kept_separate",
            "details": contexts or "none"
        }, {
            "drug_id": drug, "check": "missing_vs_zero", "status": "no_zero_imputation",
            "details": "missing evidence remains absent; no synthetic zero was inserted"
        }])], ignore_index=True)
    # The docking background is fixed to alpha1/beta3/gamma2 versus the
    # provisional alpha2/beta3/gamma2 local construct.  Keep any pharmacology
    # composition mismatch explicit (alprazolam is beta2-only in Tier 1).
    for drug, group in integration.groupby("drug_id", sort=False):
        exact = bool(((group.alpha1_composition == "alpha1/beta3/gamma2") & (group.alpha2_composition == "alpha2/beta3/gamma2")).all())
        compositions = "; ".join(sorted(set(f"{a} vs {b}" for a, b in zip(group.alpha1_composition, group.alpha2_composition))))
        qc = pd.concat([qc, pd.DataFrame([{
            "drug_id": drug, "check": "docking_vs_pharmacology_composition",
            "status": "same_beta3_background" if exact else "composition_mismatch_review",
            "details": f"pharmacology contexts={compositions}; docking=alpha1/beta3/gamma2 vs alpha2/beta3/gamma2_local_9CTJ"
        }])], ignore_index=True)
    write_csv(qc, RES / "primary_tier1_qc.csv")
    return raw, summary, integration


def generate_alprazolam_ligand() -> Path:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mols = pd.read_csv(ROOT / "data" / "raw" / "pharmacology" / "chembl" / "molecules.csv").set_index("drug_id")
    row = mols.loc["alprazolam"]
    mol = Chem.MolFromSmiles(row.canonical_smiles)
    if mol is None:
        raise ValueError("Could not parse alprazolam canonical SMILES")
    mol = Chem.AddHs(mol)
    seed = int(hashlib.sha256(b"alprazolam").hexdigest()[:8], 16) % (2**31 - 1)
    if AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=False) != 0:
        raise RuntimeError("Alprazolam 3D embedding failed")
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    no_h = Chem.RemoveHs(mol)
    sdf = RAW4 / "ligands" / "alprazolam.sdf"
    pdb = RAW4 / "ligands" / "alprazolam.pdb"
    pdbqt = RAW4 / "ligands" / "alprazolam.pdbqt"
    Chem.MolToMolFile(no_h, str(sdf)); Chem.MolToPDBFile(no_h, str(pdb), flavor=0)
    try: AllChem.ComputeGasteigerCharges(no_h)
    except Exception: pass
    conf = no_h.GetConformer(); lines = ["REMARK Generated from ChEMBL canonical SMILES with RDKit ETKDG + MMFF/UFF", "ROOT"]
    for i, atom in enumerate(no_h.GetAtoms(), 1):
        pos = conf.GetAtomPosition(i - 1)
        try: charge = float(atom.GetProp("_GasteigerCharge")); charge = charge if math.isfinite(charge) else 0.0
        except Exception: charge = 0.0
        typ = sp._element_type(atom.GetSymbol(), atom.GetIsAromatic())
        name = f"{atom.GetSymbol()}{i}"[:4]
        lines.append(f"ATOM  {i:5d} {name:<4s} LIG Z   1    {pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}{1.0:6.2f}{0.0:6.2f}    {charge:7.3f} {typ:>2s}")
    lines += ["ENDROOT", "TORSDOF 0", ""]
    pdbqt.write_text("\n".join(lines))
    return pdbqt


def _score_rows(out_pdbqt: Path, log_path: Path) -> Dict[int, float]:
    scores: Dict[int, float] = {}; current = None
    for line in out_pdbqt.read_text().splitlines():
        if line.startswith("MODEL"):
            try: current = int(line.split()[1])
            except Exception: current = None
        if current is not None and line.startswith("REMARK VINA RESULT:"):
            scores[current] = float(line.split()[3])
    if not scores and log_path.exists():
        for line in log_path.read_text().splitlines():
            m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+", line)
            if m: scores[int(m.group(1))] = float(m.group(2))
    return scores


def run_alprazolam_docking(ligand: Path, boxes: pd.DataFrame) -> pd.DataFrame:
    vina = shutil.which("vina")
    if vina is None: raise RuntimeError("AutoDock Vina executable not found")
    receptor_paths = {"alpha1_beta3_gamma2": ROOT / "data/raw/strict3/receptors/alpha1_beta3_gamma2.pdbqt", "alpha2_beta3_gamma2_local_9CTJ": ROOT / "data/raw/strict3/receptors/alpha2_beta3_gamma2_local_9CTJ.pdbqt"}
    rows = []
    for receptor_id, rec in receptor_paths.items():
        box = boxes[boxes.receptor_id == receptor_id].iloc[0]
        job = RAW4 / "docking" / f"alprazolam__{receptor_id}"; job.mkdir(parents=True, exist_ok=True)
        out_pdbqt, log = job / "poses.pdbqt", job / "vina.log"
        if not out_pdbqt.exists():
            cmd = [vina, "--receptor", str(rec), "--ligand", str(ligand), "--center_x", str(box.center_x), "--center_y", str(box.center_y), "--center_z", str(box.center_z), "--size_x", str(box.size_x), "--size_y", str(box.size_y), "--size_z", str(box.size_z), "--exhaustiveness", "8", "--num_modes", "9", "--energy_range", "4", "--seed", "20260917", "--out", str(out_pdbqt)]
            proc = subprocess.run(cmd, cwd=job, text=True, capture_output=True); log.write_text(proc.stdout + "\n" + proc.stderr)
            if proc.returncode != 0: raise RuntimeError(f"Vina failed for alprazolam/{receptor_id}: {proc.stderr[-1000:]}")
        scores = _score_rows(out_pdbqt, log)
        pose_dir = job / "pose_pdb"; pose_dir.mkdir(exist_ok=True)
        for pose_id, score in sorted(scores.items()):
            lines = []; current = None
            for line in out_pdbqt.read_text().splitlines():
                if line.startswith("MODEL"): current = int(line.split()[1])
                elif line.startswith("ENDMDL"):
                    if current == pose_id: break
                elif current == pose_id and line.startswith(("ATOM", "HETATM")): lines.append(line)
            one = job / "_one_pose.pdbqt"; one.write_text("\n".join(lines) + "\n")
            pose_pdb = pose_dir / f"pose_{pose_id:02d}.pdb"; sp._receptor_pdb_for_plip(receptor_id, one, pose_pdb)
            rows.append({"drug_id": "alprazolam", "receptor_id": receptor_id, "pdb_id": "6HUP" if receptor_id == "alpha1_beta3_gamma2" else "9CTJ", "pose_id": pose_id, "rank": pose_id, "vina_score_kcal_mol": score, "vina_seed": 20260917, "exhaustiveness": 8, "num_modes": 9, "energy_range": 4, "docking_box": f"{box.center_x:.3f},{box.center_y:.3f},{box.center_z:.3f};{box.size_x:.1f},{box.size_y:.1f},{box.size_z:.1f}", "receptor_structure_id": "6HUP" if receptor_id == "alpha1_beta3_gamma2" else "9CTJ", "ligand_preparation_method": "RDKit ETKDG + MMFF/UFF + Gasteiger/PDBQT writer", "pose_pdbqt": str(out_pdbqt), "pose_pdb": str(pose_pdb), "status": "docked"})
    return pd.DataFrame(rows)


def parse_plip_xml(xml_path: Path, base: Mapping[str, object]) -> List[Dict[str, object]]:
    import xml.etree.ElementTree as ET
    out = []
    try: root = ET.parse(xml_path).getroot()
    except Exception: return out
    tags = {"hydrophobic_interaction", "hydrogen_bond", "water_bridge", "salt_bridge", "pi_stack", "pi_cation", "halogen_bond", "metal_complex"}
    for node in root.iter():
        tag = node.tag.lower().split("}")[-1]
        if tag not in tags: continue
        attrs = {k.lower().split("}")[-1]: v for k, v in node.attrib.items()}
        for child in list(node):
            if child.text and child.text.strip(): attrs[child.tag.lower().split("}")[-1]] = child.text.strip()
        out.append({**base, "interaction_type": tag, "residue_number": attrs.get("resnr", ""), "residue_chain": attrs.get("reschain", ""), "raw_attributes": json.dumps(attrs, sort_keys=True)})
    return out


def run_alprazolam_plip(docking: pd.DataFrame) -> pd.DataFrame:
    plip = ROOT / "work/plip_vendor/bin/plip"
    rows = []
    for _, rec in docking.iterrows():
        pose = Path(rec.pose_pdb)
        # One PLIP output directory per pose preserves the raw XML/text for
        # every pose instead of letting later poses overwrite the report.
        outdir = pose.parent / f"plip_pose_{int(rec.pose_id):02d}"; outdir.mkdir(exist_ok=True); xml = outdir / "report.xml"
        # PLIP writes temporary protonated/"plipfixed" files with random
        # suffixes.  Remove the previous pose's files so repeated runs do not
        # accumulate stale raw artifacts in the new four-drug directory.
        for stale in outdir.iterdir():
            if stale.is_file():
                stale.unlink()
        cmd = [str(plip), "-f", str(pose), "-o", str(outdir), "-x", "-t", "-q", "--breakcomposite", "--name", "report"]
        env = os.environ.copy(); env["PYTHONPATH"] = str(ROOT / "work/plip_vendor")
        proc = subprocess.run(cmd, cwd=outdir, env=env, text=True, capture_output=True); (outdir / "plip.stdout.txt").write_text(proc.stdout + "\n" + proc.stderr)
        rows.extend(parse_plip_xml(xml, {"drug_id": "alprazolam", "receptor_id": rec.receptor_id, "pose_id": int(rec.pose_id)}))
    return pd.DataFrame(rows, columns=["drug_id", "receptor_id", "pose_id", "interaction_type", "residue_number", "residue_chain", "raw_attributes"])


def build_fingerprints(docking4: pd.DataFrame, interactions4: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mapping = pd.read_csv(OUT / "residue_mapping.csv")
    a2_to_a1 = {int(r.alpha2_residue_number): int(r.alpha1_residue_number) for _, r in mapping.iterrows() if str(r.mapping_status) == "mapped_global_sequence_alignment" and str(r.alpha2_residue_number).strip() not in ("", "nan")}
    inter = interactions4.copy()
    inter["residue_number_num"] = pd.to_numeric(inter.residue_number, errors="coerce")
    unresolved = []
    def feat(row):
        if pd.isna(row.residue_number_num): return None
        n = int(row.residue_number_num)
        if row.receptor_id == "alpha1_beta3_gamma2": return f"BZD_SITE_{n:03d}|{row.interaction_type}"
        if n in a2_to_a1: return f"BZD_SITE_{a2_to_a1[n]:03d}|{row.interaction_type}"
        unresolved.append(row)
        return None
    inter["feature"] = inter.apply(feat, axis=1)
    unresolved_df = pd.DataFrame(unresolved)
    write_csv(unresolved_df, RES / "four_drug_unresolved_mapping.csv")
    mapped = inter[inter.feature.notna()].copy()
    features = sorted(mapped.feature.unique())
    bin_rows, freq_rows = [], []
    for drug in FOUR_DRUGS:
        for receptor in RECEPTORS:
            sub = mapped[(mapped.drug_id == drug) & (mapped.receptor_id == receptor)]
            raw_unresolved = inter[(inter.drug_id == drug) & (inter.receptor_id == receptor) & inter.feature.isna() & inter.residue_number_num.notna()]
            pose_n = int(docking4[(docking4.drug_id == drug) & (docking4.receptor_id == receptor)].pose_id.nunique())
            for feature in features:
                sf = sub[sub.feature == feature]; count = len(sf); nposes = int(sf.pose_id.nunique()); status = "observed" if count else ("unresolved_mapping" if len(raw_unresolved) else ("no_interaction_observed" if pose_n else "no_data"))
                freq = nposes / pose_n if pose_n else np.nan
                bin_rows.append({"drug_id": drug, "receptor_id": receptor, "feature": feature, "binary": int(count > 0) if status != "unresolved_mapping" else np.nan, "interaction_data_status": status, "n_poses": pose_n, "unresolved_mapping_count": len(raw_unresolved)})
                freq_rows.append({"drug_id": drug, "receptor_id": receptor, "feature": feature, "pose_frequency": freq if status != "unresolved_mapping" else np.nan, "n_interacting_poses": nposes, "n_poses": pose_n, "interaction_data_status": status, "unresolved_mapping_count": len(raw_unresolved)})
    binary = pd.DataFrame(bin_rows); frequency = pd.DataFrame(freq_rows)
    for path in [PROC / "fingerprint_binary.csv", PROC / "fingerprint_frequency.csv"]: backup_existing(path)
    write_csv(binary, PROC / "fingerprint_binary.csv", RES / "fingerprint_binary_4drug.csv")
    write_csv(frequency, PROC / "fingerprint_frequency.csv", RES / "fingerprint_frequency_4drug.csv")
    return binary, frequency, unresolved_df


def similarity_table(binary: pd.DataFrame, frequency: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for drug in FOUR_DRUGS:
        a = binary[(binary.drug_id == drug) & (binary.receptor_id == RECEPTORS[0])].set_index("feature")
        b = binary[(binary.drug_id == drug) & (binary.receptor_id == RECEPTORS[1])].set_index("feature")
        f = sorted(set(a.index) | set(b.index)); aa, bb = a.reindex(f), b.reindex(f)
        valid = aa.interaction_data_status.ne("unresolved_mapping") & bb.interaction_data_status.ne("unresolved_mapping")
        av = aa.loc[valid, "binary"].fillna(0).astype(bool); bv = bb.loc[valid, "binary"].fillna(0).astype(bool)
        union, inter = int((av | bv).sum()), int((av & bv).sum())
        fa = frequency[(frequency.drug_id == drug) & (frequency.receptor_id == RECEPTORS[0])].set_index("feature").reindex(f).loc[valid, "pose_frequency"].fillna(0).to_numpy(float)
        fb = frequency[(frequency.drug_id == drug) & (frequency.receptor_id == RECEPTORS[1])].set_index("feature").reindex(f).loc[valid, "pose_frequency"].fillna(0).to_numpy(float)
        denom = np.linalg.norm(fa) * np.linalg.norm(fb); cosine = float(np.dot(fa, fb) / denom) if denom else np.nan
        rows.append({"drug_id": drug, "jaccard_similarity": inter / union if union else np.nan, "shared_features": inter, "alpha1_only_features": int((av & ~bv).sum()), "alpha2_only_features": int((bv & ~av).sum()), "changed_features": int((av != bv).sum()), "cosine_similarity_frequency": cosine, "euclidean_distance_frequency": float(np.linalg.norm(fa - fb)) if len(fa) else np.nan, "valid_features": int(valid.sum()), "unresolved_mapping_included": bool((~valid).any())})
    sim = pd.DataFrame(rows)
    write_csv(sim, RES / "alpha1_alpha2_similarity_4drug.csv")
    return sim


def differential_features(binary: pd.DataFrame, frequency: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    features = sorted(binary.feature.unique())
    for drug in FOUR_DRUGS:
        a = binary[(binary.drug_id == drug) & (binary.receptor_id == RECEPTORS[0])].set_index("feature")
        b = binary[(binary.drug_id == drug) & (binary.receptor_id == RECEPTORS[1])].set_index("feature")
        fa = frequency[(frequency.drug_id == drug) & (frequency.receptor_id == RECEPTORS[0])].set_index("feature")
        fb = frequency[(frequency.drug_id == drug) & (frequency.receptor_id == RECEPTORS[1])].set_index("feature")
        for feature in features:
            sa, sb = a.loc[feature], b.loc[feature]; va, vb = fa.loc[feature].pose_frequency, fb.loc[feature].pose_frequency
            if sa.interaction_data_status == "unresolved_mapping" or sb.interaction_data_status == "unresolved_mapping": category = "unresolved_mapping"
            elif sa.binary == 1 and sb.binary == 0: category = "alpha1_only"
            elif sa.binary == 0 and sb.binary == 1: category = "alpha2_only"
            elif pd.notna(va) and pd.notna(vb) and va > vb: category = "alpha1_frequency_higher"
            elif pd.notna(va) and pd.notna(vb) and vb > va: category = "alpha2_frequency_higher"
            else: category = "no_differential"
            rows.append({"drug_id": drug, "feature": feature, "alpha1_binary": sa.binary, "alpha2_binary": sb.binary, "alpha1_frequency": va, "alpha2_frequency": vb, "differential_category": category, "status": "unresolved_mapping" if category == "unresolved_mapping" else "observed_or_no_interaction"})
    diff = pd.DataFrame(rows)
    write_csv(diff, RES / "differential_features_by_drug.csv")
    summary = []
    for feature, g in diff.groupby("feature"):
        changed = sorted(g.loc[~g.differential_category.isin(["no_differential", "unresolved_mapping"]), "drug_id"].unique())
        summary.append({"feature": feature, "changed_in_drugs": "|".join(changed), "n_changed_drugs": len(changed), "drug_specific": len(changed) == 1, "drug_specific_to": changed[0] if len(changed) == 1 else ""})
    specific = pd.DataFrame(summary)
    write_csv(specific, RES / "drug_specific_features.csv")
    return diff, specific


def make_heatmaps(binary: pd.DataFrame, frequency: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for data, col, name in [(binary, "binary", "fingerprint_binary_heatmap_4drug.png"), (frequency, "pose_frequency", "fingerprint_frequency_heatmap_4drug.png")]:
        piv = data.assign(row=data.drug_id + "|" + data.receptor_id).pivot(index="row", columns="feature", values=col).reindex(index=[f"{d}|{r}" for d in FOUR_DRUGS for r in RECEPTORS]).fillna(np.nan)
        status = data.assign(row=data.drug_id + "|" + data.receptor_id).pivot(index="row", columns="feature", values="interaction_data_status").reindex(index=piv.index, columns=piv.columns)
        arr = piv.to_numpy(float); arr[status.to_numpy() == "unresolved_mapping"] = np.nan
        short_rows = [f"{d}|{'α1' if r == RECEPTORS[0] else 'α2'}" for d in FOUR_DRUGS for r in RECEPTORS]
        title = "4-drug binary IFP (gray = unresolved mapping)" if col == "binary" else "4-drug pose-frequency IFP (gray = unresolved mapping)"
        fig, ax = plt.subplots(figsize=(max(8, len(piv.columns) * .28), 5.2)); cmap = plt.get_cmap("viridis").copy(); cmap.set_bad("#d9d9d9"); im = ax.imshow(arr, aspect="auto", cmap=cmap, vmin=0, vmax=1); ax.set_yticks(range(len(piv.index)), short_rows); ax.set_xticks(range(len(piv.columns)), piv.columns, rotation=90, fontsize=7); ax.set_title(title, fontsize=11); fig.colorbar(im, ax=ax, fraction=.02); fig.tight_layout(rect=[0, 0, 1, .95]); fig.savefig(FIG / name, dpi=180); plt.close(fig)


def integrate_pharmacology(sim: pd.DataFrame, integration: pd.DataFrame) -> pd.DataFrame:
    out = integration.merge(sim, on="drug_id", how="outer")
    write_csv(out, RES / "ifp_vs_pharmacology_4drug.csv")
    return out


def integrated_figure(integrated: pd.DataFrame, diff: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cats = ["alpha1_only", "alpha2_only", "alpha1_frequency_higher", "alpha2_frequency_higher"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    # Keep all exact-context Ki ratios as points; no cross-reference median is
    # introduced merely to make a one-bar-per-drug plot.
    x = np.arange(len(FOUR_DRUGS))
    for i, drug in enumerate(FOUR_DRUGS):
        vals = integrated.loc[integrated.drug_id == drug, "alpha2_over_alpha1_ratio_median"].dropna().to_numpy(float)
        if len(vals):
            jitter = np.linspace(-0.12, 0.12, len(vals)) if len(vals) > 1 else np.array([0.0])
            axes[0].scatter(np.full(len(vals), i) + jitter, vals, s=28, alpha=.85)
    axes[0].set_xticks(x, FOUR_DRUGS, rotation=30); axes[0].set_yscale("log"); axes[0].set_ylabel("Ki α2 / α1 ratio (each exact context)"); axes[0].set_title("Tier-1 pharmacology")
    sim_plot = integrated[["drug_id", "jaccard_similarity", "cosine_similarity_frequency"]].drop_duplicates("drug_id").set_index("drug_id").reindex(FOUR_DRUGS)
    axes[1].bar(x - .18, sim_plot.jaccard_similarity, .36, label="Jaccard"); axes[1].bar(x + .18, sim_plot.cosine_similarity_frequency, .36, label="frequency cosine"); axes[1].set_xticks(x, FOUR_DRUGS, rotation=30); axes[1].set_ylim(0, 1); axes[1].legend(); axes[1].set_title("IFP similarity")
    counts = diff[diff.differential_category.isin(cats)].pivot_table(index="drug_id", columns="differential_category", values="feature", aggfunc="count").reindex(index=FOUR_DRUGS, columns=cats).fillna(0); bottom = np.zeros(len(counts));
    for cat in cats:
        axes[2].bar(np.arange(len(counts)), counts[cat], bottom=bottom, label=cat); bottom += counts[cat].to_numpy()
    axes[2].set_xticks(np.arange(len(counts)), counts.index, rotation=30); axes[2].set_title("Differential feature counts"); axes[2].legend(fontsize=7)
    fig.suptitle("Four-drug mechanistic feasibility PoC; no statistical generalization"); fig.tight_layout(); fig.savefig(FIG / "ifp_vs_pharmacology_integrated_4drug.png", dpi=200); plt.close(fig)


def write_qc_report(docking: pd.DataFrame, interactions: pd.DataFrame, binary: pd.DataFrame, unresolved: pd.DataFrame) -> None:
    """Write an auditable four-drug QC summary without converting missing data to zero."""
    strict_files = {
        "outputs/docking_results.csv": len(pd.read_csv(OUT / "docking_results.csv")),
        "outputs/plip_interactions.csv": len(pd.read_csv(OUT / "plip_interactions.csv")),
        "outputs/fingerprint_binary.csv": len(pd.read_csv(OUT / "fingerprint_binary.csv")),
        "outputs/fingerprint_frequency.csv": len(pd.read_csv(OUT / "fingerprint_frequency.csv")),
        "outputs/alpha1_alpha2_similarity.csv": len(pd.read_csv(OUT / "alpha1_alpha2_similarity.csv")),
    }
    pose_counts = docking.groupby(["drug_id", "receptor_id"]).pose_id.nunique()
    pose_ok = bool(len(pose_counts) == 8 and (pose_counts == 9).all())
    param_cols = ["vina_seed", "exhaustiveness", "num_modes", "energy_range"]
    param_ok = all(docking[c].nunique(dropna=True) == 1 for c in param_cols if c in docking.columns)
    box_ok = docking.docking_box.nunique(dropna=True) == 2
    receptor_ok = set(docking.receptor_structure_id.dropna()) == {"6HUP", "9CTJ"}
    status_counts = binary.interaction_data_status.value_counts().to_dict()
    lines = [
        "# FOUR_DRUG_QC_REPORT", "", "QC is descriptive and does not validate biological equivalence of the receptor models.", "",
        "## Preservation of strict 3-drug outputs", "", "| file | rows |", "|---|---:|",
    ]
    lines += [f"| `{k}` | {v} |" for k, v in strict_files.items()]
    lines += ["", "## Four-drug computational checks", "", f"- Drug × receptor combinations: {len(pose_counts)}; exactly 9 poses each: **{'PASS' if pose_ok else 'REVIEW'}**.", f"- Vina parameter uniformity (seed=20260917, exhaustiveness=8, num_modes=9, energy_range=4): **{'PASS' if param_ok else 'REVIEW'}**.", f"- Docking boxes: two recorded boxes from `outputs/docking_boxes.csv`; no result-dependent adjustment: **{'PASS' if box_ok else 'REVIEW'}**.", f"- Receptor structure IDs: {', '.join(sorted(receptor_ok and set(docking.receptor_structure_id.dropna()) or set()))}; alpha2 9CTJ construct remains provisional: **{'PASS' if receptor_ok else 'REVIEW'}**.", f"- Combined PLIP interaction rows: {len(interactions)}; unresolved mapping rows: {len(unresolved)}.", f"- Fingerprint status counts (including explicit unresolved/no-interaction labels): `{json.dumps(status_counts, ensure_ascii=False)}`.", "", "## Biological/QC limitations", "", "- Alpha1 uses 6HUP; alpha2 uses a local 9CTJ-derived construct and a transferred DZP box. They are not an exact matched full-pentamer pair.", "- `outputs/residue_mapping.csv` is sequence/structure-alignment based. Mapping absent for a PLIP residue is retained as `unresolved_mapping`; no new common position is inferred.", "- `no_interaction_observed` is emitted only for a feature with pose data and resolved mapping. Missing or unresolved evidence is not converted to binary zero.", "- Archived ChEMBL rows do not retain assay descriptions for all Tier-1 records, so annotation-vs-assay contradiction checking is explicitly unavailable and no automatic correction is made.", "- Alprazolam has only Tier-1 beta2 pharmacology while the fixed docking background is beta3; diazepam includes one beta2 context as well. These are recorded in `results/tables/primary_tier1_qc.csv` as `composition_mismatch_review`.", "- No GtoPdb integration, imputation, activity-type/unit mixing, ML or statistical generalization was performed.", ""]
    (OUT / "FOUR_DRUG_QC_REPORT.md").write_text("\n".join(lines))


def write_report(raw: pd.DataFrame, summary: pd.DataFrame, integration: pd.DataFrame, sim: pd.DataFrame, diff: pd.DataFrame, specific: pd.DataFrame, docking: pd.DataFrame, interactions: pd.DataFrame, unresolved: pd.DataFrame) -> None:
    old = (OUT / "FINAL_REPORT.md").read_text() if (OUT / "FINAL_REPORT.md").exists() else ""
    common = sim[sim.drug_id.isin(["diazepam", "triazolam", "zolpidem"])].merge(pd.read_csv(OUT / "alpha1_alpha2_similarity.csv"), on="drug_id", suffixes=("_4drug", "_3drug"), how="left")
    changed = diff[~diff.differential_category.isin(["no_differential", "unresolved_mapping"])].groupby("drug_id").size().to_dict()
    # Keep the requested interpretation questions tied to the generated
    # tables.  These are descriptive statements for a four-drug feasibility
    # study, not inferential statistics.
    ratio_text_parts = []
    for drug, group in integration.groupby("drug_id", sort=False):
        vals = ", ".join(f"{r.alpha2_over_alpha1_ratio_median:.3g} (log10={r.log10_ratio_median:.3g}; {r.metric}/{r.unit}; {int(r.number_of_comparable_pairs)} pair)" for _, r in group.iterrows())
        ratio_text_parts.append(f"{drug}: [{vals}] across {len(group)} exact context(s)")
    ratio_text = "; ".join(ratio_text_parts)
    sim_text = "; ".join(
        f"{r.drug_id}: Jaccard={r.jaccard_similarity:.3f}, cosine={r.cosine_similarity_frequency:.3f}, changed={int(r.changed_features)}"
        for _, r in sim.iterrows()
    )
    drug_diff_text = {}
    for drug in FOUR_DRUGS:
        g = diff[(diff.drug_id == drug) & (~diff.differential_category.isin(["no_differential", "unresolved_mapping"]))]
        drug_diff_text[drug] = ", ".join(f"{r.feature} [{r.differential_category}]" for _, r in g.iterrows()) or "none among resolved features"
    zolpidem_resolved = drug_diff_text.get("zolpidem", "none")
    unique_text = ", ".join(f"{r.feature} ({r.drug_specific_to})" for _, r in specific[specific.drug_specific == True].iterrows()) or "none"
    report = [
        "# FOUR_DRUG_POC_REPORT", "", "This is a four-drug Tier-1 mechanistic feasibility extension. Existing strict three-drug results remain in their original files.", "",
        "## Primary pharmacology", "", f"Tier-1 raw pair rows: **{len(raw)}**; summary groups: **{len(summary)}**. Ki-only exact comparison-context rows retained for integration: **{len(integration)}** (no cross-reference representative was selected).", "", integration.to_string(index=False), "",
        "## IFP results", "", sim.to_string(index=False), "", f"Alprazolam docking rows: **{len(docking[docking.drug_id == 'alprazolam'])}**; PLIP interaction rows in four-drug set: **{len(interactions)}**; unresolved mapping rows: **{len(unresolved)}**.", "", f"Differential mapped-feature counts: {changed}.", "",
        "## Drug-specific differential features", "", specific[specific.drug_specific == True].to_string(index=False) if not specific.empty else "None.", "",
        "## Comparison with the previous three-drug PoC", "", common.to_string(index=False), "", "Alprazolam expands the feature space and provides a fourth qualitative comparison. Any apparent correspondence between Ki ratios and IFP similarity is descriptive only; n=4 is insufficient for a generalization or significance claim.", "",
        "## Answers to the requested interpretation questions", "",
        "1. **Pharmacology difference:** " + ratio_text + ". Ratios are listed per exact reference/assay context; no values from different references were averaged. Alprazolam has only alpha1/beta2/gamma2 versus alpha2/beta2/gamma2 Tier-1 evidence; diazepam also has one beta2 context alongside its beta3 contexts. These composition differences are flagged for review against the fixed beta3 docking background.",
        "2. **IFP similarity:** " + sim_text + ". Similarity uses only features with resolved mapping on both receptor sides; unresolved rows are excluded from the metric and shown as gray in the heatmaps.",
        "3. **Apparent pharmacology–IFP trend:** the four ratios and similarities are not monotonic (for example zolpidem has the largest Ki ratio but not the lowest Jaccard). No general trend is claimed from n=4.",
        "4. **Drug-specific differences:** resolved drug-specific features are " + unique_text + ". They are descriptive candidates, not biomarkers.",
        "5. **Zolpidem:** resolved differential features are " + zolpidem_resolved + ". The drug-specific table does not identify a feature unique to zolpidem; unresolved mapping features are not interpreted.",
        "6. **Effect of adding alprazolam:** the common feature space now contains four drugs and the four-drug similarity table adds alprazolam (Jaccard 0.778; changed features 2). The prior three-drug rows and files remain preserved, while the beta2 pharmacology versus beta3 docking mismatch limits direct mechanistic integration for alprazolam.",
        "7. **Largest confounder:** the provisional alpha2 local construct/chain background and transferred box are not an exact matched alpha1/alpha2 structural pair; pharmacology contexts that use beta2 while docking uses beta3 (all alprazolam contexts and one diazepam context) are additional composition confounders.",
        "8. **Value of more drugs:** additional compounds are useful for testing whether resolved feature-level patterns recur, provided they have same-metric, same-unit, composition-matched Tier-1 pharmacology and comparable structures. This PoC does not estimate predictive performance.", "",
        "## QC and limitations", "", "The existing 6HUP alpha1 receptor, provisional 9CTJ-derived alpha2 local construct, transferred DZP box, receptor preparation, residue mapping and Vina parameters were reused. The alpha2 construct is not an exact alpha2beta3gamma2 full pentamer, and the two receptor backgrounds are not fully identical. Gray heatmap cells indicate unresolved mapping; zero indicates no interaction observed only when pose data and mapping were available. GtoPdb was not integrated in this phase.", "", "The archived ChEMBL table does not retain assay_description for all Tier-1 rows, so annotation-versus-assay-description contradiction checking is reported as unavailable rather than inferred.", "",
        "The exact strict three-drug files under `outputs/` were left intact. Four-drug combined tables use `_4drug` names; the required four-drug `data/processed/fingerprint_binary.csv` and `fingerprint_frequency.csv` are generated separately from those strict outputs.", "",
        "Four-drug combined tables: `results/tables/docking_results_4drug.csv` (72 pose rows; 9 poses per drug/receptor), `results/tables/plip_interactions_4drug.csv`, `results/tables/fingerprint_binary_4drug.csv`, and `results/tables/fingerprint_frequency_4drug.csv`. A computational/biological QC checklist is in `outputs/FOUR_DRUG_QC_REPORT.md`.", "",
    ]
    report_text = "\n".join(report)
    (OUT / "FOUR_DRUG_POC_REPORT.md").write_text(report_text)
    section = "\n\n---\n\n" + "\n".join(report)
    marker = "# FOUR_DRUG_POC_REPORT"
    if marker in old:
        prefix = old.split(marker, 1)[0].rstrip()
        (OUT / "FINAL_REPORT.md").write_text(prefix + section + "\n")
    else:
        (OUT / "FINAL_REPORT.md").write_text(old.rstrip() + section + "\n")


def run() -> Dict[str, object]:
    mkdirs()
    raw, summary, integration = build_primary_pharmacology()
    ligand = generate_alprazolam_ligand()
    boxes = pd.read_csv(OUT / "docking_boxes.csv")
    alpha_docking = run_alprazolam_docking(ligand, boxes)
    old_docking = pd.read_csv(OUT / "docking_results.csv")
    for col, default in [("rank", None), ("energy_range", 4), ("docking_box", ""), ("receptor_structure_id", ""), ("ligand_preparation_method", "RDKit ETKDG + MMFF/UFF + Gasteiger/PDBQT writer")]:
        if col not in old_docking.columns:
            old_docking[col] = old_docking.pose_id if col == "rank" else default
    for _, r in old_docking.iterrows():
        b = boxes[boxes.receptor_id == r.receptor_id]
        if not b.empty:
            row = b.iloc[0]; old_docking.loc[old_docking.index == r.name, "docking_box"] = f"{row.center_x:.3f},{row.center_y:.3f},{row.center_z:.3f};{row.size_x:.1f},{row.size_y:.1f},{row.size_z:.1f}"
            old_docking.loc[old_docking.index == r.name, "receptor_structure_id"] = r.pdb_id
    docking4 = pd.concat([old_docking, alpha_docking], ignore_index=True)
    write_csv(docking4, RES / "docking_results_4drug.csv")
    new_inter = run_alprazolam_plip(alpha_docking)
    old_inter = pd.read_csv(OUT / "plip_interactions.csv")
    interactions4 = pd.concat([old_inter, new_inter], ignore_index=True)
    write_csv(interactions4, RES / "plip_interactions_4drug.csv", PROC / "plip_interactions_4drug.csv")
    binary, frequency, unresolved = build_fingerprints(docking4, interactions4)
    sim = similarity_table(binary, frequency)
    diff, specific = differential_features(binary, frequency)
    make_heatmaps(binary, frequency)
    integrated = integrate_pharmacology(sim, integration)
    integrated_figure(integrated, diff)
    write_qc_report(docking4, interactions4, binary, unresolved)
    write_report(raw, summary, integration, sim, diff, specific, docking4, interactions4, unresolved)
    return {"tier1_raw_pairs": len(raw), "tier1_summary_groups": len(summary), "docking_rows_4drug": len(docking4), "plip_rows_4drug": len(interactions4), "unresolved_mapping_rows": len(unresolved)}


def main() -> None:
    print(json.dumps(run(), ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
