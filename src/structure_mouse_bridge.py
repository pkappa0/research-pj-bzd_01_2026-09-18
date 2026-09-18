"""Bridge residue-level structural features to a transparent mouse phenotype evidence table.

This phase consumes the lossless PLIP and docking outputs from the corrected
residue run.  It deliberately does not add interactions that are absent from
the raw PLIP table and does not convert qualitative animal observations into
numeric scores.  A new immutable run directory is created on every invocation.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from . import run_manager
except ImportError:  # pragma: no cover
    import run_manager


ROOT = Path(__file__).resolve().parents[1]
DRUGS = ["diazepam", "alprazolam", "triazolam", "zolpidem"]
RECEPTORS = {
    "alpha1": "alpha1_beta3_gamma2",
    "alpha2": "alpha2_beta3_gamma2_local_9CTJ",
}
SOURCE_RUN = "20260917_1629_residue_hypothesis"

# These positions are the requested BZD-site candidates, expressed using the
# alpha1/common-position numbering in the corrected alignment.  A row is used
# only when it occurs in the re-mapped raw PLIP data.
CANDIDATE_POSITIONS = {
    "BZD_GAMMA2_058", "BZD_GAMMA2_060", "BZD_GAMMA2_077",
    "BZD_SITE_100", "BZD_SITE_101", "BZD_SITE_102",
    "BZD_SITE_156", "BZD_SITE_157", "BZD_SITE_205",
    "BZD_SITE_206", "BZD_SITE_210", "BZD_SITE_211",
}


def _num(x: Any) -> float:
    try:
        y = float(x)
        return y if math.isfinite(y) else np.nan
    except Exception:
        return np.nan


def _attrs(x: Any) -> dict[str, Any]:
    try:
        v = json.loads(str(x))
        return v if isinstance(v, dict) else {}
    except Exception:
        return {}


def _stats(vals: list[float]) -> tuple[float, float, float, int]:
    x = pd.Series(vals, dtype="float64").dropna()
    if x.empty:
        return np.nan, np.nan, np.nan, 0
    return float(x.mean()), float(x.median()), float(x.std(ddof=1)) if len(x) > 1 else 0.0, int(len(x))


def _source_dir() -> Path:
    p = run_manager.RUNS / SOURCE_RUN
    if not (p / "raw" / "plip_interactions_4drug.csv").exists():
        raise FileNotFoundError(f"Expected corrected source run is missing: {p}")
    return p


def _copy_lineage(run_dir: Path, source: Path, manifest: dict[str, Any]) -> None:
    files = [
        source / "raw" / "plip_interactions_4drug.csv",
        source / "raw" / "docking_results_4drug.csv",
        source / "tables" / "corrected_residue_mapping.csv",
        source / "tables" / "recovered_interactions_qc.csv",
        source / "tables" / "residue_interaction_features_alpha1.csv",
        source / "tables" / "residue_interaction_features_alpha2.csv",
        source / "tables" / "aromatic_interaction_summary.csv",
    ]
    for f in files:
        if f.exists():
            rel = Path("data") / "raw" / "lineage" / f.name
            run_manager.copy_into_run(f, run_dir / rel)
    manifest["input_files"] = sorted(set(manifest["input_files"] + [str(p.relative_to(run_dir)) for p in (run_dir / "data" / "raw" / "lineage").glob("*")]))


def _mapping_index(mapping: pd.DataFrame) -> dict[tuple[str, str, int], dict[str, Any]]:
    """Build receptor-specific keys; chain D is present in both structures."""
    idx: dict[tuple[str, str, int], dict[str, Any]] = {}
    for _, r in mapping.iterrows():
        status = str(r.get("mapping_status", ""))
        if bool(r.get("qc_required", False)) or not status.startswith("mapped"):
            continue
        cp = str(r.get("common_position", ""))
        if not cp or cp == "nan":
            continue
        for chain_col, num_col, receptor_class in [
            ("alpha1_chain", "alpha1_residue_number", "alpha1"),
            ("alpha2_chain", "alpha2_residue_number", "alpha2"),
        ]:
            chain = str(r.get(chain_col, ""))
            try:
                n = int(float(r.get(num_col)))
            except Exception:
                continue
            if not chain or not chain.startswith(("D", "C", "E")):
                continue
            idx[(receptor_class, chain, n)] = {
                "common_position": cp,
                "mapping_status": status,
                "mapping_method": str(r.get("mapping_method", "")),
                "mapped_residue": str(r.get("alpha1_residue" if receptor_class == "alpha1" else "alpha2_residue", "")),
                "source_chain": chain,
                "source_residue_number": n,
            }
    return idx


def _assign_receptor_mapping(raw: pd.DataFrame, mapping: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    idx = _mapping_index(mapping)
    assigned: list[dict[str, Any]] = []
    qc: list[dict[str, Any]] = []
    for source_row_id, row in raw.iterrows():
        rid = str(row.receptor_id)
        cls = "alpha1" if rid == RECEPTORS["alpha1"] else "alpha2"
        chain = str(row.residue_chain)
        n = int(row.residue_number)
        # C is gamma2 in 6HUP, E is gamma2 in the 9CTJ local construct.
        if cls == "alpha1" and chain == "C":
            key = ("alpha1", "C", n)
        elif cls == "alpha2" and chain == "E":
            key = ("alpha2", "E", n)
        else:
            key = (cls, "D", n)
        m = idx.get(key)
        a = _attrs(row.raw_attributes)
        out = row.to_dict()
        out.update({
            "source_row_id": int(source_row_id),
            "receptor": cls,
            "common_position": m["common_position"] if m else "",
            "mapping_status": m["mapping_status"] if m else "unmapped_receptor_specific",
            "mapping_method": m["mapping_method"] if m else "",
            "mapped_residue": m["mapped_residue"] if m else a.get("restype", ""),
            "mapping_qc_required": bool(m is None),
            "raw_residue_name": a.get("restype", ""),
        })
        assigned.append(out)
        qc.append({
            "source_row_id": int(source_row_id),
            "drug_id": row.drug_id,
            "receptor_id": rid,
            "residue_chain": chain,
            "residue_number": n,
            "common_position": m["common_position"] if m else "",
            "mapping_status": m["mapping_status"] if m else "unmapped_receptor_specific",
            "mapping_method": m["mapping_method"] if m else "",
            "qc_required": bool(m is None),
            "raw_attributes": row.raw_attributes,
        })
    return pd.DataFrame(assigned), pd.DataFrame(qc)


def _distance_values(kind: str, attrs: dict[str, Any]) -> list[float]:
    keys = {
        "hydrophobic_interaction": ["dist"],
        "halogen_bond": ["dist"],
        "hydrogen_bond": ["dist_d-a", "dist_d_a"],
    }.get(kind, [])
    return [_num(attrs[k]) for k in keys if k in attrs and not pd.isna(_num(attrs[k]))]


def _features(assigned: pd.DataFrame, docking: pd.DataFrame) -> pd.DataFrame:
    d = assigned[assigned.common_position.isin(CANDIDATE_POSITIONS)].copy()
    rows: list[dict[str, Any]] = []
    for (receptor, cp, rname, itype), g0 in d.groupby(["receptor", "common_position", "mapped_residue", "interaction_type"], dropna=False):
        for drug in DRUGS:
            g = g0[g0.drug_id == drug]
            rid = RECEPTORS[receptor]
            nposes = int(docking[(docking.receptor_id == rid) & (docking.drug_id == drug)].pose_id.nunique())
            nint = int(g.pose_id.nunique())
            attrs = g.raw_attributes.map(_attrs).tolist() if not g.empty else []
            dist = [v for a in attrs for v in _distance_values(str(itype), a)]
            cent = [_num(a.get("centdist")) for a in attrs if not pd.isna(_num(a.get("centdist")))]
            angle = [_num(a.get("angle")) for a in attrs if not pd.isna(_num(a.get("angle")))]
            offset = [_num(a.get("offset")) for a in attrs if not pd.isna(_num(a.get("offset")))]
            dm, dmed, ds, dn = _stats(dist)
            cm, cmed, cs, cn = _stats(cent)
            am, amed, ass, an = _stats(angle)
            om, omed, os, on = _stats(offset)
            p = len({str(r.pose_id) for _, r in g.iterrows() if _attrs(r.raw_attributes).get("type") == "P"})
            t = len({str(r.pose_id) for _, r in g.iterrows() if _attrs(r.raw_attributes).get("type") == "T"})
            rows.append({
                "drug_id": drug, "receptor": receptor, "receptor_id": rid,
                "common_position": cp, "residue_name": rname,
                "interaction_type": itype, "interaction_count": int(len(g)),
                "n_interacting_poses": nint, "n_poses": nposes,
                "frequency": float(nint / nposes) if nposes else np.nan,
                "data_status": "observed" if len(g) else ("no_interaction_observed" if nposes else "no_data"),
                "distance_mean": dm, "distance_median": dmed, "distance_std": ds, "distance_n_observed": dn,
                "centdist_mean": cm, "centdist_median": cmed, "centdist_std": cs, "centdist_n_observed": cn,
                "angle_mean": am, "angle_median": amed, "angle_std": ass, "angle_n_observed": an,
                "offset_mean": om, "offset_median": omed, "offset_std": os, "offset_n_observed": on,
                "type_P_frequency": float(p / nposes) if itype == "pi_stack" and nposes else np.nan,
                "type_T_frequency": float(t / nposes) if itype == "pi_stack" and nposes else np.nan,
                "type_P_n_interacting_poses": p if itype == "pi_stack" else np.nan,
                "type_T_n_interacting_poses": t if itype == "pi_stack" else np.nan,
            })
    return pd.DataFrame(rows)


def _frequency_matrix(features: pd.DataFrame) -> pd.DataFrame:
    # Keep a row only when the feature was observed for at least one ligand.
    keys = ["receptor", "common_position", "residue_name", "interaction_type"]
    keep = features.groupby(keys).n_interacting_poses.sum()
    keep = set(keep[keep > 0].index)
    x = features[features.set_index(keys).index.isin(keep)].copy()
    rows = []
    for key, g in x.groupby(keys, sort=False):
        row = dict(zip(keys, key))
        for drug in DRUGS:
            z = g[g.drug_id == drug]
            if z.empty:
                row[drug] = 0.0  # zero is absence in lossless raw PLIP with valid poses
                row[f"{drug}_status"] = "no_interaction_observed"
            else:
                row[drug] = float(z.iloc[0].frequency) if pd.notna(z.iloc[0].frequency) else np.nan
                row[f"{drug}_status"] = str(z.iloc[0].data_status)
        row["feature_note"] = "frequency=interacting poses/valid poses; zeros are raw-absence with valid poses"
        rows.append(row)
    return pd.DataFrame(rows)


def _phenotype_raw() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    # Nishino et al. 2008 Table 1: positive/fallen mice on rotarod, 10 mice,
    # oral administration, values transcribed for all reported time points.
    nishino = {
        "diazepam": {2: [1, 3, 4, 3], 5: [2, 4, 5, 4], 10: [4, 7, 9, 8], 20: [9, 10, 10, 9]},
        "triazolam": {0.5: [2, 3, 3, 9], 1: [3, 4, 4, 2], 2: [5, 8, 7, 5], 5: [8, 9, 9, 4]},
    }
    for drug, doses in nishino.items():
        for dose, vals in doses.items():
            for time, val in zip([15, 30, 60, 90], vals):
                rows.append({
                    "drug_id": drug, "species": "mouse", "strain": "not_reported_in_extracted_table",
                    "phenotype_category": "motor_impairment", "assay": "rotarod",
                    "route": "oral", "dose_mg_per_kg": dose, "time_min": time,
                    "endpoint": "positive_or_fallen_mice", "raw_value": val,
                    "raw_unit": "positive_mice_count", "n_tested": 10,
                    "effect_size": val / 10.0, "effect_direction": "motor_impairment",
                    "evidence_tier": "P1", "study_id": "Nishino2008",
                    "reference": "Nishino et al., J Pharmacol Sci 107:349-354 (2008)",
                    "pmid": "18603828", "doi": "10.1254/jphs.08107FP",
                    "source_url": "https://www.jstage.jst.go.jp/article/jphs/107/3/107_08107FP/_pdf",
                    "condition_note": "Raw positive-count transcription; dose and time are not pooled across drugs.",
                })
    # Tanaka et al. 2008: direct rotarod result is qualitative in the accessible
    # text; the righting-reflex values are preserved as a separate endpoint.
    for drug, dose, direction in [("triazolam", 0.3, "motor_coordination_impaired"), ("zolpidem", 3.0, "no_effect_under_test_condition")]:
        rows.append({
            "drug_id": drug, "species": "mouse", "strain": "ddY",
            "phenotype_category": "motor_impairment", "assay": "rotarod",
            "route": "intraperitoneal", "dose_mg_per_kg": dose, "time_min": np.nan,
            "endpoint": "qualitative_motor_coordination_result", "raw_value": np.nan,
            "raw_unit": "qualitative", "n_tested": np.nan, "effect_size": np.nan,
            "effect_direction": direction, "evidence_tier": "P1", "study_id": "Tanaka2008",
            "reference": "Tanaka et al., J Pharmacol Sci 107:277-284 (2008)",
            "pmid": "18603831", "doi": "10.1254/jphs.FP0071991",
            "source_url": "https://www.jstage.jst.go.jp/article/jphs/107/3/107_FP0071991/_pdf",
            "condition_note": "Same-study qualitative rotarod result; no numeric fraction inferred.",
        })
    for drug, dose, mean, sd in [("zolpidem", 3.0, 205.3, 62.6), ("zolpidem", 10.0, 1109.1, 122.1), ("triazolam", 0.3, 234.8, 39.5)]:
        rows.append({
            "drug_id": drug, "species": "mouse", "strain": "ddY",
            "phenotype_category": "sedation_hypnosis", "assay": "righting_reflex_with_thiopental",
            "route": "intraperitoneal_with_thiopental", "dose_mg_per_kg": dose, "time_min": np.nan,
            "endpoint": "righting_reflex_duration", "raw_value": mean, "raw_unit": "seconds_mean",
            "n_tested": np.nan, "effect_size": np.nan, "effect_direction": "sedation_duration",
            "sd": sd, "evidence_tier": "P1", "study_id": "Tanaka2008",
            "reference": "Tanaka et al., J Pharmacol Sci 107:277-284 (2008)",
            "pmid": "18603831", "doi": "10.1254/jphs.FP0071991",
            "source_url": "https://www.jstage.jst.go.jp/article/jphs/107/3/107_FP0071991/_pdf",
            "condition_note": "Supporting endpoint; not merged with rotarod values.",
        })
    # Qualitative supporting evidence is retained without inventing a number.
    support = [
        ("diazepam", "Bourin1992", "rotarod", "single/repeated behavioral models; motor impairment reported", "P2", "https://pubmed.ncbi.nlm.nih.gov/1637802/", "1637802", ""),
        ("alprazolam", "Bourin1992", "rotarod", "myorelaxing/motor effects reported; no extracted numeric endpoint", "P2", "https://pubmed.ncbi.nlm.nih.gov/1637802/", "1637802", ""),
        ("diazepam", "Bayley1996", "rotarod", "dose-related impairment in mice", "P1", "https://pubmed.ncbi.nlm.nih.gov/22302946/", "22302946", "10.1177/026988119601000305"),
        ("zolpidem", "Bayley1996", "rotarod", "dose-related impairment in mice", "P1", "https://pubmed.ncbi.nlm.nih.gov/22302946/", "22302946", "10.1177/026988119601000305"),
        ("alprazolam", "Henauer1984", "rotarod", "rotarod neurologic deficit at 2 mg/kg; tolerance study", "P2", "https://pubmed.ncbi.nlm.nih.gov/6493004/", "6493004", "10.1016/0024-3205(84)90564-2"),
        ("alprazolam", "Lopez1988", "open_field", "low-dose activity increase; higher doses decrease activity", "P2", "https://pubmed.ncbi.nlm.nih.gov/2845448/", "2845448", "10.1016/0091-3057(88)90488-1"),
        ("triazolam", "Lopez1988", "open_field", "activity decreased dose-dependently", "P2", "https://pubmed.ncbi.nlm.nih.gov/2845448/", "2845448", "10.1016/0091-3057(88)90488-1"),
    ]
    for drug, study, assay, direction, tier, url, pmid, doi in support:
        rows.append({
            "drug_id": drug, "species": "mouse", "strain": "not_reported_in_extracted_table",
            "phenotype_category": "motor_impairment" if assay == "rotarod" else "locomotor_activity",
            "assay": assay, "route": "not_reported_in_abstract", "dose_mg_per_kg": np.nan,
            "time_min": np.nan, "endpoint": "qualitative_report", "raw_value": np.nan,
            "raw_unit": "qualitative", "n_tested": np.nan, "effect_size": np.nan,
            "effect_direction": direction, "evidence_tier": tier, "study_id": study,
            "reference": study, "pmid": pmid, "doi": doi, "source_url": url,
            "condition_note": "Qualitative source evidence retained; no numeric value imputed.",
        })
    return pd.DataFrame(rows)


def _primary_phenotype_matrix(raw: pd.DataFrame) -> pd.DataFrame:
    x = raw[(raw.study_id == "Nishino2008") & (raw.assay == "rotarod") & (raw.time_min == 60)].copy()
    rows = []
    for _, r in x.iterrows():
        row = {
            "row_id": f"Nishino2008_{r.drug_id}_{r.dose_mg_per_kg:g}mgkg_60min",
            "phenotype": "motor_impairment",
            "assay": "rotarod",
            "metric": "positive_mice_fraction",
            "study_id": r.study_id,
            "evidence_tier": r.evidence_tier,
            "condition": f"{r.dose_mg_per_kg:g} mg/kg oral, 60 min, n=10",
            "reference": r.reference,
        }
        for drug in DRUGS:
            row[drug] = float(r.effect_size) if drug == r.drug_id else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _literature_validation() -> pd.DataFrame:
    rows = [
        ("gamma2", "TYR58", "BZD_GAMMA2_058", "consistent", "Direct proximity/mutagenesis evidence for gamma2 Tyr58 in classical BZD pocket", "9351978", "", "https://pubmed.ncbi.nlm.nih.gov/9351978/"),
        ("gamma2", "PHE77", "BZD_GAMMA2_077", "consistent", "gamma2 Phe77 is required for high-affinity binding/modulation of several BZD-site ligands", "9351978", "", "https://pubmed.ncbi.nlm.nih.gov/9351978/"),
        ("gamma2", "ASN60", "BZD_GAMMA2_060", "consistent", "Cysteine-reactive ligand study located gamma2 Asn60 relative to classical BZDs", "", "10.1021/cb500186a", "https://pubs.acs.org/doi/10.1021/cb500186a"),
        ("alpha1/alpha2", "HIS102/H101", "BZD_SITE_102", "consistent", "Conserved alpha histidine is required for diazepam/BZD agonist sensitivity; numbering differs by construct", "1346133", "", "https://pubmed.ncbi.nlm.nih.gov/1346133/"),
        ("alpha1", "PHE100", "BZD_SITE_100", "partially_consistent", "Loop-A neighboring Phe100 is close to the BZD pocket; direct ligand proximity is reported for adjacent alpha residues", "17854801", "", "https://pubmed.ncbi.nlm.nih.gov/17854801/"),
        ("alpha1", "LYS156/157", "BZD_SITE_156/BZD_SITE_157", "not_established", "No direct primary residue-specific validation was identified in the searched set", "", "", ""),
        ("alpha1", "SER205/206", "BZD_SITE_205/BZD_SITE_206", "partially_consistent", "Alpha loop-C positions 206/209 influence BZD-site ligand affinity; exact numbering differs", "9380031", "", "https://pubmed.ncbi.nlm.nih.gov/9380031/"),
        ("alpha1", "TYR210/211", "BZD_SITE_210/BZD_SITE_211", "consistent", "Alpha tyrosines near 159/209-210 are reported as critical for diazepam binding/modulation and aromatic pocket structure", "9145922", "10.1124/mol.51.5.833", "https://pubmed.ncbi.nlm.nih.gov/9145922/"),
        ("alpha1/alpha2", "BZD interface", "all", "partially_consistent", "Structural review places F100/H102/Y160/Y210 and gamma2 Y58/F77 at the alpha-gamma BZD interface", "", "10.1038/s41586-018-0255-3", "https://www.nature.com/articles/s41586-018-0255-3"),
    ]
    return pd.DataFrame(rows, columns=["subunit_scope", "residue_label", "common_position", "validation_class", "evidence_summary", "pmid", "doi", "source_url"])


def _heatmap(matrix: pd.DataFrame, out: Path, title: str, value_cols: list[str]) -> None:
    if matrix.empty:
        return
    vals = matrix[value_cols].astype(float).to_numpy()
    fig_h = max(4.5, 0.32 * len(matrix) + 1.5)
    fig, ax = plt.subplots(figsize=(8.6, fig_h))
    im = ax.imshow(vals, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(value_cols)), value_cols, rotation=30, ha="right")
    labels = [f"{r.receptor}:{r.common_position}:{r.residue_name}:{r.interaction_type}" for _, r in matrix.iterrows()]
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    ax.set_title(title)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            if np.isfinite(vals[i, j]):
                ax.text(j, i, f"{vals[i,j]:.2f}", ha="center", va="center", fontsize=6, color="white" if vals[i,j] < .55 else "black")
    fig.colorbar(im, ax=ax, label="interacting pose frequency")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)


def _phenotype_heatmap(matrix: pd.DataFrame, out: Path) -> None:
    if matrix.empty:
        return
    vals = matrix[DRUGS].astype(float).to_numpy()
    fig, ax = plt.subplots(figsize=(7.5, max(3.5, 0.38 * len(matrix) + 1.3)))
    cmap = plt.cm.magma.copy(); cmap.set_bad("#d9d9d9")
    im = ax.imshow(np.ma.masked_invalid(vals), aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(len(DRUGS)), DRUGS, rotation=30, ha="right")
    ax.set_yticks(range(len(matrix)), matrix.row_id, fontsize=7)
    ax.set_title("Mouse rotarod primary evidence (60 min; study-specific dose)")
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            if np.isfinite(vals[i, j]):
                ax.text(j, i, f"{vals[i,j]:.2f}", ha="center", va="center", fontsize=7, color="white")
    fig.colorbar(im, ax=ax, label="positive/fallen mice fraction")
    fig.tight_layout(); fig.savefig(out, dpi=220); plt.close(fig)


def _integrated_figure(struct: pd.DataFrame, phen: pd.DataFrame, out: Path) -> None:
    fig = plt.figure(figsize=(14, max(5.0, 0.35 * max(len(struct), len(phen)) + 2)))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1])
    ax1 = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[0, 1])
    if not struct.empty:
        vals = struct[DRUGS].astype(float).to_numpy()
        im1 = ax1.imshow(vals, aspect="auto", cmap="viridis", vmin=0, vmax=1)
        labels = [f"{r.receptor}:{r.common_position}:{r.residue_name}:{r.interaction_type}" for _, r in struct.iterrows()]
        ax1.set_yticks(range(len(labels)), labels, fontsize=6); ax1.set_xticks(range(4), DRUGS, rotation=30, ha="right")
        ax1.set_title("Structural interaction frequency")
        fig.colorbar(im1, ax=ax1, fraction=.03, pad=.02)
    if not phen.empty:
        vals = phen[DRUGS].astype(float).to_numpy()
        cmap = plt.cm.magma.copy(); cmap.set_bad("#d9d9d9")
        im2 = ax2.imshow(np.ma.masked_invalid(vals), aspect="auto", cmap=cmap, vmin=0, vmax=1)
        ax2.set_yticks(range(len(phen)), phen.row_id, fontsize=6); ax2.set_xticks(range(4), DRUGS, rotation=30, ha="right")
        ax2.set_title("Mouse rotarod evidence (numeric cells only)")
        fig.colorbar(im2, ax=ax2, fraction=.05, pad=.04)
    fig.suptitle("Structure-to-mouse phenotype bridge: descriptive evidence")
    fig.tight_layout(); fig.savefig(out, dpi=220); plt.close(fig)


def _report(run_dir: Path, manifest: dict[str, Any], mapping_qc: pd.DataFrame, features: pd.DataFrame, matrix: pd.DataFrame, aromatic: pd.DataFrame, phen_raw: pd.DataFrame, phen: pd.DataFrame, literature: pd.DataFrame) -> None:
    completed = run_manager.iso_now()
    header = [
        f"Run ID: {manifest['run_id']}",
        f"Parent Run: {manifest.get('parent_run') or ''}",
        f"Analysis Phase: {manifest['analysis_phase']}",
        f"Started: {manifest['started_at']}",
        f"Completed: {completed}", "",
    ]
    # High-difference features are descriptive ranges within a receptor.
    high = []
    for key, g in matrix.groupby(["receptor", "common_position", "residue_name", "interaction_type"]):
        vals = pd.to_numeric(g[DRUGS].iloc[0], errors="coerce")
        high.append((*key, float(vals.max() - vals.min())))
    high_df = pd.DataFrame(high, columns=["receptor", "common_position", "residue_name", "interaction_type", "range"])
    high_text = high_df.sort_values("range", ascending=False).head(12).to_string(index=False) if not high_df.empty else "No observed candidate features"
    aromatic_text = aromatic.to_string(index=False) if not aromatic.empty else "No observed pi_stack rows"
    literature_text = literature.to_string(index=False)
    report = header + [
        "# STRUCTURE_TO_MOUSE_PHENOTYPE_REPORT", "",
        "This run links residue-level interaction frequencies to a lossless, literature-derived mouse phenotype evidence table. It is descriptive and hypothesis-generating. It does not claim that docking features cause in vivo behavior, and it does not fit a predictive or statistical model.", "",
        "## Q1. Which structural drug differences are present?", "",
        "The structural matrix contains only candidate residue/interaction rows observed in the corrected raw PLIP table. Hydrophobic interaction, pi_stack, hydrogen_bond, and halogen_bond are kept as separate features. Frequency is interacting poses divided by valid poses; a zero in the matrix means no raw interaction row for that drug-feature while valid poses existed. The largest within-receptor frequency ranges are:", "",
        high_text, "",
        "## Q2. TYR58/PHE77 pi-stacking", "",
        "TYR58 and PHE77 are reported separately for each receptor and drug. Missing geometry is not replaced by zero. The full table below retains centdist, angle, offset and PLIP type P/T frequencies:", "",
        aromatic_text, "",
        "## Q3. Primary mouse phenotype", "",
        "Acute motor impairment measured by rotarod was selected because all four drugs have at least some mouse motor/rotarod evidence in the searched literature. Numeric coverage is incomplete: Nishino et al. provide raw positive/fallen-mouse counts for diazepam and triazolam; Tanaka et al. provide qualitative rotarod results for triazolam and zolpidem; alprazolam is supported by separate mouse studies. The numeric heatmap therefore shows only exact Nishino values and leaves the other drug cells as NA. Doses are retained as study conditions and are not converted into a cross-drug potency ranking.", "",
        "## Q4. Is there a same-condition four-drug numeric comparison?", "",
        "No. The reviewed evidence does not provide a defensible single mouse assay/strain/route/time/dose design with numeric rotarod values for all four drugs. Accordingly, no missing values were imputed and no mg/kg-only ranking was made.", "",
        "## Q5. Literature validation of high-difference residues", "",
        literature_text, "",
        "The classifications distinguish direct primary residue evidence from interface-level or adjacent-numbering evidence. They do not establish that a particular docking interaction caused a behavioral effect.", "",
        "## Q6. Integrated interpretation and limits", "",
        "The output supports a residue-level hypothesis test: within a fixed receptor construct, drugs differ in which BZD-site contacts recur across poses and in the geometry of observed pi-stacking. The mouse table provides separately sourced phenotype context, with study, route, assay, endpoint, and raw values preserved. Because the structural and animal records are not matched experimental conditions, the integrated figure is an evidence map rather than a correlation or prediction. No alpha1-vs-alpha2 global similarity or ML analysis is used as the primary result.", "",
        "## QC", "",
        f"Raw PLIP rows re-assigned: {len(mapping_qc)}; receptor-specific mapping unresolved: {int(mapping_qc.qc_required.sum())}; candidate feature rows (drug-level): {len(features)}; frequency-matrix rows: {len(matrix)}; raw phenotype evidence rows: {len(phen_raw)}; numeric primary phenotype rows: {len(phen)}.",
        "The receptor-specific key prevents chain D in alpha1 and alpha2 structures from overwriting each other. Existing corrected mapping and raw attributes are copied into data/raw/lineage. No previous run was modified.", "",
        "## Files", "",
        "results/figures/residue_ligand_frequency_heatmap.png; results/tables/residue_ligand_frequency_matrix.csv; results/tables/aromatic_geometry_by_ligand.csv; data/processed/mouse_phenotype_raw.csv; results/tables/mouse_primary_phenotype_matrix.csv; results/figures/mouse_primary_phenotype_heatmap.png; results/figures/structural_vs_mouse_phenotype.png; results/tables/residue_literature_validation.csv",
    ]
    text = "\n".join(report) + "\n"
    (run_dir / "reports" / "STRUCTURE_TO_MOUSE_PHENOTYPE_REPORT.md").write_text(text)
    (run_dir / "FINAL_REPORT.md").write_text(text)
    qc = "\n".join(header + ["# QC_REPORT", "", f"Receptor-specific mapping QC unresolved: {int(mapping_qc.qc_required.sum())} / {len(mapping_qc)}.", "No unobserved interaction feature was added. Missing numeric phenotype values remain NA.", ""])
    (run_dir / "QC_REPORT.md").write_text(qc)


def run() -> dict[str, Any]:
    source = _source_dir()
    run_dir, manifest = run_manager.create_run(
        short_phase_name="structure_mouse_bridge",
        analysis_phase="structure-to-mouse phenotype evidence bridge",
        analysis_purpose="Re-map corrected BZD-interface PLIP contacts by receptor-specific chain keys and place them beside raw mouse phenotype evidence.",
        parent_run=SOURCE_RUN,
        drugs=DRUGS,
        receptors=RECEPTORS,
        docking_parameters={"docking_executed": False, "source_run": SOURCE_RUN, "frequency_definition": "interacting poses / valid poses"},
        input_files=[f"runs/{SOURCE_RUN}/raw/plip_interactions_4drug.csv", f"runs/{SOURCE_RUN}/raw/docking_results_4drug.csv", f"runs/{SOURCE_RUN}/tables/corrected_residue_mapping.csv"],
        notes=["New immutable run; no prior run modified.", "Receptor-specific mapping key prevents chain-D collision between alpha1 and alpha2 contexts.", "Mouse phenotype records are raw literature evidence; no imputation or cross-study pooling."],
        source_run=SOURCE_RUN,
        source_files=["plip_interactions_4drug.csv", "docking_results_4drug.csv", "corrected_residue_mapping.csv", "recovered_interactions_qc.csv"],
        code_paths=[Path(__file__), ROOT / "src" / "run_manager.py"],
    )
    _copy_lineage(run_dir, source, manifest)
    run_manager.write_manifest(run_dir, manifest)
    raw = pd.read_csv(source / "raw" / "plip_interactions_4drug.csv")
    docking = pd.read_csv(source / "raw" / "docking_results_4drug.csv")
    mapping = pd.read_csv(source / "tables" / "corrected_residue_mapping.csv")
    mapping.to_csv(run_dir / "results" / "tables" / "corrected_residue_mapping.csv", index=False)
    assigned, mapping_qc = _assign_receptor_mapping(raw, mapping)
    assigned.to_csv(run_dir / "data" / "processed" / "reassigned_interactions.csv", index=False)
    mapping_qc.to_csv(run_dir / "results" / "tables" / "mapping_reassignment_qc.csv", index=False)
    features = _features(assigned, docking)
    features.to_csv(run_dir / "results" / "tables" / "candidate_residue_features.csv", index=False)
    matrix = _frequency_matrix(features)
    matrix.to_csv(run_dir / "results" / "tables" / "residue_ligand_frequency_matrix.csv", index=False)
    aromatic = features[(features.interaction_type == "pi_stack") & features.common_position.isin({"BZD_GAMMA2_058", "BZD_GAMMA2_077", "BZD_SITE_100", "BZD_SITE_101", "BZD_SITE_102", "BZD_SITE_156", "BZD_SITE_157", "BZD_SITE_205", "BZD_SITE_206", "BZD_SITE_210", "BZD_SITE_211"})].copy()
    aromatic.to_csv(run_dir / "results" / "tables" / "aromatic_geometry_by_ligand.csv", index=False)
    phen = _phenotype_raw()
    phen.to_csv(run_dir / "data" / "processed" / "mouse_phenotype_raw.csv", index=False)
    phen_matrix = _primary_phenotype_matrix(phen)
    phen_matrix.to_csv(run_dir / "results" / "tables" / "mouse_primary_phenotype_matrix.csv", index=False)
    literature = _literature_validation()
    literature.to_csv(run_dir / "results" / "tables" / "residue_literature_validation.csv", index=False)
    _heatmap(matrix, run_dir / "results" / "figures" / "residue_ligand_frequency_heatmap.png", "BZD-site residue × ligand interaction frequency", DRUGS)
    _phenotype_heatmap(phen_matrix, run_dir / "results" / "figures" / "mouse_primary_phenotype_heatmap.png")
    _integrated_figure(matrix, phen_matrix, run_dir / "results" / "figures" / "structural_vs_mouse_phenotype.png")
    _report(run_dir, manifest, mapping_qc, features, matrix, aromatic, phen, phen_matrix, literature)
    finalized = run_manager.finalize_run(run_dir, manifest, qc_status="COMPLETED_WITH_QC", extra_notes=["Structural matrix includes only candidate features observed in lossless raw PLIP; absent drug-feature rows are labeled no_interaction_observed."])
    return {"run_id": finalized["run_id"], "run_dir": str(run_dir), "mapping_qc_rows": len(mapping_qc), "mapping_unresolved": int(mapping_qc.qc_required.sum()), "feature_rows": len(features), "frequency_matrix_rows": len(matrix), "phenotype_rows": len(phen), "numeric_phenotype_rows": len(phen_matrix)}


def main() -> None:
    print(json.dumps(run(), ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
