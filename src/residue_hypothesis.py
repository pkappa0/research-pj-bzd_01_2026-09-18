"""Residue-level recovery of BZD-interface interactions.

This phase is deliberately separate from the multi-receptor feature matrix.  It
extends the existing alpha1-chain-D/alpha2-chain-D mapping with the gamma2
chain C/E mapping needed for the BZD interface in the selected structures.
Hydrophobic contacts and PLIP pi stacking remain separate interaction types.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping

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
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}


def _safe_int(value: object) -> int | None:
    try:
        if pd.isna(value):
            return None
        return int(float(value))
    except Exception:
        return None


def _attrs(raw: object) -> dict[str, Any]:
    try:
        value = json.loads(str(raw))
        return value if isinstance(value, dict) else {}
    except Exception:
        return {}


def _num(value: object) -> float:
    try:
        value = float(value)
        return value if math.isfinite(value) else np.nan
    except Exception:
        return np.nan


def _stats(values: Iterable[object]) -> tuple[float, float, float, int]:
    x = pd.to_numeric(pd.Series(list(values)), errors="coerce").dropna()
    if x.empty:
        return np.nan, np.nan, np.nan, 0
    return float(x.mean()), float(x.median()), float(x.std(ddof=1)) if len(x) > 1 else 0.0, int(len(x))


def _chain_records(path: Path, chain_id: str) -> list[Any]:
    from Bio.PDB import PDBParser

    structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))[0][chain_id]
    return [res for res in structure if res.id[0] == " " and "CA" in res]


def _aligned_chain_mapping(
    source_records: list[Any], target_records: list[Any], *, prefix: str, source_chain: str, target_chain: str,
    start: int = 50, end: int = 220,
) -> tuple[pd.DataFrame, dict[tuple[str, int], dict[str, Any]], dict[str, Any]]:
    from Bio.Align import PairwiseAligner
    from Bio.SeqUtils import seq1

    source_seq = "".join(seq1(r.resname, custom_map={"MSE": "M"}) for r in source_records)
    target_seq = "".join(seq1(r.resname, custom_map={"MSE": "M"}) for r in target_records)
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(source_seq, target_seq)
    alignment = alignments[0]
    optimal_count = len(alignments)
    source_to_target: dict[int, int] = {}
    for block_a, block_b in zip(*alignment.aligned):
        for ia, ib in zip(range(block_a[0], block_a[1]), range(block_b[0], block_b[1])):
            source_to_target[ia] = ib
    source_by_num = {int(r.id[1]): (i, r) for i, r in enumerate(source_records)}
    rows: list[dict[str, Any]] = []
    index: dict[tuple[str, int], dict[str, Any]] = {}
    for number, (source_index, source_res) in source_by_num.items():
        if not (start <= number <= end):
            continue
        target_index = source_to_target.get(source_index)
        target_res = target_records[target_index] if target_index is not None and target_index < len(target_records) else None
        status = "mapped_gamma2_sequence_alignment" if target_res is not None and optimal_count == 1 else "unavailable_or_alignment_ambiguous"
        common = f"{prefix}_{number:03d}"
        row = {
            "common_position": common,
            "interface_subunit": "gamma2",
            "alpha1_interface_subunit": "gamma2",
            "alpha2_interface_subunit": "gamma2",
            "alpha1_chain": source_chain,
            "alpha1_residue_number": number,
            "alpha1_residue": source_res.resname,
            "alpha2_chain": target_chain,
            "alpha2_residue_number": int(target_res.id[1]) if target_res is not None else "",
            "alpha2_residue": target_res.resname if target_res is not None else "",
            "mapping_method": "global sequence alignment (Biopython PairwiseAligner) + PDB residue lookup",
            "alignment_score": float(alignment.score),
            "alignment_n_optimal": int(optimal_count),
            "residue_identity": bool(target_res is not None and source_res.resname == target_res.resname),
            "mapping_status": status,
            "qc_required": status != "mapped_gamma2_sequence_alignment",
            "mapping_scope": "BZD_interface_gamma2_C_to_E",
        }
        rows.append(row)
        index[(source_chain, number)] = row
        if target_res is not None:
            index[(target_chain, int(target_res.id[1]))] = row
    meta = {
        "alignment_score": float(alignment.score),
        "alignment_n_optimal": int(optimal_count),
        "source_residue_count": len(source_records),
        "target_residue_count": len(target_records),
    }
    return pd.DataFrame(rows), index, meta


def _prior_mapping(path: Path) -> tuple[pd.DataFrame, dict[tuple[str, int], dict[str, Any]]]:
    prior = pd.read_csv(path)
    index: dict[tuple[str, int], dict[str, Any]] = {}
    for _, row in prior.iterrows():
        n1 = _safe_int(row.get("alpha1_residue_number")); n2 = _safe_int(row.get("alpha2_residue_number"))
        if n1 is not None and str(row.get("alpha1_chain", "")):
            value = row.to_dict()
            value["interface_subunit"] = "alpha1"
            value["alpha1_interface_subunit"] = "alpha1"
            value["alpha2_interface_subunit"] = "alpha2"
            index[(str(row.alpha1_chain), n1)] = value
        if n2 is not None and str(row.get("alpha2_chain", "")):
            value = row.to_dict()
            value["interface_subunit"] = "alpha2"
            value["alpha1_interface_subunit"] = "alpha1"
            value["alpha2_interface_subunit"] = "alpha2"
            index[(str(row.alpha2_chain), n2)] = value
    return prior, index


def _corrected_mapping(run_dir: Path, source_run: Path) -> tuple[pd.DataFrame, dict[tuple[str, int], dict[str, Any]], dict[str, Any]]:
    prior, prior_index = _prior_mapping(source_run / "raw" / "residue_mapping.csv")
    gamma1 = _chain_records(ROOT / "data/raw/structures/strict3/6HUP.pdb", "C")
    gamma2 = _chain_records(ROOT / "data/raw/structures/strict3/9CTJ.pdb", "E")
    gamma_df, gamma_index, gamma_meta = _aligned_chain_mapping(gamma1, gamma2, prefix="BZD_GAMMA2", source_chain="C", target_chain="E")

    prior_out = prior.copy()
    prior_out["interface_subunit"] = "alpha1/alpha2"
    prior_out["alpha1_interface_subunit"] = "alpha1"
    prior_out["alpha2_interface_subunit"] = "alpha2"
    prior_out["mapping_scope"] = "BZD_interface_alpha1_D_to_alpha2_D_prior"
    prior_out["alignment_n_optimal"] = np.nan
    prior_out["residue_identity"] = prior_out.apply(lambda r: str(r.get("alpha1_residue", "")) == str(r.get("alpha2_residue", "")), axis=1)
    corrected = pd.concat([prior_out, gamma_df], ignore_index=True, sort=False)
    corrected = corrected.drop_duplicates(subset=["interface_subunit", "alpha1_chain", "alpha1_residue_number", "alpha2_chain", "alpha2_residue_number"], keep="first")
    corrected = corrected.sort_values(["interface_subunit", "alpha1_residue_number"], na_position="last").reset_index(drop=True)
    corrected.to_csv(run_dir / "tables" / "corrected_residue_mapping.csv", index=False)

    mapping_index = dict(prior_index)
    mapping_index.update(gamma_index)
    return corrected, mapping_index, gamma_meta


def _recover_interactions(run_dir: Path, source_run: Path, mapping_index: dict[tuple[str, int], dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(source_run / "raw" / "plip_interactions_4drug.csv")
    records: list[dict[str, Any]] = []
    recovered: list[dict[str, Any]] = []
    old_unresolved_path = source_run / "tables" / "unresolved_feature_qc.csv"
    if old_unresolved_path.exists():
        old_unresolved = pd.read_csv(old_unresolved_path)
        old_unresolved_keys = set(zip(old_unresolved.receptor_id.astype(str), old_unresolved.pose_id.astype(str), old_unresolved.interaction_type.astype(str), old_unresolved.residue_chain.astype(str), old_unresolved.residue_number.astype(str)))
    else:
        # A later residue-hypothesis run may contain only the lossless raw
        # source files. Chain C/E is still explicitly known to have been
        # unresolved in the original feature-engineering run.
        old_unresolved_keys = set()
    for source_row_id, row in raw.iterrows():
        attrs = _attrs(row.raw_attributes)
        number = _safe_int(row.residue_number)
        key = (str(row.residue_chain), number) if number is not None else None
        mapping = mapping_index.get(key) if key else None
        prior_key = (str(row.receptor_id), str(row.pose_id), str(row.interaction_type), str(row.residue_chain), str(row.residue_number))
        was_unresolved = prior_key in old_unresolved_keys or str(row.residue_chain) in {"C", "E"}
        if mapping and str(mapping.get("mapping_status", "")).startswith("mapped") and not bool(mapping.get("qc_required", False)):
            recovery_status = "recovered_bzd_interface" if str(row.residue_chain) in {"C", "E"} else "already_mapped_prior"
            recovery_reason = "gamma2_C_to_E_alignment" if recovery_status == "recovered_bzd_interface" else "prior_alpha_chain_mapping"
            recovered_row = {
                "source_row_id": int(source_row_id),
                **row.to_dict(),
                "common_position": mapping.get("common_position", ""),
                "interface_subunit": mapping.get("interface_subunit", "alpha1"),
                "residue_name": attrs.get("restype", ""),
                "mapping_status": mapping.get("mapping_status", ""),
                "mapping_method": mapping.get("mapping_method", ""),
                "qc_required": bool(mapping.get("qc_required", False)),
                "recovery_status": recovery_status,
                "recovery_reason": recovery_reason,
                "originally_unresolved": bool(was_unresolved),
            }
            recovered.append(recovered_row)
        records.append({
            "source_row_id": int(source_row_id),
            **row.to_dict(),
            "residue_name": attrs.get("restype", ""),
            "common_position": mapping.get("common_position", "") if mapping else "",
            "interface_subunit": mapping.get("interface_subunit", "") if mapping else "",
            "mapping_status": mapping.get("mapping_status", "unmapped") if mapping else "unmapped",
            "mapping_method": mapping.get("mapping_method", "") if mapping else "",
            "qc_required": bool(mapping.get("qc_required", True)) if mapping else True,
            "original_mapping_status": "unresolved" if was_unresolved else "mapped_prior_or_not_in_prior_qc",
            "originally_unresolved": bool(was_unresolved),
            "recovery_status": recovery_status if mapping and str(mapping.get("mapping_status", "")).startswith("mapped") and not bool(mapping.get("qc_required", False)) else "unresolved_after_correction",
            "recovery_reason": recovery_reason if mapping and str(mapping.get("mapping_status", "")).startswith("mapped") and not bool(mapping.get("qc_required", False)) else "no_unambiguous_chain_residue_mapping",
        })
    qc = pd.DataFrame(records)
    recovered_df = pd.DataFrame(recovered)
    qc.to_csv(run_dir / "tables" / "recovered_interactions_qc.csv", index=False)
    return recovered_df, qc


def _distance_values(interaction_type: str, attrs: Mapping[str, Any]) -> list[float]:
    keys = {
        "hydrophobic_interaction": ["dist"],
        "halogen_bond": ["dist"],
        "hydrogen_bond": ["dist-d-a", "dist_d-a"],
        "water_bridge": ["dist-a-w", "dist_a_w", "dist"],
        "salt_bridge": ["dist"],
    }.get(interaction_type, [])
    return [_num(attrs[k]) for k in keys if k in attrs and not pd.isna(_num(attrs[k]))]


def _feature_rows(recovered: pd.DataFrame, docking: pd.DataFrame, receptor_name: str) -> pd.DataFrame:
    receptor_id = RECEPTORS[receptor_name]
    d = recovered[recovered.receptor_id == receptor_id].copy()
    all_features = d[["common_position", "interface_subunit", "residue_name", "interaction_type"]].drop_duplicates().to_dict("records") if not d.empty else []
    # Keep the two priority aromatic interaction slots explicit even when a
    # receptor has no observed pi_stack row for that residue.  The resulting
    # row is marked no_interaction_observed rather than being dropped.
    existing = {(x["common_position"], x["interaction_type"]) for x in all_features}
    for position, residue in (("BZD_GAMMA2_058", "TYR"), ("BZD_GAMMA2_077", "PHE")):
        if (position, "pi_stack") not in existing:
            all_features.append({"common_position": position, "interface_subunit": "gamma2", "residue_name": residue, "interaction_type": "pi_stack"})
    pose_counts = docking[docking.receptor_id == receptor_id].groupby("drug_id").pose_id.nunique().to_dict()
    rows: list[dict[str, Any]] = []
    for drug in DRUGS:
        n_poses = int(pose_counts.get(drug, 0))
        for feature in all_features:
            g = d[(d.drug_id == drug) & (d.common_position == feature["common_position"]) & (d.interaction_type == feature["interaction_type"])]
            attrs = g.raw_attributes.map(_attrs) if not g.empty else pd.Series(dtype=object)
            unique_poses = int(g.pose_id.nunique())
            distances = [v for a in attrs for v in _distance_values(feature["interaction_type"], a)]
            dm, dmed, dstd, dn = _stats(distances)
            cent = [_num(a.get("centdist")) for a in attrs if not pd.isna(_num(a.get("centdist")))]
            angle = [_num(a.get("angle")) for a in attrs if not pd.isna(_num(a.get("angle")))]
            offset = [_num(a.get("offset")) for a in attrs if not pd.isna(_num(a.get("offset")))]
            cm, cmed, cstd, cn = _stats(cent)
            am, amed, astd, an = _stats(angle)
            om, omed, ostd, on = _stats(offset)
            p_poses = int(g[["pose_id", "raw_attributes"]].assign(_type=g.raw_attributes.map(lambda x: _attrs(x).get("type", ""))).query("_type == 'P'").pose_id.nunique()) if not g.empty else 0
            t_poses = int(g[["pose_id", "raw_attributes"]].assign(_type=g.raw_attributes.map(lambda x: _attrs(x).get("type", ""))).query("_type == 'T'").pose_id.nunique()) if not g.empty else 0
            rows.append({
                "drug_id": drug,
                "receptor": receptor_name,
                "receptor_id": receptor_id,
                "common_position": feature["common_position"],
                "interface_subunit": feature["interface_subunit"],
                "residue_name": feature["residue_name"],
                "interaction_type": feature["interaction_type"],
                "interaction_count": int(len(g)),
                "n_interacting_poses": unique_poses,
                "n_poses": n_poses,
                "frequency": float(unique_poses / n_poses) if n_poses else np.nan,
                "data_status": "observed" if len(g) else ("no_interaction_observed" if n_poses else "no_data"),
                "distance_mean": dm, "distance_median": dmed, "distance_std": dstd, "distance_n_observed": dn,
                "centdist_mean": cm, "centdist_median": cmed, "centdist_std": cstd, "centdist_n_observed": cn,
                "angle_mean": am, "angle_median": amed, "angle_std": astd, "angle_n_observed": an,
                "offset_mean": om, "offset_median": omed, "offset_std": ostd, "offset_n_observed": on,
                "type_P_frequency": float(p_poses / n_poses) if feature["interaction_type"] == "pi_stack" and n_poses else np.nan,
                "type_T_frequency": float(t_poses / n_poses) if feature["interaction_type"] == "pi_stack" and n_poses else np.nan,
                "type_P_n_interacting_poses": p_poses if feature["interaction_type"] == "pi_stack" else np.nan,
                "type_T_n_interacting_poses": t_poses if feature["interaction_type"] == "pi_stack" else np.nan,
            })
    return pd.DataFrame(rows)


def _comparison(features: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    group_cols = ["receptor", "receptor_id", "common_position", "interface_subunit", "residue_name", "interaction_type"]
    for key, g in features.groupby(group_cols, dropna=False):
        g = g.set_index("drug_id").reindex(DRUGS).reset_index()
        observed = g[g.data_status == "observed"]
        row = dict(zip(group_cols, key))
        row.update({
            "n_drugs_with_interaction": int(len(observed)),
            "frequency_range": float(observed.frequency.max() - observed.frequency.min()) if not observed.empty else np.nan,
            "frequency_high_drug": str(observed.loc[observed.frequency.idxmax(), "drug_id"]) if not observed.empty else "",
            "frequency_low_drug": str(observed.loc[observed.frequency.idxmin(), "drug_id"]) if not observed.empty else "",
            "frequency_values": ";".join(f"{r.drug_id}={r.frequency:g}" if pd.notna(r.frequency) else f"{r.drug_id}=NA" for _, r in g.iterrows()),
            "distance_mean_range": _range(g.distance_mean),
            "centdist_mean_range": _range(g.centdist_mean),
            "angle_mean_range": _range(g.angle_mean),
            "offset_mean_range": _range(g.offset_mean),
            "profile_note": "frequency=interacting poses / available poses; geometry statistics use raw PLIP rows; hydrophobic and pi_stack remain separate",
        })
        rows.append(row)
    return pd.DataFrame(rows)


def _range(values: Iterable[object]) -> float:
    x = pd.to_numeric(pd.Series(list(values)), errors="coerce").dropna()
    return float(x.max() - x.min()) if not x.empty else np.nan


def _report(run_dir: Path, manifest: Mapping[str, Any], corrected: pd.DataFrame, recovered: pd.DataFrame, qc: pd.DataFrame, features: dict[str, pd.DataFrame], comparisons: pd.DataFrame, aromatic: pd.DataFrame, gamma_meta: Mapping[str, Any]) -> None:
    completed = run_manager.iso_now()
    header = [
        f"Run ID: {manifest['run_id']}",
        f"Parent Run: {manifest.get('parent_run') or ''}",
        f"Analysis Phase: {manifest['analysis_phase']}",
        f"Started: {manifest['started_at']}",
        f"Completed: {completed}",
        "",
    ]
    recovered_n = int((qc.recovery_status == "recovered_bzd_interface").sum())
    unresolved_n = int((qc.recovery_status == "unresolved_after_correction").sum())
    site = aromatic[aromatic.priority_site]
    site_text = site[["receptor", "drug_id", "site_label", "interaction_type", "frequency", "distance_mean", "centdist_mean", "angle_mean", "offset_mean", "type_P_frequency", "type_T_frequency"]].to_string(index=False) if not site.empty else "No TYR58/PHE77 rows"
    comp_site = comparisons[(comparisons.common_position.isin(["BZD_GAMMA2_058", "BZD_GAMMA2_077"])) & (comparisons.interaction_type == "pi_stack")]
    comp_text = comp_site.to_string(index=False) if not comp_site.empty else "No site-level pi_stack comparison rows"
    report = header + [
        "# RESIDUE_HYPOTHESIS_REPORT", "",
        "This run tests residue-level, within-receptor interaction hypotheses. It does not perform global alpha1-vs-alpha2 similarity analysis, ML, imputation, activity integration, or docking.", "",
        "## Recovery and mapping", "",
        f"Raw PLIP rows: **{len(qc)}**; previously unresolved rows: **{int(qc.originally_unresolved.sum())}**; recovered through gamma2 C→E mapping: **{recovered_n}**; unresolved after correction: **{unresolved_n}**.",
        f"Gamma2 alignment score: **{gamma_meta['alignment_score']}**; optimal alignments: **{gamma_meta['alignment_n_optimal']}**. The selected alignment is unique in the Biopython PairwiseAligner result.",
        "The prior alpha1 chain D→alpha2 chain D mapping is retained. The new BZD-interface mapping is gamma2 chain C in 6HUP to gamma2 chain E in the provisional 9CTJ local construct. Mapping is sequence/structure based; no PDB-number-only identity was assumed.", "",
        "## Feature policy", "",
        "`hydrophobic_interaction` and `pi_stack` are separate features. Pi stacking retains centdist, angle, offset, and PLIP type P/T frequencies; `frequency` and type P/T frequencies are unique interacting poses divided by available poses, while geometry statistics use raw PLIP interaction rows. Hydrophobic contacts retain PLIP distance statistics. Hydrogen bonds, halogen bonds, and other recovered interaction types remain in the same residue-level tables.", "",
        "## TYR58/PHE77 aromatic interaction summary", "",
        site_text, "",
        "The priority-site rows above are descriptive within each receptor. A missing row for a drug is represented as `no_interaction_observed` when pose data exist; it is not treated as missing data or as a zero geometry.", "",
        "## Drug-level residue comparison", "",
        comp_text, "",
        "Frequency ranges compare the four drugs on the same receptor and feature. Geometry ranges are descriptive and are not used to infer a statistical effect.", "",
        "## QC boundaries", "",
        f"Feature rows: alpha1={len(features['alpha1'])}, alpha2={len(features['alpha2'])}; corrected mapping rows: {len(corrected)}. Recovered interaction rows are retained in tables/recovered_interactions_qc.csv with source_row_id and raw_attributes.",
        "The alpha2 receptor remains a provisional local construct derived from 9CTJ. Therefore alpha1-vs-alpha2 comparisons are secondary context only; the primary question here is drug-to-drug variation within alpha1 and within alpha2.",
        "No interaction was collapsed across interaction types, residues, drugs, poses, or geometry fields. No values were imputed.", "",
        "## Files", "",
        "tables/corrected_residue_mapping.csv; tables/recovered_interactions_qc.csv; tables/residue_interaction_features_alpha1.csv; tables/residue_interaction_features_alpha2.csv; tables/residue_level_drug_comparison.csv; tables/aromatic_interaction_summary.csv",
    ]
    text = "\n".join(report) + "\n"
    (run_dir / "RESIDUE_HYPOTHESIS_REPORT.md").write_text(text)
    (run_dir / "FINAL_REPORT.md").write_text(text)
    qc_text = "\n".join(header + [
        "# QC_REPORT", "",
        f"Raw rows: {len(qc)}; recovered BZD-interface rows: {recovered_n}; unresolved after correction: {unresolved_n}.",
        f"Gamma2 C→E alignment optimal count: {gamma_meta['alignment_n_optimal']}; residue mapping rows: {len(corrected)}.",
        "Hydrophobic and pi_stack remain separate. Unresolved rows, if any, are retained with raw attributes.", "",
    ])
    (run_dir / "QC_REPORT.md").write_text(qc_text)
    (run_dir / "reports" / "RESIDUE_HYPOTHESIS_REPORT.md").write_text(text)


def run() -> dict[str, Any]:
    latest = run_manager.latest_run()
    if latest:
        latest_manifest = run_manager.RUNS / latest / "run_manifest.json"
        try:
            previous = json.loads(latest_manifest.read_text())
            if previous.get("qc_status") == "FAILED" and previous.get("source_run"):
                latest = previous["source_run"]
        except Exception:
            pass
    source_run = latest if latest and (run_manager.RUNS / latest / "raw" / "plip_interactions_4drug.csv").exists() else "legacy_unversioned_4drug_poc"
    source_dir = run_manager.RUNS / source_run if source_run != "legacy_unversioned_4drug_poc" else None
    source_files = [
        f"runs/{source_run}/raw/plip_interactions_4drug.csv" if source_dir else "results/tables/plip_interactions_4drug.csv",
        f"runs/{source_run}/raw/docking_results_4drug.csv" if source_dir else "results/tables/docking_results_4drug.csv",
        f"runs/{source_run}/raw/residue_mapping.csv" if source_dir else "outputs/residue_mapping.csv",
        "data/raw/structures/strict3/6HUP.pdb",
        "data/raw/structures/strict3/9CTJ.pdb",
    ]
    run_dir, manifest = run_manager.create_run(
        short_phase_name="residue_hypothesis",
        analysis_phase="BZD interface residue hypothesis validation",
        analysis_purpose="Recover chain-mapped BZD-interface interactions and compare residue-level drug profiles within each receptor.",
        parent_run=source_run,
        drugs=DRUGS,
        receptors=RECEPTORS,
        docking_parameters={"docking_executed": False, "source": source_run, "pose_counts": "copied existing docking results"},
        input_files=source_files,
        notes=[
            "Gamma2 chain C (6HUP) to chain E (9CTJ) mapping added by unique global sequence alignment.",
            "Hydrophobic and pi_stack are retained as separate interaction features.",
            "Alpha2 remains a provisional local construct; primary comparisons are within receptor.",
        ],
        source_run=source_run,
        source_files=source_files,
        code_paths=[Path(__file__), ROOT / "src" / "run_manager.py"],
    )
    raw_dir = run_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    if source_dir:
        for name in ("plip_interactions_4drug.csv", "docking_results_4drug.csv", "residue_mapping.csv"):
            run_manager.copy_into_run(source_dir / "raw" / name, raw_dir / name)
    else:
        for source in (ROOT / "results/tables/plip_interactions_4drug.csv", ROOT / "results/tables/docking_results_4drug.csv", ROOT / "outputs/residue_mapping.csv"):
            run_manager.copy_into_run(source, raw_dir / source.name)
    for source in (ROOT / "data/raw/structures/strict3/6HUP.pdb", ROOT / "data/raw/structures/strict3/9CTJ.pdb"):
        run_manager.copy_into_run(source, raw_dir / "structures" / "strict3" / source.name)
    manifest["input_files"] = sorted(set(manifest["input_files"] + [str(p.relative_to(run_dir)) for p in raw_dir.rglob("*") if p.is_file()]))
    run_manager.write_manifest(run_dir, manifest)

    corrected, mapping_index, gamma_meta = _corrected_mapping(run_dir, source_dir or ROOT)
    recovered, qc = _recover_interactions(run_dir, source_dir or ROOT, mapping_index)
    docking = pd.read_csv((source_dir / "raw/docking_results_4drug.csv") if source_dir else ROOT / "results/tables/docking_results_4drug.csv")
    features = {name: _feature_rows(recovered, docking, name) for name in ("alpha1", "alpha2")}
    features["alpha1"].to_csv(run_dir / "tables/residue_interaction_features_alpha1.csv", index=False)
    features["alpha2"].to_csv(run_dir / "tables/residue_interaction_features_alpha2.csv", index=False)
    comparison = _comparison(pd.concat([features["alpha1"], features["alpha2"]], ignore_index=True))
    comparison.to_csv(run_dir / "tables/residue_level_drug_comparison.csv", index=False)
    aromatic = pd.concat([features["alpha1"], features["alpha2"]], ignore_index=True)
    aromatic = aromatic[aromatic.residue_name.isin(AROMATIC)].copy()
    aromatic["site_label"] = aromatic.apply(lambda r: f"{r.interface_subunit}_TYR58" if r.common_position == "BZD_GAMMA2_058" and r.residue_name == "TYR" else (f"{r.interface_subunit}_PHE77" if r.common_position == "BZD_GAMMA2_077" and r.residue_name == "PHE" else f"{r.interface_subunit}_{r.residue_name}{str(r.common_position).split('_')[-1]}"), axis=1)
    aromatic["priority_site"] = aromatic.common_position.isin(["BZD_GAMMA2_058", "BZD_GAMMA2_077"])
    aromatic.to_csv(run_dir / "tables/aromatic_interaction_summary.csv", index=False)
    _report(run_dir, manifest, corrected, recovered, qc, features, comparison, aromatic, gamma_meta)
    finalized = run_manager.finalize_run(run_dir, manifest, qc_status="COMPLETED_WITH_QC")
    return {"run_id": finalized["run_id"], "run_dir": str(run_dir), "raw_rows": len(qc), "recovered_bzd_interface_rows": int((qc.recovery_status == "recovered_bzd_interface").sum()), "unresolved_after_correction": int((qc.recovery_status == "unresolved_after_correction").sum()), "alpha1_feature_rows": len(features["alpha1"]), "alpha2_feature_rows": len(features["alpha2"])}


def main() -> None:
    print(json.dumps(run(), ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
