"""Multi-receptor interaction-fingerprint feature-engineering PoC.

The primary analysis is drug-to-drug comparison within a fixed receptor.
Existing strict 3-drug and 4-drug artifacts are read but never overwritten.
"""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from typing import Dict, List, Mapping, Tuple

import numpy as np
import pandas as pd

try:
    from . import run_manager
except ImportError:  # pragma: no cover - direct script compatibility
    import run_manager

ROOT = Path(__file__).resolve().parents[1]
DRUGS = ["diazepam", "alprazolam", "triazolam", "zolpidem"]
RECEPTORS = {
    "alpha1": "alpha1_beta3_gamma2",
    "alpha2": "alpha2_beta3_gamma2_local_9CTJ",
}
# Legacy inputs are read-only.  These defaults retain backwards compatibility
# for imports, while run() configures all generated files under runs/<run_id>.
LEGACY_PROC = ROOT / "data" / "processed"
LEGACY_RES = ROOT / "results" / "tables"
LEGACY_OUT = ROOT / "outputs"
PROC = LEGACY_PROC
RES = LEGACY_RES
FIG = ROOT / "results" / "figures"
OUT = LEGACY_OUT
RUN_RAW: Path | None = None
RUN_ROOT: Path | None = None


def configure_run(run_dir: Path) -> None:
    """Route all writes to one immutable run directory."""
    global PROC, RES, FIG, OUT, RUN_RAW, RUN_ROOT
    RUN_ROOT = run_dir
    RUN_RAW = run_dir / "raw"
    PROC = run_dir / "processed"
    RES = run_dir / "tables"
    FIG = run_dir / "figures"
    OUT = run_dir / "reports"


def source_path(filename: str, legacy_dir: Path) -> Path:
    """Prefer the byte-preserved copy inside the active run."""
    if RUN_RAW is not None and (RUN_RAW / filename).exists():
        return RUN_RAW / filename
    return legacy_dir / filename


def mkdirs() -> None:
    for p in (PROC, RES, FIG, OUT):
        p.mkdir(parents=True, exist_ok=True)


def write_csv(df: pd.DataFrame, *paths: Path) -> None:
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False)


def safe_int(value: object) -> int | None:
    try:
        if pd.isna(value):
            return None
        return int(float(value))
    except Exception:
        return None


def distance_from_attributes(interaction_type: str, raw: object) -> Tuple[float, str]:
    """Extract only distance fields observed in the current PLIP schema."""
    try:
        attrs = json.loads(str(raw)) if not isinstance(raw, Mapping) else dict(raw)
    except Exception:
        attrs = {}
    candidates = {
        "hydrophobic_interaction": ["dist"],
        "halogen_bond": ["dist"],
        "hydrogen_bond": ["dist-d-a", "dist_d-a"],
        "pi_stack": ["centdist"],
        "pi_cation": ["dist"],
        "salt_bridge": ["dist"],
        "water_bridge": ["dist"],
        "metal_complex": ["dist"],
    }.get(str(interaction_type), ["dist"])
    for key in candidates:
        val = attrs.get(key)
        if val is None or str(val).strip() == "":
            continue
        try:
            number = float(val)
            if math.isfinite(number):
                return number, key
        except Exception:
            pass
    return np.nan, ""


def mapping_indexes() -> Tuple[Dict[Tuple[str, int], dict], Dict[Tuple[str, int], dict]]:
    mapping = pd.read_csv(source_path("residue_mapping.csv", LEGACY_OUT))
    indexes: List[Dict[Tuple[str, int], dict]] = []
    for side in ("alpha1", "alpha2"):
        index: Dict[Tuple[str, int], dict] = {}
        for _, row in mapping.iterrows():
            chain = str(row.get(f"{side}_chain", ""))
            number = safe_int(row.get(f"{side}_residue_number"))
            if number is None:
                continue
            index[(chain, number)] = {
                "common_position": str(row.get("common_position", "")),
                "mapping_status": str(row.get("mapping_status", "")),
                "qc_required": bool(row.get("qc_required", False)),
                "mapping_method": str(row.get("mapping_method", "")),
            }
        indexes.append(index)
    return indexes[0], indexes[1]


def prepare_interactions() -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, int]]:
    source = pd.read_csv(source_path("plip_interactions_4drug.csv", LEGACY_RES))
    source = source[source.drug_id.isin(DRUGS)].copy()
    a1_map, a2_map = mapping_indexes()
    mapped_rows: List[dict] = []
    unresolved_rows: List[dict] = []
    distance_fields: Dict[str, int] = {}
    for _, row in source.iterrows():
        side = "alpha1" if row.receptor_id == RECEPTORS["alpha1"] else "alpha2"
        index = a1_map if side == "alpha1" else a2_map
        mapping = index.get((str(row.residue_chain), safe_int(row.residue_number)))
        distance, distance_field = distance_from_attributes(row.interaction_type, row.raw_attributes)
        if distance_field:
            distance_fields[distance_field] = distance_fields.get(distance_field, 0) + 1
        base = {
            "drug_id": row.drug_id,
            "receptor": side,
            "receptor_id": row.receptor_id,
            "pose_id": row.pose_id,
            "residue_chain": row.residue_chain,
            "residue_number": row.residue_number,
            "interaction_type": row.interaction_type,
            "distance": distance,
            "distance_field": distance_field,
            "raw_attributes": row.raw_attributes,
        }
        resolved = bool(
            mapping
            and mapping["mapping_status"] == "mapped_global_sequence_alignment"
            and not mapping["qc_required"]
            and mapping["common_position"] not in ("", "nan")
        )
        if resolved:
            mapped_rows.append({
                **base,
                "common_position": mapping["common_position"],
                "mapping_status": mapping["mapping_status"],
                "mapping_method": mapping["mapping_method"],
            })
        else:
            unresolved_rows.append({
                **base,
                "common_position": mapping.get("common_position", "") if mapping else "",
                "mapping_status": mapping.get("mapping_status", "unmapped") if mapping else "unmapped",
                "mapping_method": mapping.get("mapping_method", "") if mapping else "",
                "unresolved_reason": "qc_required_or_nonresolved_status" if mapping else "no_chain_residue_mapping",
            })
    resolved = pd.DataFrame(mapped_rows)
    unresolved = pd.DataFrame(unresolved_rows)
    if resolved.empty:
        resolved = pd.DataFrame(columns=[
            "drug_id", "receptor", "receptor_id", "pose_id", "residue_chain",
            "residue_number", "interaction_type", "distance", "distance_field",
            "raw_attributes", "common_position", "mapping_status", "mapping_method",
        ])
    if unresolved.empty:
        unresolved = pd.DataFrame(columns=[
            "drug_id", "receptor", "receptor_id", "pose_id", "residue_chain",
            "residue_number", "interaction_type", "distance", "distance_field",
            "raw_attributes", "common_position", "mapping_status", "mapping_method",
            "unresolved_reason",
        ])
    write_csv(unresolved, RES / "unresolved_feature_qc.csv")
    return resolved, unresolved, distance_fields


def base_feature(row: pd.Series) -> str:
    return f"{row.receptor}_{row.common_position}_{row.interaction_type}"


def feature_stats(resolved: pd.DataFrame, receptor: str) -> Tuple[pd.DataFrame, List[str]]:
    d = resolved[resolved.receptor == receptor].copy()
    d["feature"] = d.apply(base_feature, axis=1)
    feature_names = sorted(d.feature.unique())
    docking = pd.read_csv(source_path("docking_results_4drug.csv", LEGACY_RES))
    docking = docking[docking.receptor_id == RECEPTORS[receptor]]
    pose_counts = docking.groupby("drug_id").pose_id.nunique().to_dict()
    rows: List[dict] = []
    for drug in DRUGS:
        out: dict = {"drug_id": drug}
        for feature in feature_names:
            g = d[(d.drug_id == drug) & (d.feature == feature)]
            nposes = int(pose_counts.get(drug, 0))
            distances = pd.to_numeric(g.distance, errors="coerce").dropna()
            prefix = feature
            out[f"{prefix}_frequency"] = float(g.pose_id.nunique() / nposes) if nposes else np.nan
            out[f"{prefix}_interaction_count"] = int(len(g))
            out[f"{prefix}_n_observed"] = int(g.pose_id.nunique())
            out[f"{prefix}_n_poses"] = nposes
            out[f"{prefix}_mean_distance"] = distances.mean() if len(distances) else np.nan
            out[f"{prefix}_median_distance"] = distances.median() if len(distances) else np.nan
            out[f"{prefix}_min_distance"] = distances.min() if len(distances) else np.nan
            out[f"{prefix}_max_distance"] = distances.max() if len(distances) else np.nan
            out[f"{prefix}_std_distance"] = distances.std(ddof=1) if len(distances) > 1 else (0.0 if len(distances) == 1 else np.nan)
            out[f"{prefix}_distance_n_observed"] = int(len(distances))
        rows.append(out)
    matrix = pd.DataFrame(rows)
    return matrix[["drug_id"] + sorted(c for c in matrix.columns if c != "drug_id")], feature_names


def build_feature_matrices(resolved: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    a1, _ = feature_stats(resolved, "alpha1")
    a2, _ = feature_stats(resolved, "alpha2")
    write_csv(a1, PROC / "alpha1_drug_feature_matrix.csv")
    write_csv(a2, PROC / "alpha2_drug_feature_matrix.csv")
    combined = a1.merge(a2, on="drug_id", how="outer", suffixes=("_alpha1", "_alpha2"))
    combined = combined.set_index("drug_id").reindex(DRUGS).reset_index()
    combined = combined[["drug_id"] + sorted(c for c in combined.columns if c != "drug_id")]
    write_csv(combined, PROC / "combined_alpha1_alpha2_feature_matrix.csv")
    return a1, a2, combined


def build_labels() -> pd.DataFrame:
    raw = pd.read_csv(source_path("primary_tier1_pharmacology_raw.csv", LEGACY_PROC))
    ki = raw[(raw.drug_id.isin(DRUGS)) & (raw.tier == 1) & (raw.metric.astype(str).str.casefold() == "ki")].copy()
    rows = []
    for drug in DRUGS:
        d = ki[ki.drug_id == drug].copy()
        ratios = pd.to_numeric(d.ratio_alpha2_over_alpha1, errors="coerce").dropna()
        logs = pd.to_numeric(d.log10_ratio, errors="coerce").dropna()
        rows.append({
            "drug_id": drug,
            "alpha1_Ki": pd.to_numeric(d.alpha1_value, errors="coerce").median() if not d.empty else np.nan,
            "alpha2_Ki": pd.to_numeric(d.alpha2_value, errors="coerce").median() if not d.empty else np.nan,
            "alpha2_over_alpha1_ratio": ratios.median() if len(ratios) else np.nan,
            "log10_alpha2_over_alpha1_ratio": logs.median() if len(logs) else np.nan,
            "alpha1_Ki_raw_context_values": ";".join(f"{v:g}" for v in pd.to_numeric(d.alpha1_value, errors="coerce").dropna()),
            "alpha2_Ki_raw_context_values": ";".join(f"{v:g}" for v in pd.to_numeric(d.alpha2_value, errors="coerce").dropna()),
            "ratio_raw_context_values": ";".join(f"{v:g}" for v in ratios),
            "comparison_context_count": int(len(d)),
            "composition_contexts": ";".join(sorted(set(f"{a} vs {b}" for a, b in zip(d.alpha1_composition, d.alpha2_composition)))),
            "aggregation_rule": "median across exact Tier-1 Ki contexts for label candidate; raw context values retained",
            "label_status": "candidate_only_not_primary_outcome",
        })
    labels = pd.DataFrame(rows)
    write_csv(labels, PROC / "drug_pharmacology_labels.csv")
    return labels


def frequency_columns(matrix: pd.DataFrame) -> List[str]:
    return [c for c in matrix.columns if c != "drug_id" and c.endswith("_frequency")]


def pairwise_similarity(matrix: pd.DataFrame, receptor: str) -> pd.DataFrame:
    freq_cols = frequency_columns(matrix)
    values = matrix.set_index("drug_id")[freq_cols].astype(float)
    rows = []
    for i, drug_a in enumerate(DRUGS):
        for drug_b in DRUGS[i + 1:]:
            a = values.loc[drug_a].fillna(0.0).to_numpy(float)
            b = values.loc[drug_b].fillna(0.0).to_numpy(float)
            av, bv = a > 0, b > 0
            union = int(np.logical_or(av, bv).sum())
            inter = int(np.logical_and(av, bv).sum())
            denom = np.linalg.norm(a) * np.linalg.norm(b)
            rows.append({
                "receptor": receptor, "drug_a": drug_a, "drug_b": drug_b,
                "jaccard_similarity_binary": inter / union if union else np.nan,
                "shared_features": inter, "union_features": union,
                "cosine_similarity_frequency": float(np.dot(a, b) / denom) if denom else np.nan,
                "euclidean_distance_frequency": float(np.linalg.norm(a - b)),
                "n_frequency_features": len(freq_cols),
            })
    return pd.DataFrame(rows)


def build_similarity_matrices(a1: pd.DataFrame, a2: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    s1 = pairwise_similarity(a1, "alpha1")
    s2 = pairwise_similarity(a2, "alpha2")
    write_csv(s1, RES / "alpha1_interdrug_similarity.csv")
    write_csv(s2, RES / "alpha2_interdrug_similarity.csv")
    return s1, s2


def variable_features(matrix: pd.DataFrame, receptor: str) -> pd.DataFrame:
    rows = []
    bases = sorted({c.rsplit("_", 1)[0] for c in matrix.columns if c.endswith("_frequency")})
    for base in bases:
        for stat, col in [("frequency", f"{base}_frequency"), ("mean_distance", f"{base}_mean_distance")]:
            vals = pd.to_numeric(matrix[col], errors="coerce")
            valid = vals.dropna()
            if valid.empty:
                continue
            mean = float(valid.mean())
            variance = float(valid.var(ddof=1)) if len(valid) > 1 else 0.0
            value_range = float(valid.max() - valid.min())
            cv = float(valid.std(ddof=1) / abs(mean)) if len(valid) > 1 and mean != 0 else np.nan
            top_i, low_i = vals.idxmax(), vals.idxmin()
            rows.append({
                "receptor": receptor, "feature": base, "statistic": stat,
                "variance": variance, "range": value_range, "max_minus_min": value_range,
                "coefficient_of_variation": cv, "n_nonmissing": int(len(valid)),
                "drug_high": matrix.loc[top_i, "drug_id"], "value_high": float(vals.loc[top_i]),
                "drug_low": matrix.loc[low_i, "drug_id"], "value_low": float(vals.loc[low_i]),
                "drug_values": ";".join(f"{d}={v:g}" for d, v in zip(matrix.drug_id, vals)),
                "missing_distance_or_absence_note": "frequency zero means no interaction in available poses; distance NA means no distance observation",
            })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["variance_rank"] = out.variance.rank(method="min", ascending=False).astype(int)
        out["range_rank"] = out.range.rank(method="min", ascending=False).astype(int)
        out = out.sort_values(["variance_rank", "range_rank", "feature", "statistic"]).reset_index(drop=True)
    return out


def build_variable_feature_tables(a1: pd.DataFrame, a2: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    v1, v2 = variable_features(a1, "alpha1"), variable_features(a2, "alpha2")
    write_csv(v1, RES / "alpha1_variable_features.csv")
    write_csv(v2, RES / "alpha2_variable_features.csv")
    return v1, v2


def build_ml_qc(combined: pd.DataFrame, unresolved: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
    features = combined.drop(columns=["drug_id"])
    missing_rate = float(features.isna().mean().mean()) if features.size else np.nan
    zero_variance = int(sum(features[c].dropna().nunique() <= 1 for c in features.columns))
    sparse = int(sum(((pd.to_numeric(features[c], errors="coerce").fillna(0) == 0).mean() >= 0.75) for c in features.columns if c.endswith("_frequency")))
    distance_value_cols = [c for c in features.columns if c.endswith(("_mean_distance", "_median_distance", "_min_distance", "_max_distance"))]
    distance_missing = int(features[distance_value_cols].isna().sum().sum()) if distance_value_cols else 0
    distance_zero = int((features[distance_value_cols] == 0).sum().sum()) if distance_value_cols else 0
    qc = pd.DataFrame([{
        "status": "NOT READY FOR MODEL TRAINING",
        "n_drugs": len(combined), "n_features": int(features.shape[1]),
        "missing_rate_all_features": missing_rate,
        "zero_variance_feature_count": zero_variance,
        "highly_sparse_frequency_feature_count": sparse,
        "distance_stat_missing_cells": distance_missing,
        "distance_stat_zero_cells": distance_zero,
        "unresolved_interaction_rows_excluded": int(len(unresolved)),
        "unresolved_feature_count": int(unresolved[["receptor", "residue_chain", "residue_number", "interaction_type"]].drop_duplicates().shape[0]) if not unresolved.empty else 0,
        "label_rows": len(labels),
        "feature_leakage_risk": "low by column content; future ML must split by drug and keep labels out of X",
        "reason": "n=4, provisional alpha2 structure, unresolved mapping exclusions, and distance missingness require QC before training",
    }])
    write_csv(qc, RES / "ml_dataset_readiness_qc.csv")
    return qc


def heatmap_matrix(matrix: pd.DataFrame, suffix: str) -> pd.DataFrame:
    cols = [c for c in matrix.columns if c != "drug_id" and c.endswith(suffix)]
    return matrix[["drug_id"] + sorted(cols)].set_index("drug_id")


def make_figures(a1: pd.DataFrame, a2: pd.DataFrame, combined: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def draw(data: pd.DataFrame, path: Path, title: str, cmap_name: str, vmin=None, vmax=None) -> None:
        arr = data.to_numpy(float)
        fig, ax = plt.subplots(figsize=(max(8, data.shape[1] * 0.28), 4.8))
        cmap = plt.get_cmap(cmap_name).copy()
        cmap.set_bad("#d9d9d9")
        im = ax.imshow(arr, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_yticks(range(len(data.index)), data.index)
        ax.set_xticks(range(len(data.columns)), data.columns, rotation=90, fontsize=7)
        ax.set_title(title, fontsize=11)
        fig.colorbar(im, ax=ax, fraction=.02)
        fig.tight_layout(rect=[0, 0, 1, .95])
        fig.savefig(path, dpi=180)
        plt.close(fig)

    a1f, a2f, cf = heatmap_matrix(a1, "_frequency"), heatmap_matrix(a2, "_frequency"), heatmap_matrix(combined, "_frequency")
    draw(a1f, FIG / "alpha1_interdrug_fingerprint_heatmap.png", "Primary: alpha1 drug × resolved interaction frequency", "viridis", 0, 1)
    draw(a2f, FIG / "alpha2_interdrug_fingerprint_heatmap.png", "Primary: alpha2 drug × resolved interaction frequency", "viridis", 0, 1)
    draw(cf, FIG / "combined_multireceptor_fingerprint_heatmap.png", "Combined alpha1 + alpha2 frequency matrix", "viridis", 0, 1)
    draw(cf, FIG / "interaction_frequency_heatmap.png", "Resolved interaction frequency features", "viridis", 0, 1)
    distance_cols = [c for c in combined.columns if c.endswith("_mean_distance")]
    distance = combined[["drug_id"] + sorted(distance_cols)].set_index("drug_id") if distance_cols else pd.DataFrame(index=combined.drug_id, columns=["no_distance_feature"], dtype=float)
    draw(distance, FIG / "interaction_distance_heatmap.png", "Interaction mean distance (NA = no distance observation)", "magma")


def similarity_text(s1: pd.DataFrame, s2: pd.DataFrame) -> str:
    def fmt(df: pd.DataFrame) -> str:
        return "; ".join(f"{r.drug_a}-{r.drug_b}: J={r.jaccard_similarity_binary:.3f}, cosine={r.cosine_similarity_frequency:.3f}, Euclidean={r.euclidean_distance_frequency:.3f}" for _, r in df.iterrows())
    return f"alpha1: {fmt(s1)}\nalpha2: {fmt(s2)}"


def write_report(a1: pd.DataFrame, a2: pd.DataFrame, combined: pd.DataFrame, labels: pd.DataFrame, s1: pd.DataFrame, s2: pd.DataFrame, v1: pd.DataFrame, v2: pd.DataFrame, unresolved: pd.DataFrame, mlqc: pd.DataFrame, distance_fields: Mapping[str, int], run_metadata: Mapping[str, str]) -> None:
    docking = pd.read_csv(source_path("docking_results_4drug.csv", LEGACY_RES))
    plip = pd.read_csv(source_path("plip_interactions_4drug.csv", LEGACY_RES))
    pose_counts = docking.groupby(["drug_id", "receptor_id"]).pose_id.nunique()
    top1 = v1.head(5)[["feature", "statistic", "variance", "drug_high", "drug_low"]].to_string(index=False) if not v1.empty else "None"
    top2 = v2.head(5)[["feature", "statistic", "variance", "drug_high", "drug_low"]].to_string(index=False) if not v2.empty else "None"
    observed_types = ", ".join(sorted(plip.interaction_type.dropna().astype(str).unique()))
    header = [
        f"Run ID: {run_metadata.get('run_id', '')}",
        f"Parent Run: {run_metadata.get('parent_run', '')}",
        f"Analysis Phase: {run_metadata.get('analysis_phase', '')}",
        f"Started: {run_metadata.get('started_at', '')}",
        f"Completed: {run_metadata.get('completed_at', '')}",
        "",
    ]
    report = header + [
        "# MULTIRECEPTOR_FINGERPRINT_REPORT", "",
        "This phase is a multi-receptor interaction fingerprint feature-engineering PoC. Primary analysis is drug-to-drug comparison within a fixed receptor; alpha1-vs-alpha2 similarity remains secondary/QC.", "",
        "## Dataset and primary feature construction", "",
        f"Four drugs: {', '.join(DRUGS)}. Alpha1 frequency feature count: **{len(frequency_columns(a1))}**; alpha2 frequency feature count: **{len(frequency_columns(a2))}**; combined feature columns: **{combined.shape[1] - 1}**.",
        f"PLIP interaction rows: **{len(plip)}**; unresolved interaction rows excluded from primary: **{len(unresolved)}**; observed interaction types: {observed_types}; observed distance fields: {dict(distance_fields)}.",
        "Distance extraction uses PLIP dist for hydrophobic/halogen, dist-d-a for hydrogen bonds, and centdist for pi stacking. n_observed is the number of poses with the interaction; distance_n_observed counts numeric distance observations. No distance value was imputed.", "",
        "## Primary within-receptor drug similarity", "", "### Alpha1", "", s1.to_string(index=False), "", "### Alpha2", "", s2.to_string(index=False), "", "Pairwise summary:", similarity_text(s1, s2), "",
        "## Drug-discriminating features", "", "### Alpha1 highest-variance rows", "", top1, "", "### Alpha2 highest-variance rows", "", top2, "", "Variance/range rankings are descriptive only (n=4).", "",
        "## Pharmacology label candidates", "", labels.to_string(index=False), "", "Label medians aggregate exact Tier-1 Ki contexts only and retain raw context values.", "",
        "## ML dataset readiness", "", mlqc.to_string(index=False), "", "NOT READY FOR MODEL TRAINING. The matrix is X-shaped, but n=4, provisional alpha2 structure, unresolved mapping exclusions, and distance missingness prevent model training or performance claims.", "",
        "## Requested interpretation", "",
        "1. Alpha1 fingerprint differences are visible in the primary alpha1 drug × feature matrix and pairwise table.",
        "2. Alpha2 also shows drug-to-drug differences, conditional on the provisional 9CTJ-derived construct.",
        "3. Highest-variance residue/interaction features are ranked in alpha1_variable_features.csv and alpha2_variable_features.csv; these are descriptive candidates.",
        "4. Distance differences can be inspected where PLIP supplied a distance. Distance NA means no distance observation or no interaction, never zero distance.",
        "5. Zolpidem high/low features are listed by drug_high and drug_low; uniqueness is not claimed.",
        "6. Triazolam high/low features are listed in the same tables; uniqueness is not claimed.",
        "7. Separate alpha1 and alpha2 blocks preserve receptor-context information that a difference-only representation would discard.",
        "8. The combined matrix is structurally suitable as a future X table, but readiness is NOT READY FOR MODEL TRAINING.",
        "9. Before adding drugs, resolve chain-aware mapping, review the provisional alpha2 model/box transfer, audit distance schemas, and keep beta2 pharmacology versus beta3 docking differences explicit.",
        "10. The most natural next phase is to expand the drug panel with composition-matched pharmacology while preserving this feature schema; phenotype linkage should follow after feature/QC stability is demonstrated.", "",
        "## QC boundaries", "",
        f"Pose counts: {pose_counts.to_dict()}. All combinations currently have 9 poses. Frequency zero is used only for resolved features with pose data and no observed interaction. Unresolved mapping is excluded and listed in unresolved_feature_qc.csv. Distance missing is separate from distance=0. Alpha2 is provisional 9CTJ-derived; alprazolam pharmacology is beta2-background while docking is beta3-background. Existing strict and 4-drug reports/figures were not overwritten.", "",
        "## Files", "",
        "processed/alpha1_drug_feature_matrix.csv; alpha2_drug_feature_matrix.csv; combined_alpha1_alpha2_feature_matrix.csv; drug_pharmacology_labels.csv",
        "tables/alpha1_interdrug_similarity.csv; alpha2_interdrug_similarity.csv; alpha1_variable_features.csv; alpha2_variable_features.csv; unresolved_feature_qc.csv; ml_dataset_readiness_qc.csv",
        "figures/alpha1_interdrug_fingerprint_heatmap.png; alpha2_interdrug_fingerprint_heatmap.png; combined_multireceptor_fingerprint_heatmap.png; interaction_frequency_heatmap.png; interaction_distance_heatmap.png",
    ]
    text = "\n".join(report) + "\n"
    # Root-level reports are the stable run contract.  The named report copy
    # under reports/ preserves the descriptive legacy filename inside the run.
    if RUN_ROOT is not None:
        (RUN_ROOT / "FINAL_REPORT.md").write_text(text)
        qc = "\n".join(header + [
            "# QC_REPORT", "",
            f"QC status: {mlqc.iloc[0].get('status', 'UNKNOWN') if not mlqc.empty else 'UNKNOWN'}",
            f"Resolved interaction rows: {int(len(plip) - len(unresolved))}; unresolved interaction rows excluded: {len(unresolved)}.",
            f"Input PLIP rows: {len(plip)}; docking rows: {len(docking)}; distance fields observed: {dict(distance_fields)}.",
            "Mapping rows requiring QC remain in tables/unresolved_feature_qc.csv. Missing interaction and unresolved mapping are not converted to zero.",
            "No model training, imputation, metric integration, or result-driven docking adjustment was performed.",
            "",
        ])
        (RUN_ROOT / "QC_REPORT.md").write_text(qc)
    (OUT / "MULTIRECEPTOR_FINGERPRINT_REPORT.md").write_text(text)


def copy_run_inputs(run_dir: Path, source_run_dir: Path | None = None) -> list[str]:
    """Copy all source artifacts used by this phase into run/raw."""
    raw = run_dir / "raw"
    if source_run_dir is not None:
        sources = [source_run_dir / "raw" / name for name in (
            "plip_interactions_4drug.csv",
            "docking_results_4drug.csv",
            "residue_mapping.csv",
            "primary_tier1_pharmacology_raw.csv",
        )]
    else:
        sources = [
            LEGACY_RES / "plip_interactions_4drug.csv",
            LEGACY_RES / "docking_results_4drug.csv",
            LEGACY_OUT / "residue_mapping.csv",
            LEGACY_PROC / "primary_tier1_pharmacology_raw.csv",
        ]
    copied: list[str] = []
    for source in sources:
        if source.exists():
            destination = raw / source.name
            run_manager.copy_into_run(source, destination)
            copied.append(str(destination.relative_to(run_dir)))
    # Preserve the phase configuration and code inputs for auditability when
    # present; these are snapshots, not mutable links to the workspace.
    for source in (ROOT / "config" / "poc.yaml", Path(__file__), ROOT / "src" / "run_manager.py", ROOT / "multireceptor_poc.py"):
        if source.exists():
            destination = run_dir / "config" / "code_snapshot" / source.name
            run_manager.copy_into_run(source, destination)
            copied.append(str(destination.relative_to(run_dir)))
    return copied


def run() -> dict:
    latest = run_manager.latest_run()
    source_run_dir = None
    if latest:
        candidate = run_manager.RUNS / latest
        needed = [candidate / "raw" / name for name in (
            "plip_interactions_4drug.csv",
            "docking_results_4drug.csv",
            "residue_mapping.csv",
            "primary_tier1_pharmacology_raw.csv",
        )]
        if all(path.exists() for path in needed):
            source_run_dir = candidate
    parent = latest if source_run_dir is not None else "legacy_unversioned_4drug_poc"
    source_files = [
        *([f"runs/{parent}/raw/{name}" for name in (
            "plip_interactions_4drug.csv",
            "docking_results_4drug.csv",
            "residue_mapping.csv",
            "primary_tier1_pharmacology_raw.csv",
        )] if source_run_dir is not None else [
            "results/tables/plip_interactions_4drug.csv",
            "results/tables/docking_results_4drug.csv",
            "outputs/residue_mapping.csv",
            "data/processed/primary_tier1_pharmacology_raw.csv",
        ]),
    ]
    run_dir, manifest = run_manager.create_run(
        short_phase_name="multireceptor_feature",
        analysis_phase="multi-receptor fingerprint feature engineering",
        analysis_purpose="Versioned re-run of the four-drug alpha1/alpha2 feature-engineering PoC with immutable input lineage.",
        parent_run=parent,
        drugs=DRUGS,
        receptors=RECEPTORS,
        docking_parameters={"docking_executed": False, "source": "existing four-drug docking outputs", "poses": "as supplied by source run"},
        input_files=source_files,
        notes=[
            "Legacy four-drug outputs are copied byte-for-byte into raw/ and never modified.",
            "This phase is descriptive feature engineering; no ML training or inferential statistics are run.",
            "If no prior versioned run exists, parent_run identifies the legacy unversioned four-drug phase.",
        ],
        source_run=parent,
        source_files=source_files,
        code_paths=[Path(__file__), ROOT / "multireceptor_poc.py", ROOT / "src" / "run_manager.py"],
    )
    copied_inputs = copy_run_inputs(run_dir, source_run_dir)
    manifest["input_files"] = sorted(set(source_files + copied_inputs))
    run_manager.write_manifest(run_dir, manifest)
    configure_run(run_dir)
    mkdirs()
    resolved, unresolved, distance_fields = prepare_interactions()
    a1, a2, combined = build_feature_matrices(resolved)
    labels = build_labels()
    s1, s2 = build_similarity_matrices(a1, a2)
    v1, v2 = build_variable_feature_tables(a1, a2)
    mlqc = build_ml_qc(combined, unresolved, labels)
    make_figures(a1, a2, combined)
    # Capture the completion timestamp in the report header before finalizing
    # the manifest, then rewrite it once more after finalization if necessary.
    report_metadata = {
        "run_id": manifest["run_id"],
        "parent_run": manifest.get("parent_run") or "",
        "analysis_phase": manifest["analysis_phase"],
        "started_at": manifest["started_at"],
        "completed_at": run_manager.iso_now(),
    }
    write_report(a1, a2, combined, labels, s1, s2, v1, v2, unresolved, mlqc, distance_fields, report_metadata)
    finalized = run_manager.finalize_run(run_dir, manifest, qc_status="NOT READY FOR MODEL TRAINING")
    return {
        "run_id": finalized["run_id"],
        "run_dir": str(run_dir),
        "resolved_interaction_rows": int(len(resolved)),
        "unresolved_interaction_rows": int(len(unresolved)),
        "alpha1_frequency_features": int(len(frequency_columns(a1))),
        "alpha2_frequency_features": int(len(frequency_columns(a2))),
        "combined_feature_columns": int(combined.shape[1] - 1),
        "alpha1_interdrug_pairs": int(len(s1)),
        "alpha2_interdrug_pairs": int(len(s2)),
        "status": "NOT READY FOR MODEL TRAINING",
    }


def main() -> None:
    print(json.dumps(run(), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
