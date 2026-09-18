"""Grouped and row-normalized views of the existing residue frequency matrix.

The source matrix is copied from an immutable parent run.  This phase only
changes ordering/visual scaling; it never rewrites the raw frequency values.
"""

from __future__ import annotations

import json
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
PARENT_RUN = "20260917_1755_grouped_heatmap"
DRUGS = ["diazepam", "alprazolam", "triazolam", "zolpidem"]
CATEGORY_ORDER = ["Aromatic", "Hydrophobic / aliphatic", "Hydrophilic / polar", "Charged", "Other"]
CATEGORY_RANK = {x: i for i, x in enumerate(CATEGORY_ORDER)}


def residue_class(residue: str) -> str:
    residue = str(residue).upper()
    if residue in {"PHE", "TYR", "TRP", "HIS"}:
        return "Aromatic"
    if residue in {"ALA", "VAL", "LEU", "ILE", "MET", "PRO"}:
        return "Hydrophobic / aliphatic"
    if residue in {"SER", "THR", "ASN", "GLN", "CYS"}:
        return "Hydrophilic / polar"
    if residue in {"LYS", "ARG", "ASP", "GLU"}:
        return "Charged"
    return "Other"


def _position_number(position: str) -> str:
    try:
        return str(int(str(position).split("_")[-1]))
    except Exception:
        return str(position)


def _load_matrix(source: Path) -> pd.DataFrame:
    # The grouped-heatmap parent keeps the unmodified source matrix under its
    # lineage directory.  Falling back to the parent bridge run supports
    # reproducibility if this module is run before a grouped view exists.
    path = source / "data" / "raw" / "lineage" / "residue_ligand_frequency_matrix.csv"
    if not path.exists():
        path = source / "results" / "tables" / "residue_ligand_frequency_matrix.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    d = pd.read_csv(path)
    required = {"receptor", "common_position", "residue_name", "interaction_type", *DRUGS}
    missing = required - set(d.columns)
    if missing:
        raise ValueError(f"frequency matrix missing columns: {sorted(missing)}")
    d["residue_class"] = d["residue_name"].map(residue_class)
    d["class_rank"] = d["residue_class"].map(CATEGORY_RANK)
    d["position_number"] = d["common_position"].map(_position_number)
    d["position_sort"] = pd.to_numeric(d["common_position"].str.extract(r"(\\d+)$")[0], errors="coerce")
    d["subunit_label"] = d.apply(lambda r: "gamma2" if str(r.common_position).startswith("BZD_GAMMA2_") else str(r.receptor), axis=1)
    d["residue_label"] = d.apply(lambda r: f"{r.subunit_label} {r.residue_name}{r.position_number}", axis=1)
    d["feature_label"] = d.apply(lambda r: f"{r.receptor} | {r.residue_label} | {r.residue_class} | {r.interaction_type}", axis=1)
    # Receptor is the primary block key. Within each block the requested
    # chemical-class order is used, then residue/common position and type.
    d["receptor_rank"] = d["receptor"].map({"alpha1": 0, "alpha2": 1}).fillna(9)
    d["interaction_sort"] = d["interaction_type"].astype(str)
    d = d.sort_values(["receptor_rank", "class_rank", "position_sort", "common_position", "interaction_sort"], kind="stable").reset_index(drop=True)
    d["row_order"] = np.arange(len(d))
    return d


def _normalized(d: pd.DataFrame) -> pd.DataFrame:
    out = d.copy()
    x = out[DRUGS].astype(float)
    out[DRUGS] = x.sub(x.mean(axis=1), axis=0)
    out["normalization"] = "row_mean_centered: raw frequency minus feature-row mean across four drugs"
    return out


def _plot(d: pd.DataFrame, out: Path, *, normalized: bool) -> None:
    plot_rows: list[pd.Series | None] = []
    labels: list[str] = []
    boundaries: list[tuple[float, str]] = []
    current_receptor = None
    current_class = None
    for _, row in d.iterrows():
        if current_receptor is not None and row.receptor != current_receptor:
            plot_rows.append(None); labels.append("")
            boundaries.append((len(plot_rows) - 0.5, f"{row.receptor}"))
            current_class = None
        if current_class is not None and row.residue_class != current_class:
            plot_rows.append(None); labels.append("")
        if current_receptor != row.receptor:
            current_receptor = row.receptor
            boundaries.append((len(plot_rows) - 0.5, current_receptor))
        current_class = row.residue_class
        plot_rows.append(row); labels.append(str(row.feature_label))
    # First pass above records receptor anchors; class spans are computed from
    # the displayed rows for labelled horizontal separators.
    arr = np.full((len(plot_rows), len(DRUGS)), np.nan, dtype=float)
    for i, row in enumerate(plot_rows):
        if row is not None:
            arr[i, :] = row[DRUGS].astype(float).to_numpy()
    fig_h = max(6.5, 0.34 * len(plot_rows) + 1.8)
    fig, ax = plt.subplots(figsize=(10.8, fig_h))
    cmap = plt.cm.RdBu_r.copy() if normalized else plt.cm.viridis.copy()
    cmap.set_bad("white")
    if normalized:
        vmax = float(np.nanmax(np.abs(arr))) if np.isfinite(arr).any() else 1.0
        vmax = max(vmax, 0.05)
        im = ax.imshow(np.ma.masked_invalid(arr), aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
        color_label = "row mean-centered frequency"
    else:
        im = ax.imshow(np.ma.masked_invalid(arr), aspect="auto", cmap=cmap, vmin=0, vmax=1)
        color_label = "raw interaction frequency"
    ax.set_xticks(range(len(DRUGS)), DRUGS, rotation=30, ha="right")
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    ax.tick_params(axis="y", length=0)
    # Separator rows are white; heavy lines make alpha1/alpha2 blocks explicit.
    for i, row in enumerate(plot_rows):
        if row is None:
            ax.axhline(i - 0.5, color="#444444", linewidth=0.8)
    ax.set_title("BZD-site residue × ligand interaction frequency — grouped by residue class" if not normalized else "BZD-site residue × ligand frequency — row mean-centered pattern")
    ax.set_ylabel("alpha1/alpha2 blocks; residue class → position → interaction type")
    fig.colorbar(im, ax=ax, label=color_label, fraction=0.03, pad=0.02)
    fig.subplots_adjust(left=0.46, right=0.92, bottom=0.09, top=0.96)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def _interpretation(d: pd.DataFrame) -> pd.DataFrame:
    x = d.copy()
    x["range"] = x[DRUGS].max(axis=1) - x[DRUGS].min(axis=1)
    x["high_drug"] = x[DRUGS].idxmax(axis=1)
    x["low_drug"] = x[DRUGS].idxmin(axis=1)
    x["pattern"] = x.apply(lambda r: ";".join(f"{drug}={float(r[drug]):.3f}" for drug in DRUGS), axis=1)
    groups: list[tuple[str, pd.DataFrame]] = []
    aromatic_pi = x[(x.residue_class == "Aromatic") & (x.interaction_type == "pi_stack")]
    if len(aromatic_pi) < 3:
        aromatic_pi = x[x.residue_class == "Aromatic"]
    groups.append(("aromatic_pi_stack_features", aromatic_pi))
    groups.append(("hydrophobic_interaction_features", x[x.interaction_type == "hydrophobic_interaction"]))
    polar_hbond = x[(x.residue_class == "Hydrophilic / polar") & (x.interaction_type == "hydrogen_bond")]
    if len(polar_hbond) < 3:
        polar_hbond = x[x.residue_class == "Hydrophilic / polar"]
    groups.append(("hydrophilic_or_hbond_features", polar_hbond))
    rows = []
    for group, g in groups:
        g = g.sort_values(["range", "receptor", "position_sort"], ascending=[False, True, True]).head(3)
        for _, r in g.iterrows():
            rows.append({
                "interpretation_group": group,
                "receptor": r.receptor,
                "common_position": r.common_position,
                "residue_name": r.residue_name,
                "residue_class": r.residue_class,
                "interaction_type": r.interaction_type,
                "frequency_range": float(r["range"]),
                "highest_frequency_drug": r.high_drug,
                "lowest_frequency_drug": r.low_drug,
                "drug_frequency_pattern": r.pattern,
                "interpretation_note": "descriptive within-row pattern; no total-interaction ranking or statistical inference",
            })
    return pd.DataFrame(rows)


def _report(run_dir: Path, manifest: dict[str, Any], d: pd.DataFrame, interpretation: pd.DataFrame) -> None:
    header = [
        f"Run ID: {manifest['run_id']}",
        f"Parent Run: {manifest.get('parent_run') or ''}",
        f"Analysis Phase: {manifest['analysis_phase']}",
        f"Started: {manifest['started_at']}",
        f"Completed: {run_manager.iso_now()}", "",
    ]
    report = header + [
        "# GROUPED_HEATMAP_REPORT", "",
        "The original residue_ligand_frequency_heatmap was not modified. This additional view uses the parent run's raw frequency matrix and changes only row ordering, labeling and the optional row mean-centered color scale.", "",
        "## Residue-class visualization", "",
        "Residues are classified for visualization only: Aromatic (PHE/TYR/TRP/HIS), Hydrophobic / aliphatic (ALA/VAL/LEU/ILE/MET/PRO), Hydrophilic / polar (SER/THR/ASN/GLN/CYS), Charged (LYS/ARG/ASP/GLU), and Other. Interaction type remains an independent feature dimension.", "",
        "Rows are grouped into separate alpha1 and alpha2 blocks. Within each block, residue class, common position, and interaction type determine order; interaction types for the same residue are adjacent. The primary heatmap retains raw 0–1 frequencies.", "",
        "## Row-normalized view", "",
        "The normalized heatmap uses row mean-centering across the four drugs: each raw frequency minus that row's four-drug mean. It is a visualization of relative pattern only; the raw matrix is preserved unchanged and should be used for numeric reporting.", "",
        "## Residue-specific pattern differences", "",
        interpretation.to_string(index=False) if not interpretation.empty else "No candidate rows available.", "",
        "These entries are selected by within-row frequency range. They are not ranked by total interaction count and do not imply statistical significance, subtype causality, or pharmacological potency.", "",
        "## Key-site visibility", "",
        "The row labels explicitly retain gamma2 TYR58, gamma2 PHE77, alpha-side HIS102, LYS156/157 positions where present, and SER205/206 positions where present. A row appears only when the corresponding feature was observed in the lossless parent matrix; no unobserved feature was added.", "",
        "## Files", "",
        "results/figures/residue_ligand_frequency_heatmap_grouped.png; results/figures/residue_ligand_frequency_heatmap_row_normalized.png; results/tables/residue_ligand_frequency_matrix_grouped.csv; results/tables/grouped_heatmap_top_differences.csv",
    ]
    text = "\n".join(report) + "\n"
    (run_dir / "reports" / "GROUPED_HEATMAP_REPORT.md").write_text(text)
    (run_dir / "FINAL_REPORT.md").write_text(text)
    (run_dir / "QC_REPORT.md").write_text("\n".join(header + ["# QC_REPORT", "", f"Rows preserved from parent matrix: {len(d)}.", "Raw frequencies were not modified. Row-normalized values are derived only for visualization.", ""]))


def run() -> dict[str, Any]:
    source = run_manager.RUNS / PARENT_RUN
    d = _load_matrix(source)
    run_dir, manifest = run_manager.create_run(
        short_phase_name="grouped_heatmap",
        analysis_phase="residue-class grouped fingerprint visualization",
        analysis_purpose="Improve interpretability of the existing residue-by-ligand frequency heatmap without modifying raw frequencies or the parent run.",
        parent_run=PARENT_RUN,
        drugs=DRUGS,
        receptors={"alpha1": "alpha1_beta3_gamma2", "alpha2": "alpha2_beta3_gamma2_local_9CTJ"},
        docking_parameters={"docking_executed": False, "source_run": PARENT_RUN, "row_normalization": "mean-centered across four drugs"},
        input_files=[f"runs/{PARENT_RUN}/data/raw/lineage/residue_ligand_frequency_matrix.csv"],
        notes=["Existing heatmap and parent run are immutable.", "Residue class is visualization-only and remains separate from interaction type.", "No unobserved candidate feature was added."],
        source_run=PARENT_RUN,
        source_files=["data/raw/lineage/residue_ligand_frequency_matrix.csv", "data/raw/lineage/candidate_residue_features.csv"],
        code_paths=[Path(__file__), ROOT / "src" / "run_manager.py"],
    )
    lineage = run_dir / "data" / "raw" / "lineage"
    matrix_source = source / "data" / "raw" / "lineage" / "residue_ligand_frequency_matrix.csv"
    if not matrix_source.exists():
        matrix_source = source / "results" / "tables" / "residue_ligand_frequency_matrix.csv"
    run_manager.copy_into_run(matrix_source, lineage / "residue_ligand_frequency_matrix.csv")
    feature_source = source / "data" / "raw" / "lineage" / "candidate_residue_features.csv"
    if not feature_source.exists():
        feature_source = source / "results" / "tables" / "candidate_residue_features.csv"
    if feature_source.exists():
        run_manager.copy_into_run(feature_source, lineage / "candidate_residue_features.csv")
    manifest["input_files"] = sorted(set(manifest["input_files"] + [str(p.relative_to(run_dir)) for p in lineage.iterdir()]))
    run_manager.write_manifest(run_dir, manifest)
    norm = _normalized(d)
    grouped = d[["receptor", "common_position", "residue_name", "residue_class", "class_rank", "position_number", "interaction_type", *DRUGS, *[f"{drug}_status" for drug in DRUGS], "feature_label", "row_order"]].copy()
    grouped.to_csv(run_dir / "results" / "tables" / "residue_ligand_frequency_matrix_grouped.csv", index=False)
    # Save the derived row-centered values in a separate table; raw grouped CSV remains untouched.
    norm.to_csv(run_dir / "results" / "tables" / "residue_ligand_frequency_matrix_row_normalized.csv", index=False)
    interpretation = _interpretation(d)
    interpretation.to_csv(run_dir / "results" / "tables" / "grouped_heatmap_top_differences.csv", index=False)
    _plot(d, run_dir / "results" / "figures" / "residue_ligand_frequency_heatmap_grouped.png", normalized=False)
    _plot(norm, run_dir / "results" / "figures" / "residue_ligand_frequency_heatmap_row_normalized.png", normalized=True)
    _report(run_dir, manifest, d, interpretation)
    finalized = run_manager.finalize_run(run_dir, manifest, qc_status="COMPLETED_WITH_QC")
    return {"run_id": finalized["run_id"], "run_dir": str(run_dir), "rows": len(d), "interpretation_rows": len(interpretation)}


def main() -> None:
    print(json.dumps(run(), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
