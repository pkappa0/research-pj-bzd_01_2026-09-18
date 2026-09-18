"""Residue-class-first heatmaps for the existing raw frequency matrix.

This is a visualization-only phase.  The input frequencies are copied from the
parent run and remain unchanged.  Rows are ordered by residue class, then
residue/site, receptor, and interaction type so alpha1/alpha2 contacts for the
same site are directly adjacent.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from . import run_manager
except ImportError:  # pragma: no cover
    import run_manager


ROOT = Path(__file__).resolve().parents[1]
PARENT_RUN = "20260917_1925_residue_class_heatmap"
DRUGS = ["diazepam", "alprazolam", "triazolam", "zolpidem"]
CATEGORIES = ["Aromatic", "Hydrophobic / aliphatic", "Hydrophilic / polar", "Charged", "Other"]
CATEGORY_RANK = {c: i for i, c in enumerate(CATEGORIES)}
INTERACTION_SHORT = {
    "hydrophobic_interaction": "hydrophobic",
    "pi_stack": "pi-stack",
    "hydrogen_bond": "H-bond",
    "halogen_bond": "halogen",
    "salt_bridge": "salt-bridge",
    "water_bridge": "water-bridge",
}


def _class(residue: str) -> str:
    r = str(residue).upper()
    if r in {"PHE", "TYR", "TRP", "HIS"}:
        return "Aromatic"
    if r in {"ALA", "VAL", "LEU", "ILE", "MET", "PRO"}:
        return "Hydrophobic / aliphatic"
    if r in {"SER", "THR", "ASN", "GLN", "CYS"}:
        return "Hydrophilic / polar"
    if r in {"LYS", "ARG", "ASP", "GLU"}:
        return "Charged"
    return "Other"


def _position_number(value: str) -> int:
    try:
        return int(str(value).split("_")[-1])
    except Exception:
        return 9999


def _site_label(row: pd.Series) -> str:
    n = _position_number(row.common_position)
    if str(row.common_position).startswith("BZD_GAMMA2_"):
        return f"gamma2 {row.residue_name}{n}"
    return f"alpha {row.residue_name}{n}"


def _receptor_label(receptor: str) -> str:
    return "α1" if str(receptor) == "alpha1" else "α2"


def _load(parent: Path) -> pd.DataFrame:
    source = parent / "data" / "raw" / "lineage" / "residue_ligand_frequency_matrix.csv"
    if not source.exists():
        source = parent / "results" / "tables" / "residue_ligand_frequency_matrix.csv"
    d = pd.read_csv(source)
    required = {"receptor", "common_position", "residue_name", "interaction_type", *DRUGS}
    missing = required - set(d.columns)
    if missing:
        raise ValueError(f"missing source columns: {sorted(missing)}")
    d["residue_class"] = d.residue_name.map(_class)
    d["class_rank"] = d.residue_class.map(CATEGORY_RANK)
    d["position_sort"] = d.common_position.map(_position_number)
    d["receptor_rank"] = d.receptor.map({"alpha1": 0, "alpha2": 1}).fillna(9)
    d["interaction_sort"] = d.interaction_type.map(lambda x: INTERACTION_SHORT.get(str(x), str(x)))
    d["site_label"] = d.apply(_site_label, axis=1)
    d["receptor_label"] = d.receptor.map(_receptor_label)
    d["row_label"] = d.apply(lambda r: f"{r.site_label} | {r.receptor_label} | {INTERACTION_SHORT.get(str(r.interaction_type), str(r.interaction_type))}", axis=1)
    d = d.sort_values(["class_rank", "position_sort", "common_position", "residue_name", "receptor_rank", "interaction_sort"], kind="stable").reset_index(drop=True)
    return d


def _mean_center(d: pd.DataFrame) -> pd.DataFrame:
    out = d.copy()
    vals = out[DRUGS].astype(float)
    out[DRUGS] = vals.sub(vals.mean(axis=1), axis=0)
    out["normalization"] = "row_mean_centered: raw frequency minus four-drug row mean"
    return out


def _display_rows(d: pd.DataFrame) -> list[dict]:
    """Insert labelled class rows and unlabelled site separators for plotting."""
    rows: list[dict] = []
    last_class = None
    last_site = None
    for _, r in d.iterrows():
        cls = str(r.residue_class)
        site = (cls, str(r.common_position), str(r.residue_name))
        if cls != last_class:
            rows.append({"kind": "class", "label": cls})
            last_class = cls
            last_site = None
        if site != last_site:
            if last_site is not None:
                rows.append({"kind": "site", "label": ""})
            last_site = site
        rows.append({"kind": "feature", "label": str(r.row_label), "values": [float(r[x]) for x in DRUGS]})
    return rows


def _plot(d: pd.DataFrame, out: Path, normalized: bool) -> None:
    rows = _display_rows(d)
    arr = np.full((len(rows), len(DRUGS)), np.nan)
    labels: list[str] = []
    for i, r in enumerate(rows):
        labels.append(r["label"])
        if r["kind"] == "feature":
            arr[i, :] = r["values"]
    fig_h = max(6.0, 0.34 * len(rows) + 1.8)
    fig, ax = plt.subplots(figsize=(10.4, fig_h))
    cmap = plt.cm.RdBu_r.copy() if normalized else plt.cm.viridis.copy()
    cmap.set_bad("white")
    if normalized:
        vmax = float(np.nanmax(np.abs(arr))) if np.isfinite(arr).any() else 1.0
        vmax = max(vmax, 0.05)
        image = ax.imshow(np.ma.masked_invalid(arr), aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
        color_label = "row mean-centered frequency"
    else:
        image = ax.imshow(np.ma.masked_invalid(arr), aspect="auto", cmap=cmap, vmin=0, vmax=1)
        color_label = "raw interaction frequency"
    ax.set_xticks(range(len(DRUGS)), DRUGS, rotation=30, ha="right")
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ax.tick_params(axis="y", length=0)
    ticklabels = ax.get_yticklabels()
    for i, r in enumerate(rows):
        if r["kind"] == "class":
            ax.axhline(i - 0.5, color="#111111", linewidth=2.0)
            ticklabels[i].set_fontweight("bold")
            ticklabels[i].set_color("#111111")
        elif r["kind"] == "site":
            ax.axhline(i - 0.5, color="#555555", linewidth=0.9)
    ax.set_ylabel("residue class → residue/site → receptor → interaction type")
    ax.set_title("Residue-class grouped BZD-site interaction frequency" if not normalized else "Residue-class grouped row mean-centered pattern")
    fig.colorbar(image, ax=ax, label=color_label, fraction=0.03, pad=0.02)
    fig.subplots_adjust(left=0.44, right=0.92, bottom=0.09, top=0.96)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def _report(run_dir: Path, manifest: dict, d: pd.DataFrame) -> None:
    observed_classes = [x for x in CATEGORIES if x in set(d.residue_class)]
    empty_classes = [x for x in CATEGORIES if x not in set(d.residue_class)]
    header = [
        f"Run ID: {manifest['run_id']}",
        f"Parent Run: {manifest.get('parent_run') or ''}",
        f"Analysis Phase: {manifest['analysis_phase']}",
        f"Started: {manifest['started_at']}",
        f"Completed: {run_manager.iso_now()}", "",
    ]
    text = "\n".join(header + [
        "# RESIDUE_CLASS_GROUPED_HEATMAP_REPORT", "",
        "This run adds a residue-class-first visualization. The previous receptor-first grouped heatmaps and the original frequency heatmap remain unchanged.", "",
        "## Row hierarchy", "",
        "Rows are ordered as residue class → residue/site → receptor → interaction type. Alpha1 and alpha2 are therefore adjacent for the same common-position site; there is no alpha1 block or alpha2 block. Class boundaries use thick separators and site boundaries use medium separators.", "",
        f"Observed classes in the parent matrix: {', '.join(observed_classes)}.",
        f"Classes with no observed candidate feature rows: {', '.join(empty_classes) if empty_classes else 'none'}. Empty classes were not populated with artificial rows.", "",
        "## Two heatmaps", "",
        "The raw heatmap displays the unchanged 0–1 interaction frequency. The row-normalized heatmap displays each feature's raw frequency minus its four-drug mean. The normalized view is for relative pattern reading only; the raw matrix remains the numeric source.", "",
        "## Labels", "",
        "Labels use the compact form `gamma2 TYR58 | α1 | pi-stack` or `alpha HIS102 | α2 | hydrophobic`. Residue class is shown as a heading rather than repeated in every row.", "",
        "## QC", "",
        f"Input feature rows: {len(d)}; raw frequency columns: {', '.join(DRUGS)}; raw values copied without transformation into results/tables/residue_class_grouped_frequency_matrix.csv.", "",
        "## Files", "",
        "results/figures/residue_class_grouped_frequency_heatmap.png; results/figures/residue_class_grouped_row_normalized_heatmap.png; results/tables/residue_class_grouped_frequency_matrix.csv; results/tables/residue_class_grouped_row_normalized_matrix.csv",
    ]) + "\n"
    (run_dir / "reports" / "RESIDUE_CLASS_GROUPED_HEATMAP_REPORT.md").write_text(text)
    (run_dir / "FINAL_REPORT.md").write_text(text)
    (run_dir / "QC_REPORT.md").write_text("\n".join(header + ["# QC_REPORT", "", f"Rows retained: {len(d)}.", "Input raw frequency values were preserved; normalization was used only for the second figure.", ""]))


def run() -> dict:
    parent = run_manager.RUNS / PARENT_RUN
    d = _load(parent)
    run_dir, manifest = run_manager.create_run(
        short_phase_name="residue_class_heatmap",
        analysis_phase="residue-class-first interaction fingerprint visualization",
        analysis_purpose="Reorder residue frequency heatmaps by chemical residue class, then residue/site, receptor, and interaction type.",
        parent_run=PARENT_RUN,
        drugs=DRUGS,
        receptors={"alpha1": "alpha1_beta3_gamma2", "alpha2": "alpha2_beta3_gamma2_local_9CTJ"},
        docking_parameters={"docking_executed": False, "source_run": PARENT_RUN, "row_normalization": "row mean-centered across four drugs"},
        input_files=[f"runs/{PARENT_RUN}/data/raw/lineage/residue_ligand_frequency_matrix.csv"],
        notes=["Existing heatmaps and prior runs are immutable.", "Receptor is secondary within each residue/site; alpha1 and alpha2 are adjacent.", "Residue class is visualization-only and interaction type remains separate."],
        source_run=PARENT_RUN,
        source_files=["data/raw/lineage/residue_ligand_frequency_matrix.csv"],
        code_paths=[Path(__file__), ROOT / "src" / "run_manager.py"],
    )
    lineage = run_dir / "data" / "raw" / "lineage"
    source = parent / "data" / "raw" / "lineage" / "residue_ligand_frequency_matrix.csv"
    run_manager.copy_into_run(source, lineage / "residue_ligand_frequency_matrix.csv")
    manifest["input_files"] = sorted(set(manifest["input_files"] + [str(p.relative_to(run_dir)) for p in lineage.iterdir()]))
    run_manager.write_manifest(run_dir, manifest)
    grouped = d[["receptor", "common_position", "residue_name", "residue_class", "site_label", "receptor_label", "interaction_type", "row_label", *DRUGS, *[f"{x}_status" for x in DRUGS]]].copy()
    grouped.to_csv(run_dir / "results" / "tables" / "residue_class_grouped_frequency_matrix.csv", index=False)
    norm = _mean_center(d)
    normalized = norm[["receptor", "common_position", "residue_name", "residue_class", "site_label", "receptor_label", "interaction_type", "row_label", *DRUGS, "normalization"]].copy()
    normalized.to_csv(run_dir / "results" / "tables" / "residue_class_grouped_row_normalized_matrix.csv", index=False)
    _plot(d, run_dir / "results" / "figures" / "residue_class_grouped_frequency_heatmap.png", normalized=False)
    _plot(norm, run_dir / "results" / "figures" / "residue_class_grouped_row_normalized_heatmap.png", normalized=True)
    _report(run_dir, manifest, d)
    finalized = run_manager.finalize_run(run_dir, manifest, qc_status="COMPLETED_WITH_QC")
    return {"run_id": finalized["run_id"], "run_dir": str(run_dir), "rows": len(d), "observed_classes": sorted(set(d.residue_class), key=lambda x: CATEGORY_RANK[x])}


def main() -> None:
    print(json.dumps(run(), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
