#!/usr/bin/env python3
"""Run the general interaction-fingerprint pipeline in an immutable run dir."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from src import pipeline as pipeline_module
from src import run_manager


ROOT = Path(__file__).resolve().parent


def _resolve(path: str, config_path: Path) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else (ROOT / candidate)


def _copy_inputs(run_dir: Path, cfg: dict, config_path: Path) -> tuple[dict, list[str], list[str]]:
    run_cfg = yaml.safe_load(yaml.safe_dump(cfg))
    raw = run_dir / "raw"
    source_files: list[str] = []
    copied_files: list[str] = []
    for key, value in cfg.get("inputs", {}).items():
        source = _resolve(str(value), config_path)
        source_files.append(str(source.relative_to(ROOT) if source.is_relative_to(ROOT) else source))
        if source.exists() and source.is_file():
            destination = raw / source.name
            run_manager.copy_into_run(source, destination)
            run_cfg["inputs"][key] = str(destination)
            copied_files.append(str(destination.relative_to(run_dir)))
        else:
            # Keep missing paths inside the run namespace so the existing QC
            # code reports them without stopping execution.
            run_cfg["inputs"][key] = str(raw / source.name)

    run_cfg["outputs"] = {
        "processed": str(run_dir / "processed"),
        "figures": str(run_dir / "figures"),
        "tables": str(run_dir / "tables"),
    }
    return run_cfg, source_files, copied_files


def _write_reports(run_dir: Path, manifest: dict, result: dict) -> None:
    completed = run_manager.iso_now()
    header = "\n".join([
        f"Run ID: {manifest['run_id']}",
        f"Parent Run: {manifest.get('parent_run') or ''}",
        f"Analysis Phase: {manifest['analysis_phase']}",
        f"Started: {manifest['started_at']}",
        f"Completed: {completed}",
        "",
    ])
    qc_rows = result.get("qc", []) if isinstance(result, dict) else []
    qc_text = header + "# QC_REPORT\n\n" + json.dumps(qc_rows, ensure_ascii=False, indent=2, default=str) + "\n"
    final_text = header + "# FINAL_REPORT\n\n" + (
        "General interaction-fingerprint pipeline completed in a versioned run directory.\n\n"
        f"QC rows: {len(qc_rows)}. Missing inputs remain explicit in tables/qc_summary.csv.\n"
        "No input data were synthesized by the run wrapper.\n"
    )
    (run_dir / "QC_REPORT.md").write_text(qc_text)
    (run_dir / "FINAL_REPORT.md").write_text(final_text)


def run(config_path: Path, phase: str) -> dict:
    cfg = yaml.safe_load(config_path.read_text()) or {}
    latest = run_manager.latest_run()
    # A parent denotes inherited run data, not merely the most recently
    # completed unrelated phase.  Configurations that point to a run namespace
    # are therefore linked; legacy workspace inputs retain a legacy parent.
    configured_inputs = [str(v) for v in cfg.get("inputs", {}).values()]
    parent = latest if latest and any("runs/" in value or "/runs/" in value for value in configured_inputs) else "legacy_unversioned_pipeline"
    run_dir, manifest = run_manager.create_run(
        short_phase_name=phase,
        analysis_phase=cfg.get("analysis", {}).get("phase", phase),
        analysis_purpose="General interaction-fingerprint pipeline with immutable run-level outputs and input QC.",
        parent_run=parent,
        drugs=[],
        receptors={},
        docking_parameters={"docking_executed": False, "config_analysis": cfg.get("analysis", {})},
        input_files=[str(v) for v in cfg.get("inputs", {}).values()],
        notes=["The requested run-directory policy is enforced by this wrapper.", "Existing legacy outputs are not modified."],
        source_run=parent,
        source_files=[str(v) for v in cfg.get("inputs", {}).values()],
        code_paths=[Path(__file__), ROOT / "src" / "pipeline.py", ROOT / "src" / "run_manager.py"],
    )
    run_manager.copy_into_run(config_path, run_dir / "config" / config_path.name)
    run_cfg, source_files, copied_files = _copy_inputs(run_dir, cfg, config_path)
    manifest["input_files"] = sorted(set(source_files + copied_files + [f"config/{config_path.name}"]))
    run_manager.write_manifest(run_dir, manifest)
    try:
        result = pipeline_module.run(run_cfg)
        _write_reports(run_dir, manifest, result)
        finalized = run_manager.finalize_run(run_dir, manifest, qc_status="COMPLETED_WITH_QC")
    except Exception as exc:
        (run_dir / "QC_REPORT.md").write_text(
            f"Run ID: {manifest['run_id']}\nParent Run: {manifest.get('parent_run') or ''}\n"
            f"Analysis Phase: {manifest['analysis_phase']}\nStarted: {manifest['started_at']}\n"
            f"Completed: {run_manager.iso_now()}\n\n# QC_REPORT\n\nPipeline error: {exc}\n"
        )
        run_manager.finalize_run(run_dir, manifest, qc_status="FAILED")
        raise
    result = dict(result or {})
    result.update({"run_id": finalized["run_id"], "run_dir": str(run_dir)})
    print(json.dumps(result, ensure_ascii=False, indent=2, default=str))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/poc.yaml")
    parser.add_argument("--phase", default="pipeline", help="short phase name used in run_id")
    args = parser.parse_args()
    run(Path(args.config).resolve(), args.phase)
