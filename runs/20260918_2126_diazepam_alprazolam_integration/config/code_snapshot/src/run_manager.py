"""Immutable run-directory management for analysis phases.

The manager is intentionally small and filesystem based.  It creates a new
run directory for every invocation, records lineage and code provenance, and
never mutates a previous run.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "runs"


def now_local() -> datetime:
    return datetime.now().astimezone()


def iso_now() -> str:
    return now_local().isoformat(timespec="seconds")


def git_commit() -> str | None:
    """Return the current commit when this workspace is under git."""
    try:
        value = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=3,
        ).strip()
        return value or None
    except Exception:
        return None


def code_version(paths: Iterable[Path]) -> str:
    """Create a deterministic short fingerprint of the analysis code."""
    digest = hashlib.sha256()
    for path in sorted((Path(p) for p in paths), key=str):
        if path.exists() and path.is_file():
            digest.update(str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path).encode())
            digest.update(path.read_bytes())
    return f"workspace_code_sha256:{digest.hexdigest()[:16]}"


def latest_run() -> str | None:
    pointer = RUNS / "LATEST_RUN.txt"
    if not pointer.exists():
        return None
    value = pointer.read_text().strip()
    return value or None


def new_run_id(short_phase_name: str, when: datetime | None = None) -> str:
    """Generate a local-time run id, avoiding collisions without overwriting."""
    when = when or now_local()
    base = f"{when:%Y%m%d_%H%M}_{short_phase_name}"
    RUNS.mkdir(parents=True, exist_ok=True)
    candidate = base
    suffix = 2
    while (RUNS / candidate).exists():
        candidate = f"{base}_{suffix}"
        suffix += 1
    return candidate


def create_run(
    *,
    short_phase_name: str,
    analysis_phase: str,
    analysis_purpose: str,
    parent_run: str | None,
    drugs: list[str],
    receptors: Mapping[str, str],
    docking_parameters: Mapping[str, Any],
    input_files: list[str],
    notes: list[str],
    source_run: str | None = None,
    source_files: list[str] | None = None,
    code_paths: Iterable[Path] = (),
) -> tuple[Path, dict[str, Any]]:
    run_id = new_run_id(short_phase_name)
    run_dir = RUNS / run_id
    # Keep the requested top-level directories and conventional nested aliases
    # available for future phases without placing duplicate result files.
    for name in ("config", "data", "raw", "processed", "results", "tables", "figures", "logs", "reports"):
        (run_dir / name).mkdir(parents=True, exist_ok=True)
    for name in ("raw", "processed"):
        (run_dir / "data" / name).mkdir(parents=True, exist_ok=True)
    for name in ("tables", "figures"):
        (run_dir / "results" / name).mkdir(parents=True, exist_ok=True)

    manifest: dict[str, Any] = {
        "run_id": run_id,
        "started_at": iso_now(),
        "completed_at": None,
        "analysis_phase": analysis_phase,
        "analysis_purpose": analysis_purpose,
        "parent_run": parent_run,
        "drugs": drugs,
        "receptors": dict(receptors),
        "docking_parameters": dict(docking_parameters),
        "input_files": list(input_files),
        "output_files": [],
        "git_commit": git_commit(),
        "code_version": code_version(code_paths),
        "notes": list(notes),
        "qc_status": "RUNNING",
    }
    if source_run is not None:
        manifest["source_run"] = source_run
    if source_files is not None:
        manifest["source_files"] = list(source_files)
    write_manifest(run_dir, manifest)
    return run_dir, manifest


def write_manifest(run_dir: Path, manifest: Mapping[str, Any]) -> None:
    (run_dir / "run_manifest.json").write_text(
        json.dumps(dict(manifest), ensure_ascii=False, indent=2, sort_keys=False) + "\n"
    )


def copy_into_run(source: Path, destination: Path) -> str:
    """Copy one input into a run, preserving bytes and returning run-relative path."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return str(destination)


def finalize_run(
    run_dir: Path,
    manifest: dict[str, Any],
    *,
    qc_status: str,
    extra_notes: list[str] | None = None,
) -> dict[str, Any]:
    manifest["completed_at"] = iso_now()
    manifest["qc_status"] = qc_status
    if extra_notes:
        manifest.setdefault("notes", []).extend(extra_notes)
    manifest["output_files"] = sorted(
        str(path.relative_to(run_dir))
        for path in run_dir.rglob("*")
        if path.is_file() and path.name != "run_manifest.json"
    )
    write_manifest(run_dir, manifest)
    RUNS.mkdir(parents=True, exist_ok=True)
    (RUNS / "LATEST_RUN.txt").write_text(f"{manifest['run_id']}\n")
    return manifest

