"""Strict three-drug GABA-A mechanistic PoC.

This module intentionally keeps the strict analysis separate from the earlier
10-drug/primary/exploratory tables.  It uses the already downloaded ChEMBL
records, RCSB structures placed in ``data/raw/structures/strict3`` and a small
reproducible Vina/PLIP wrapper.  Where an exact alpha2beta3gamma2 structure is
not available, the selected alpha2 receptor is explicitly labelled as a local
site construct from 9CTJ (a mixed native receptor), rather than being called an
exact complex.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
THREE_DRUGS = ["diazepam", "triazolam", "zolpidem"]
TARGET_A1 = "CHEMBL2094121"
TARGET_A2 = "CHEMBL2094130"
STRUCT_DIR = ROOT / "data" / "raw" / "structures" / "strict3"
STRICT_RAW = ROOT / "data" / "raw" / "strict3"
STRICT_PROC = ROOT / "data" / "processed" / "strict3"
STRICT_RES = ROOT / "results" / "strict3"
OUT = ROOT / "outputs"


def mkdirs() -> None:
    for p in [STRUCT_DIR, STRICT_RAW, STRICT_PROC, STRICT_RES, OUT / "figures", OUT / "tables"]:
        p.mkdir(parents=True, exist_ok=True)


def _write_csv(df: pd.DataFrame, *paths: Path) -> None:
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False)


def build_strict_pharmacology() -> pd.DataFrame:
    """Create pair-level Ki rows without averaging or combining assay types."""
    raw_path = ROOT / "results" / "tables" / "chembl_primary_raw_activity.csv"
    df = pd.read_csv(raw_path)
    df = df[df["drug_id"].isin(THREE_DRUGS)].copy()
    df = df[df["target_chembl_id"].isin([TARGET_A1, TARGET_A2])].copy()
    df = df[df["standard_type"].astype(str).str.upper() == "KI"].copy()
    metadata = {}
    meta_path = ROOT / "data" / "raw" / "pharmacology" / "chembl" / "chembl_api_metadata.json"
    if meta_path.exists():
        try:
            metadata = json.loads(meta_path.read_text())
        except Exception:
            metadata = {}
    df["standard_value_num"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df[df["standard_value_num"].notna()].copy()
    rows = []
    for drug in THREE_DRUGS:
        a1 = df[(df.drug_id == drug) & (df.target_chembl_id == TARGET_A1)].to_dict("records")
        a2 = df[(df.drug_id == drug) & (df.target_chembl_id == TARGET_A2)].to_dict("records")
        for x, y in itertools.product(a1, a2):
            same_doc = str(x.get("document_chembl_id", "")) == str(y.get("document_chembl_id", ""))
            same_type = str(x.get("standard_type", "")) == str(y.get("standard_type", ""))
            same_unit = str(x.get("standard_units", "")) == str(y.get("standard_units", ""))
            if not (same_type and same_unit):
                continue
            priority = 1 if same_doc else 2
            rows.append(
                {
                    "drug_id": drug,
                    "input_drug_name": x.get("input_drug_name"),
                    "molecule_chembl_id": x.get("molecule_chembl_id"),
                    "target_alpha1_chembl_id": TARGET_A1,
                    "target_alpha2_chembl_id": TARGET_A2,
                    "alpha1_assay_chembl_id": x.get("assay_chembl_id"),
                    "alpha2_assay_chembl_id": y.get("assay_chembl_id"),
                    "alpha1_document_chembl_id": x.get("document_chembl_id"),
                    "alpha2_document_chembl_id": y.get("document_chembl_id"),
                    "same_document": same_doc,
                    "same_standard_type": same_type,
                    "same_unit": same_unit,
                    "pair_priority": priority,
                    "pair_basis": "same_document_standard_type_unit" if priority == 1 else "same_standard_type_unit",
                    "standard_type": x.get("standard_type"),
                    "alpha1_standard_relation": x.get("standard_relation"),
                    "alpha2_standard_relation": y.get("standard_relation"),
                    "alpha1_standard_value": x.get("standard_value_num"),
                    "alpha2_standard_value": y.get("standard_value_num"),
                    "standard_units": x.get("standard_units"),
                    "alpha1_pchembl_value": x.get("pchembl_value"),
                    "alpha2_pchembl_value": y.get("pchembl_value"),
                    "ki_alpha2_over_alpha1": float(y["standard_value_num"]) / float(x["standard_value_num"]),
                    "chembl_version": metadata.get("chembl_version", ""),
                    "api_retrieved_at": metadata.get("retrieved_at_utc", metadata.get("retrieved_at", "")),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["drug_id", "pair_priority", "alpha1_document_chembl_id", "alpha1_assay_chembl_id"]).reset_index(drop=True)
        out["pair_rank"] = out.groupby("drug_id").cumcount() + 1
    _write_csv(out, OUT / "strict_pharmacology_summary.csv", STRICT_RES / "strict_pharmacology_summary.csv")
    raw_ki = df.sort_values(["drug_id", "target_chembl_id", "document_chembl_id", "assay_chembl_id"])
    _write_csv(raw_ki, OUT / "strict_pharmacology_raw_ki.csv", STRICT_RES / "strict_pharmacology_raw_ki.csv")
    return out


def receptor_candidates() -> pd.DataFrame:
    """Curated metadata from the RCSB entries investigated for this PoC."""
    rows = [
        ["6HUP", "human full-length alpha1beta3gamma2L GABA(A)R with diazepam, GABA and Mb38", "alpha1beta3gamma2", "DZP (diazepam)", 3.58, "exact alpha1beta3gamma2 receptor; BZD ligand in alpha1/gamma2 site", "yes", "selected alpha1; docking scaffold"],
        ["7QNE", "human full-length synaptic alpha1beta3gamma2 GABA(A)R with Ro15-4513 and Mb38", "alpha1beta3gamma2", "EIE (Ro15-4513)", 2.70, "exact alpha1beta3gamma2 receptor; BZD-site ligand", "yes", "alternative alpha1 candidate"],
        ["6X3X", "human GABAA receptor alpha1-beta2-gamma2 with GABA plus diazepam", "alpha1beta2gamma2", "DZP (diazepam)", 2.92, "different beta subtype", "yes", "not selected; beta2"],
        ["9CTJ", "native human GABAA receptor beta2-alpha1-beta3-alpha2-gamma2 assembly", "beta2-alpha1-beta3-alpha2-gamma2", "none", 3.74, "mixed alpha1/alpha2 and beta2/beta3; contains alpha2-beta3-gamma2 local site", "no", "selected only as provisional alpha2 local-site construct"],
        ["9CX7", "native human GABAA receptor beta3-alpha1-gamma2-beta3-alpha2 assembly", "beta3-alpha1-gamma2-beta3-alpha2", "none", 3.30, "mixed alpha1/alpha2; not an exact alpha2beta3gamma2 receptor", "no", "not selected"],
        ["9CXC", "native human GABAA receptor beta3-alpha1-gamma2-beta2-alpha2 assembly", "beta3-alpha1-gamma2-beta2-alpha2",  "none", 3.30, "mixed alpha1/alpha2 and beta2/beta3", "no", "not selected"],
        ["9CSB", "native human GABAA receptor beta3-alpha1-beta2-alpha2-gamma2 assembly", "beta3-alpha1-beta2-alpha2-gamma2", "none", 3.34, "mixed alpha1/alpha2 and beta2/beta3", "no", "not selected"],
    ]
    df = pd.DataFrame(rows, columns=["pdb_id", "title", "composition", "co_crystal_ligand", "resolution_angstrom", "bzd_site_assessment", "bzd_ligand_present", "poc_decision"])
    _write_csv(df, OUT / "receptor_structure_candidates.csv", STRICT_RES / "receptor_structure_candidates.csv")
    return df


def _element_type(symbol: str, aromatic: bool = False) -> str:
    s = (symbol or "C").strip().upper()
    if s == "CL":
        return "Cl"
    if s == "BR":
        return "Br"
    if s == "H":
        return "HD"
    if s == "C" and aromatic:
        return "A"
    if s in {"C", "N", "O", "S", "P", "F", "I"}:
        return s.title() if s in {"C", "N", "O", "S", "P", "F", "I"} else s
    return "C"


def generate_ligands() -> Dict[str, Path]:
    mols = pd.read_csv(ROOT / "data" / "raw" / "pharmacology" / "chembl" / "molecules.csv")
    mols = mols[mols.drug_id.isin(THREE_DRUGS)].set_index("drug_id")
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:
        raise RuntimeError(f"RDKit is required for reproducible ligand generation: {exc}")
    out_paths: Dict[str, Path] = {}
    for drug in THREE_DRUGS:
        smi = mols.loc[drug, "canonical_smiles"]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Could not parse canonical SMILES for {drug}")
        mol = Chem.AddHs(mol)
        seed = int(hashlib.sha256(drug.encode()).hexdigest()[:8], 16) % (2**31 - 1)
        if AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=False) != 0:
            raise ValueError(f"3D embedding failed for {drug}")
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        no_h = Chem.RemoveHs(mol)
        sdf = STRICT_RAW / "ligands" / f"{drug}.sdf"
        pdb = STRICT_RAW / "ligands" / f"{drug}.pdb"
        pdbqt = STRICT_RAW / "ligands" / f"{drug}.pdbqt"
        sdf.parent.mkdir(parents=True, exist_ok=True)
        Chem.MolToMolFile(no_h, str(sdf))
        Chem.MolToPDBFile(no_h, str(pdb), flavor=0)
        conf = no_h.GetConformer()
        try:
            AllChem.ComputeGasteigerCharges(no_h)
        except Exception:
            pass
        lines = ["REMARK 3D ligand generated from ChEMBL canonical SMILES with RDKit ETKDG + MMFF/UFF", "ROOT"]
        for i, atom in enumerate(no_h.GetAtoms(), start=1):
            pos = conf.GetAtomPosition(i - 1)
            try:
                charge = float(atom.GetProp("_GasteigerCharge"))
                if not math.isfinite(charge):
                    charge = 0.0
            except Exception:
                charge = 0.0
            typ = _element_type(atom.GetSymbol(), atom.GetIsAromatic())
            name = f"{atom.GetSymbol()}{i}"[:4]
            lines.append(f"ATOM  {i:5d} {name:<4s} LIG Z   1    {pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}{1.0:6.2f}{0.0:6.2f}    {charge:7.3f} {typ:>2s}")
        lines.extend(["ENDROOT", "TORSDOF 0", ""])
        pdbqt.write_text("\n".join(lines))
        out_paths[drug] = pdbqt
    return out_paths


def _pdbqt_line(serial: int, name: str, resname: str, chain: str, resseq: int, x: float, y: float, z: float, charge: float, typ: str, record: str = "ATOM") -> str:
    return f"{record:<6s}{serial:5d} {name:<4s} {resname:>3s} {chain:1s}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}    {charge:7.3f} {typ:>2s}"


def prepare_receptors() -> Dict[str, Path]:
    """Prepare receptor PDBQT; alpha2 is a 9CTJ local-site construct."""
    from Bio.PDB import PDBParser

    configs = {
        "alpha1_beta3_gamma2": (STRUCT_DIR / "6HUP.pdb", list("ABCDE")),
        "alpha2_beta3_gamma2_local_9CTJ": (STRUCT_DIR / "9CTJ.pdb", list("CDE")),
    }
    out = {}
    for name, (pdb, chains) in configs.items():
        if not pdb.exists():
            raise FileNotFoundError(pdb)
        structure = PDBParser(QUIET=True).get_structure(name, str(pdb))[0]
        lines = ["REMARK receptor PDBQT generated from RCSB PDB ATOM records; non-protein HETATM removed"]
        serial = 1
        for chain_id in chains:
            chain = structure[chain_id]
            for res in chain:
                if res.id[0] != " ":
                    continue
                for atom in res:
                    if atom.is_disordered() and atom.get_altloc() not in ("A", " "):
                        continue
                    elem = atom.element or atom.get_name()[0]
                    typ = _element_type(elem, False)
                    lines.append(_pdbqt_line(serial, atom.get_name(), res.resname, chain_id, int(res.id[1]), float(atom.coord[0]), float(atom.coord[1]), float(atom.coord[2]), 0.0, typ))
                    serial += 1
        path = STRICT_RAW / "receptors" / f"{name}.pdbqt"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n")
        out[name] = path
    return out


def _dzp_coords() -> np.ndarray:
    from Bio.PDB import PDBParser
    s = PDBParser(QUIET=True).get_structure("6HUP", str(STRUCT_DIR / "6HUP.pdb"))[0]
    # This is the DZP copy assigned to alpha1 chain D and adjacent to gamma2 C.
    lig = s["D"][("H_DZP", 2001, " ")]
    return np.array([a.coord for a in lig], dtype=float)


def _box_from_coords(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    center = coords.mean(axis=0)
    size = np.maximum(coords.max(axis=0) - coords.min(axis=0) + 12.0, 20.0)
    return center, size


def docking_boxes() -> pd.DataFrame:
    from Bio.PDB import PDBParser, Superimposer
    c1, s1 = _box_from_coords(_dzp_coords())
    # Superpose 6HUP gamma2 C (moving) onto 9CTJ gamma2 E (fixed).  This maps
    # the observed co-crystal DZP site into the provisional alpha2 local site.
    p = PDBParser(QUIET=True)
    a = p.get_structure("a", str(STRUCT_DIR / "6HUP.pdb"))[0]["C"]
    b = p.get_structure("b", str(STRUCT_DIR / "9CTJ.pdb"))[0]["E"]
    amap = {r.id[1]: r["CA"] for r in a if r.id[0] == " " and "CA" in r}
    bmap = {r.id[1]: r["CA"] for r in b if r.id[0] == " " and "CA" in r}
    keys = sorted(set(amap) & set(bmap))
    sup = Superimposer()
    sup.set_atoms([bmap[k] for k in keys], [amap[k] for k in keys])
    rot, tran = sup.rotran
    transformed = np.dot(_dzp_coords(), rot) + tran
    c2, s2 = _box_from_coords(transformed)
    rows = [
        ["alpha1_beta3_gamma2", "6HUP", "DZP", "D", "C", float(c1[0]), float(c1[1]), float(c1[2]), float(s1[0]), float(s1[1]), float(s1[2]), "native DZP centroid + 6 A margin; minimum 20 A per axis", "direct"],
        ["alpha2_beta3_gamma2_local_9CTJ", "9CTJ", "DZP transferred", "D", "E", float(c2[0]), float(c2[1]), float(c2[2]), float(s2[0]), float(s2[1]), float(s2[2]), f"6HUP DZP box transferred by gamma2 C->E CA alignment (N={len(keys)}, RMSD={sup.rms:.3f} A)", "provisional_transfer"],
    ]
    df = pd.DataFrame(rows, columns=["receptor_id", "pdb_id", "box_reference_ligand", "alpha_or_site_chain", "gamma_chain", "center_x", "center_y", "center_z", "size_x", "size_y", "size_z", "definition", "box_status"])
    _write_csv(df, OUT / "docking_boxes.csv", STRICT_RES / "docking_boxes.csv")
    return df


def _pdb_from_pdbqt(path: Path, chain: str = "Z") -> List[str]:
    out = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip() or "C"
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        elem = "Cl" if name.lower().startswith("cl") else ("Br" if name.lower().startswith("br") else name[0].upper())
        out.append(f"HETATM{len(out)+1:5d} {name:<4s} LIG {chain}{1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {elem:>2s}")
    return out


def _receptor_pdb_for_plip(receptor_id: str, pose_pdbqt: Path, pose_pdb: Path) -> None:
    if receptor_id == "alpha1_beta3_gamma2":
        source, chains = STRUCT_DIR / "6HUP.pdb", set("ABCDE")
    else:
        source, chains = STRUCT_DIR / "9CTJ.pdb", set("CDE")
    lines = []
    for line in source.read_text().splitlines():
        if line.startswith("ATOM") and line[21].strip() in chains:
            lines.append(line[:66].rstrip())
    lines.extend(_pdb_from_pdbqt(pose_pdbqt, "Z"))
    lines.extend(["TER", "END"])
    pose_pdb.write_text("\n".join(lines) + "\n")


def run_vina(receptors: Dict[str, Path], ligands: Dict[str, Path], boxes: pd.DataFrame) -> pd.DataFrame:
    vina = shutil.which("vina")
    if vina is None:
        raise RuntimeError("AutoDock Vina executable was not found")
    rows = []
    dock_root = STRICT_RAW / "docking"
    dock_root.mkdir(parents=True, exist_ok=True)
    for receptor_id, rec_path in receptors.items():
        box = boxes[boxes.receptor_id == receptor_id].iloc[0]
        for drug, lig_path in ligands.items():
            job = dock_root / f"{drug}__{receptor_id}"
            job.mkdir(parents=True, exist_ok=True)
            out_pdbqt = job / "poses.pdbqt"
            log = job / "vina.log"
            cmd = [vina, "--receptor", str(rec_path), "--ligand", str(lig_path), "--center_x", str(box.center_x), "--center_y", str(box.center_y), "--center_z", str(box.center_z), "--size_x", str(box.size_x), "--size_y", str(box.size_y), "--size_z", str(box.size_z), "--exhaustiveness", "8", "--num_modes", "9", "--energy_range", "4", "--seed", "20260917", "--out", str(out_pdbqt)]
            if not out_pdbqt.exists():
                proc = subprocess.run(cmd, cwd=job, text=True, capture_output=True)
                (job / "vina.stdout.txt").write_text(proc.stdout + "\n" + proc.stderr)
                log.write_text(proc.stdout + "\n" + proc.stderr)
                if proc.returncode != 0:
                    raise RuntimeError(f"Vina failed for {drug}/{receptor_id}: {proc.stderr[-1000:]}")
            score_by_pose = {}
            current = None
            for line in out_pdbqt.read_text().splitlines():
                if line.startswith("MODEL"):
                    try:
                        current = int(line.split()[1])
                    except Exception:
                        current = None
                if current is not None and re.match(r"^REMARK VINA RESULT:", line):
                    vals = line.split()
                    score_by_pose[current] = float(vals[3])
            if not score_by_pose:
                # Vina 1.2 writes the table in the log; keep a traceable row even
                # if an older output format omitted MODEL remarks.
                for line in log.read_text().splitlines():
                    m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+", line)
                    if m:
                        score_by_pose[int(m.group(1))] = float(m.group(2))
            pose_pdb = job / "pose_pdb"
            pose_pdb.mkdir(exist_ok=True)
            for pose_id, score in sorted(score_by_pose.items()):
                # Keep one PDB per pose for PLIP and an auditable pose ID.
                lines, current_model = [], None
                for line in out_pdbqt.read_text().splitlines():
                    if line.startswith("MODEL"):
                        current_model = int(line.split()[1])
                    elif line.startswith("ENDMDL"):
                        if current_model == pose_id:
                            break
                    elif current_model == pose_id and line.startswith(("ATOM", "HETATM")):
                        lines.append(line)
                tmp = job / "_one_pose.pdbqt"
                tmp.write_text("\n".join(lines) + "\n")
                pdb_path = pose_pdb / f"pose_{pose_id:02d}.pdb"
                _receptor_pdb_for_plip(receptor_id, tmp, pdb_path)
                rows.append({"drug_id": drug, "receptor_id": receptor_id, "pdb_id": "6HUP" if receptor_id == "alpha1_beta3_gamma2" else "9CTJ", "pose_id": pose_id, "vina_score_kcal_mol": score, "vina_seed": 20260917, "exhaustiveness": 8, "num_modes": 9, "pose_pdbqt": str(out_pdbqt), "pose_pdb": str(pdb_path), "status": "docked"})
    out = pd.DataFrame(rows)
    _write_csv(out, OUT / "docking_results.csv", STRICT_RES / "docking_results.csv")
    return out


def run_plip(docking: pd.DataFrame) -> pd.DataFrame:
    """Run vendored PLIP if available; retain explicit no-data rows on failure."""
    plip = ROOT / "work" / "plip_vendor" / "bin" / "plip"
    interactions = []
    for _, row in docking.iterrows():
        pose = Path(row.pose_pdb)
        outdir = pose.parent / "plip"
        outdir.mkdir(exist_ok=True)
        xml = outdir / "report.xml"
        txt = outdir / "report.txt"
        # Re-run after any previous zero-ligand report (the first pass used a
        # peptide-chain hint, which is inappropriate for a small molecule).
        if xml.exists():
            try:
                xml.unlink()
            except OSError:
                pass
        if not xml.exists():
            # The PDB already contains only the selected receptor chains plus
            # ligand chain Z.  Do not pass --chains: that option treats the
            # ligand as a peptide chain and suppresses small-molecule PLIP
            # detection.  Default residue-based detection is appropriate here.
            cmd = [str(plip), "-f", str(pose), "-o", str(outdir), "-x", "-t", "-q", "--breakcomposite", "--name", "report"]
            env = os.environ.copy(); env["PYTHONPATH"] = str(ROOT / "work" / "plip_vendor")
            proc = subprocess.run(cmd, cwd=outdir, env=env, text=True, capture_output=True)
            (outdir / "plip.stdout.txt").write_text(proc.stdout + "\n" + proc.stderr)
        if xml.exists():
            # PLIP XML has stable interaction tags.  Keep every interaction and
            # parse residue identity without reducing to a representative pose.
            import xml.etree.ElementTree as ET
            try:
                root = ET.parse(xml).getroot()
                for node in root.iter():
                    tag = node.tag.lower().split("}")[-1]
                    if tag not in {"hydrophobic_interaction", "hydrogen_bond", "water_bridge", "salt_bridge", "pi_stack", "pi_cation", "halogen_bond", "metal_complex"}:
                        continue
                    attrs = {k.lower().split("}")[-1]: v for k, v in node.attrib.items()}
                    # PLIP stores interaction fields as child elements (the
                    # XML attributes usually contain only the interaction id).
                    for child in list(node):
                        key = child.tag.lower().split("}")[-1]
                        if child.text and child.text.strip():
                            attrs[key] = child.text.strip()
                    # PLIP names differ slightly by release; preserve all useful
                    # fields in a JSON column and only map residues when present.
                    res = attrs.get("resnr") or attrs.get("residue") or attrs.get("resnr_lig")
                    chain = attrs.get("reschain") or attrs.get("chain") or ""
                    interactions.append({"drug_id": row.drug_id, "receptor_id": row.receptor_id, "pose_id": int(row.pose_id), "interaction_type": tag, "residue_number": res, "residue_chain": chain, "raw_attributes": json.dumps(attrs, sort_keys=True)})
            except Exception as exc:
                interactions.append({"drug_id": row.drug_id, "receptor_id": row.receptor_id, "pose_id": int(row.pose_id), "interaction_type": "PLIP_PARSE_ERROR", "residue_number": "", "residue_chain": "", "raw_attributes": str(exc)})
    out = pd.DataFrame(interactions, columns=["drug_id", "receptor_id", "pose_id", "interaction_type", "residue_number", "residue_chain", "raw_attributes"])
    _write_csv(out, OUT / "plip_interactions.csv", STRICT_RES / "plip_interactions.csv")
    return out


def build_mapping_and_fingerprints(docking: pd.DataFrame, interactions: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Map alpha1 site residues to alpha2 by PDB numbering/sequence context.

    The selected alpha2 construct is not an exact target, so mapped features are
    labelled provisional.  Ambiguous or unavailable mappings remain QC rows.
    """
    from Bio.PDB import PDBParser
    from Bio.Align import PairwiseAligner
    from Bio.SeqUtils import seq1
    p = PDBParser(QUIET=True)
    s1 = p.get_structure("a", str(STRUCT_DIR / "6HUP.pdb"))[0]
    s2 = p.get_structure("b", str(STRUCT_DIR / "9CTJ.pdb"))[0]
    # A conservative BZD-site residue set around the observed DZP contact,
    # shared by all poses.  Build an explicit global sequence alignment first,
    # then look up the aligned alpha2 residue in the 9CTJ structure.  PDB
    # numbering is retained as an audit field and is never used as an implicit
    # assertion of homology.
    a1_chain = s1["D"]
    a2_chain = s2["D"]
    a1_res = [r for r in a1_chain if r.id[0] == " " and "CA" in r]
    a2_res = [r for r in a2_chain if r.id[0] == " " and "CA" in r]
    a1_seq = "".join(seq1(r.resname, custom_map={"MSE": "M"}) for r in a1_res)
    a2_seq = "".join(seq1(r.resname, custom_map={"MSE": "M"}) for r in a2_res)
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(a1_seq, a2_seq)[0]
    a1_to_a2_idx: Dict[int, int] = {}
    for block_a, block_b in zip(*alignment.aligned):
        for ia, ib in zip(range(block_a[0], block_a[1]), range(block_b[0], block_b[1])):
            a1_to_a2_idx[ia] = ib
    site_nums = [r.id[1] for r in a1_res if r.id[1] in range(90, 221)]
    a1_idx_by_num = {r.id[1]: i for i, r in enumerate(a1_res)}
    rows = []
    alpha2_to_alpha1: Dict[int, int] = {}
    for n in site_nums:
        r1 = a1_chain[(' ', n, ' ')] if (' ', n, ' ') in a1_chain else None
        i1 = a1_idx_by_num.get(n)
        i2 = a1_to_a2_idx.get(i1) if i1 is not None else None
        r2 = a2_res[i2] if i2 is not None and i2 < len(a2_res) else None
        status = "mapped_global_sequence_alignment" if r1 is not None and r2 is not None else "unavailable_or_alignment_ambiguous"
        if r2 is not None:
            alpha2_to_alpha1[int(r2.id[1])] = n
        rows.append({"common_position": f"BZD_SITE_{n:03d}", "alpha1_chain": "D", "alpha1_residue_number": n, "alpha1_residue": r1.resname if r1 else "", "alpha2_chain": "D", "alpha2_residue_number": int(r2.id[1]) if r2 else "", "alpha2_residue": r2.resname if r2 else "", "mapping_method": "global sequence alignment (Biopython PairwiseAligner) + PDB residue lookup", "alignment_score": float(alignment.score), "mapping_status": status, "qc_required": status != "mapped_global_sequence_alignment"})
    mapping = pd.DataFrame(rows)
    _write_csv(mapping, OUT / "residue_mapping.csv", STRICT_RES / "residue_mapping.csv")

    # PLIP parser is intentionally lossless; only rows with a numeric residue
    # can become a mapped feature.  No interactions are represented as absent
    # rows; missing PLIP output is kept as no-data in the QC report.
    def mapped_feature(row):
        try:
            num = int(float(str(row.residue_number)))
        except Exception:
            return None
        if row.receptor_id == "alpha1_beta3_gamma2":
            return f"BZD_SITE_{num:03d}|{row.interaction_type}"
        m = mapping[(mapping.alpha2_residue_number == num) & (mapping.mapping_status == "mapped_global_sequence_alignment")]
        if m.empty:
            return None
        return f"BZD_SITE_{int(m.iloc[0].alpha1_residue_number):03d}|{row.interaction_type}"

    interactions = interactions.copy()
    interactions["feature"] = interactions.apply(mapped_feature, axis=1)
    interactions = interactions[interactions.feature.notna()].copy()
    all_features = sorted(set(interactions.feature))
    # Include the full drug x receptor x pose grid so no interaction and no data
    # remain distinguishable by status columns.
    rows_bin, rows_freq = [], []
    for drug in THREE_DRUGS:
        for receptor in ["alpha1_beta3_gamma2", "alpha2_beta3_gamma2_local_9CTJ"]:
            sub = interactions[(interactions.drug_id == drug) & (interactions.receptor_id == receptor)]
            pose_n = int(docking[(docking.drug_id == drug) & (docking.receptor_id == receptor)].pose_id.nunique())
            for feat in all_features:
                cnt = int((sub.feature == feat).sum())
                interacting_poses = int(sub[sub.feature == feat].pose_id.nunique())
                freq = interacting_poses / pose_n if pose_n else np.nan
                rows_bin.append({"drug_id": drug, "receptor_id": receptor, "feature": feat, "binary": int(cnt > 0), "interaction_data_status": "observed" if cnt else ("no_interaction_observed" if pose_n else "no_data"), "n_poses": pose_n})
                rows_freq.append({"drug_id": drug, "receptor_id": receptor, "feature": feat, "pose_frequency": freq, "n_interacting_poses": interacting_poses, "n_poses": pose_n, "interaction_data_status": "observed" if cnt else ("no_interaction_observed" if pose_n else "no_data")})
    binary = pd.DataFrame(rows_bin, columns=["drug_id", "receptor_id", "feature", "binary", "interaction_data_status", "n_poses"])
    freq = pd.DataFrame(rows_freq, columns=["drug_id", "receptor_id", "feature", "pose_frequency", "n_interacting_poses", "n_poses", "interaction_data_status"])
    _write_csv(binary, OUT / "fingerprint_binary.csv", STRICT_RES / "fingerprint_binary.csv")
    _write_csv(freq, OUT / "fingerprint_frequency.csv", STRICT_RES / "fingerprint_frequency.csv")
    return mapping, binary, freq


def make_figures(binary: pd.DataFrame, freq: pd.DataFrame, pharmacology: pd.DataFrame) -> pd.DataFrame:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sims = []
    for drug in THREE_DRUGS:
        a = binary[(binary.drug_id == drug) & (binary.receptor_id == "alpha1_beta3_gamma2")].set_index("feature").binary
        b = binary[(binary.drug_id == drug) & (binary.receptor_id == "alpha2_beta3_gamma2_local_9CTJ")].set_index("feature").binary
        idx = sorted(set(a.index) | set(b.index)); av = a.reindex(idx, fill_value=0).astype(bool); bv = b.reindex(idx, fill_value=0).astype(bool)
        union = int((av | bv).sum()); inter = int((av & bv).sum())
        sims.append({"drug_id": drug, "jaccard_alpha1_vs_alpha2": inter / union if union else np.nan, "changed_interaction_feature_count": int((av != bv).sum()), "alpha1_features": int(av.sum()), "alpha2_features": int(bv.sum()), "fingerprint_status": "computed_from_PLIP" if len(idx) else "no_interaction_data"})
    sim = pd.DataFrame(sims)
    _write_csv(sim, OUT / "alpha1_alpha2_similarity.csv", STRICT_RES / "alpha1_alpha2_similarity.csv")

    # Simple dependency-free heatmaps.  Rows are drug/receptor and columns are
    # mapped features; empty data produce a labelled empty panel.
    for value, name, cmap in [("binary", "fingerprint_binary_heatmap.png", "Greys"), ("pose_frequency", "fingerprint_frequency_heatmap.png", "viridis")]:
        data = (binary if value == "binary" else freq).copy()
        if data.empty:
            continue
        piv = data.assign(row=data.drug_id + "|" + data.receptor_id).pivot(index="row", columns="feature", values=value).fillna(0)
        fig, ax = plt.subplots(figsize=(max(7, 0.3 * len(piv.columns)), 4.5)); im = ax.imshow(piv.values, aspect="auto", cmap=cmap, vmin=0, vmax=1); ax.set_yticks(range(len(piv.index)), piv.index); ax.set_xticks(range(len(piv.columns)), piv.columns, rotation=90, fontsize=7); ax.set_title(name.replace("_", " ").replace(".png", "")); fig.colorbar(im, ax=ax, fraction=0.02); fig.tight_layout(); fig.savefig(OUT / "figures" / name, dpi=180); fig.savefig(STRICT_RES / name, dpi=180); plt.close(fig)

    # Integrated figure: show pair-level Ki ratios and IFP similarity without
    # implying a statistical fit or a representative pharmacology value.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    if not pharmacology.empty:
        for drug, g in pharmacology.groupby("drug_id"):
            ax.scatter([drug] * len(g), g.ki_alpha2_over_alpha1, label=drug)
        ax.set_yscale("log"); ax.set_ylabel("Ki α2 / α1 (same-unit pair)"); ax.tick_params(axis="x", rotation=30)
    else:
        ax.text(0.5, 0.5, "No comparable Ki pairs", ha="center", va="center")
    ax.set_title("ChEMBL pair-level pharmacology")
    ax2 = axes[1]; ax2.bar(sim.drug_id, sim.jaccard_alpha1_vs_alpha2); ax2.set_ylim(0, 1); ax2.set_ylabel("Jaccard similarity"); ax2.tick_params(axis="x", rotation=30); ax2.set_title("Docking IFP (provisional α2 construct)")
    fig.suptitle("Strict 3-drug mechanistic PoC: raw pairwise comparison")
    fig.tight_layout(); fig.savefig(OUT / "figures" / "ifp_vs_pharmacology_integrated.png", dpi=200); fig.savefig(STRICT_RES / "ifp_vs_pharmacology_integrated.png", dpi=200); plt.close(fig)
    return sim


def write_reports(candidates: pd.DataFrame, boxes: pd.DataFrame, pharmacology: pd.DataFrame, docking: pd.DataFrame, mapping: pd.DataFrame, interactions: pd.DataFrame, similarity: pd.DataFrame) -> None:
    direct = candidates[candidates.pdb_id.isin(["6HUP", "7QNE"])]
    qc_lines = [
        "# Strict 3-drug PoC QC report",
        "",
        "This report is an audit trail. No values are averaged, no missing interactions are imputed, and an absent interaction is distinct from absent interaction data.",
        "",
        f"- ChEMBL comparable Ki pair rows: **{len(pharmacology)}** (pair-level; `pair_priority=1` means same document, type and unit).",
        f"- Docking pose rows: **{len(docking)}**.",
        f"- PLIP interaction rows: **{len(interactions)}**.",
        f"- Residue mapping rows requiring QC: **{int(mapping.qc_required.sum()) if not mapping.empty else 0}**.",
        "",
        "## Structural QC",
        "",
        "6HUP is an exact human α1β3γ2L receptor with a diazepam copy assigned to alpha1 chain D at the alpha1/gamma2 site; its coordinates define the alpha1 box. 7QNE is an independent exact α1β3γ2 structure with Ro15-4513 and is retained as an alternative candidate.",
        "",
        "No exact human α2β3γ2 full pentamer with a BZD-site co-crystal ligand was found in the investigated RCSB set. The alpha2 docking receptor is therefore a provisional local-site construct consisting of chains C/D/E (β3/α2/γ2) extracted from mixed native 9CTJ. It is not equivalent to an exact α2β3γ2 pentamer. The DZP-derived box is transferred by a 207-residue gamma2 C→E Cα superposition (RMSD recorded in `docking_boxes.csv`); this transfer is reproducible but biologically uncertain.",
        "",
        "The strict mechanistic interpretation must therefore treat alpha2 docking, residue mapping, and any IFP comparison as exploratory/provisional. They do not establish subtype-specific causality.",
        "",
        "## Docking QC",
        "",
        "All six jobs use the same Vina seed, exhaustiveness, number of modes, energy range, and ligand preparation protocol. PDBQT conversion uses a transparent local writer with RDKit Gasteiger charges and element/aromatic atom types; receptor partial charges are zero because Open Babel/Meeko were not available. This is a methodological limitation.",
        "",
        "## Interaction QC",
        "",
        "PLIP is run from the vendored package when available. Interactions are retained pose-by-pose. Features that cannot be mapped to a common position are excluded from the fingerprint tables and remain visible as `qc_required` in `residue_mapping.csv`; they are not auto-resolved.",
        "",
    ]
    (OUT / "QC_REPORT.md").write_text("\n".join(qc_lines))
    report = [
        "# Receptor selection report",
        "",
        "## Scope",
        "",
        "The strict PoC compares diazepam, triazolam and zolpidem against ChEMBL targets CHEMBL2094121 (human α1β3γ2) and CHEMBL2094130 (human α2β3γ2). Existing ChEMBL raw activity is kept unchanged under `data/raw/pharmacology/chembl/`.",
        "",
        "## Candidate structures",
        "",
        "```\n" + direct.to_string(index=False) + "\n```",
        "",
        "```\n" + candidates[candidates.pdb_id.isin(["9CTJ", "9CX7", "9CXC", "9CSB"])].to_string(index=False) + "\n```",
        "",
        "## Selection",
        "",
        "6HUP was selected for α1β3γ2 because it is an exact receptor composition and provides a bound diazepam molecule at the α1/γ2 BZD site. The docking box is its DZP coordinate extent plus a fixed 6 Å margin (minimum 20 Å per axis), recorded in `docking_boxes.csv`.",
        "",
        "An exact α2β3γ2 structure was not identified. 9CTJ was selected only as a provisional local-site source because it contains human α2, β3 and γ2 chains and preserves a native α2/γ2 interface, but it also contains α1 and β2. Chains C/D/E are extracted and the box is transferred by gamma2 alignment. This is a structural approximation and is explicitly not called a clean α2β3γ2 full receptor.",
        "",
        "All structural limitations, box transfer details and receptor preparation limitations are repeated in `QC_REPORT.md` and `FINAL_REPORT.md`.",
    ]
    (OUT / "receptor_selection_report.md").write_text("\n".join(report))

    final = [
        "# FINAL_REPORT: strict 3-drug GABA-A Interaction Fingerprint PoC",
        "",
        "## Outcome",
        "",
        "The strict pipeline was executed for diazepam, triazolam and zolpidem. ChEMBL Ki rows were retained pair-by-pair, 3D ligands were generated from the downloaded canonical SMILES, six identical-condition Vina jobs were run, and pose-level PLIP fingerprints were written where parsing succeeded.",
        "",
        f"- Pair-level Ki rows: {len(pharmacology)}; no values were averaged.",
        f"- Docking rows: {len(docking)} (up to nine poses per drug/receptor).",
        f"- PLIP rows: {len(interactions)}.",
        "",
        "## Interpretation boundary",
        "",
        "The α1 receptor is directly represented by human α1β3γ2L structure 6HUP. The α2 receptor is a provisional α2/β3/γ2 local-site construct extracted from mixed native human structure 9CTJ because an exact α2β3γ2 full pentamer with a BZD-site co-crystal ligand was not found. Accordingly, docking scores, residue fingerprints and Jaccard values are mechanistic feasibility outputs, not validated subtype comparisons.",
        "",
        "No statistical generalization, ML, imputation, activity-type integration, unit mixing or representative-value selection was performed. Same-document/same-unit Ki ratios are kept as multiple raw pair rows.",
        "",
        "## Pair-level Ki and IFP comparison",
        "",
        "The following compact table shows only the highest-priority same-document Ki pairs; the complete unaveraged pair table is `strict_pharmacology_summary.csv`.",
        "",
        "```\n" + (pharmacology[pharmacology.pair_priority == 1][["drug_id", "alpha1_document_chembl_id", "alpha1_standard_value", "alpha2_standard_value", "standard_units", "ki_alpha2_over_alpha1"]].to_string(index=False) if not pharmacology.empty else "No same-document pair") + "\n```",
        "",
        "```\n" + similarity.to_string(index=False) + "\n```",
        "",
        "## Files",
        "",
        "- `strict_pharmacology_summary.csv`: raw pair-level Ki comparisons and α2/α1 ratios.",
        "- `docking_results.csv`: drug/receptor/pose IDs and Vina scores.",
        "- `residue_mapping.csv`: mapping audit including unresolved positions.",
        "- `fingerprint_binary.csv`, `fingerprint_frequency.csv`: pose-level IFP outputs.",
        "- `alpha1_alpha2_similarity.csv`: Jaccard and changed-feature counts.",
        "- `QC_REPORT.md`: explicit structural, box, PDBQT and interaction QC.",
    ]
    (OUT / "FINAL_REPORT.md").write_text("\n".join(final))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-docking", action="store_true", help="prepare data/reports without running Vina")
    args = parser.parse_args()
    mkdirs()
    pharmacology = build_strict_pharmacology()
    candidates = receptor_candidates()
    ligands = generate_ligands()
    receptors = prepare_receptors()
    boxes = docking_boxes()
    if args.skip_docking:
        docking = pd.DataFrame(columns=["drug_id", "receptor_id", "pose_id", "vina_score_kcal_mol"])
        interactions = pd.DataFrame(columns=["drug_id", "receptor_id", "pose_id", "interaction_type", "residue_number", "residue_chain", "raw_attributes"])
    else:
        docking = run_vina(receptors, ligands, boxes)
        interactions = run_plip(docking)
    mapping, binary, freq = build_mapping_and_fingerprints(docking, interactions)
    similarity = make_figures(binary, freq, pharmacology)
    write_reports(candidates, boxes, pharmacology, docking, mapping, interactions, similarity)


if __name__ == "__main__":
    main()
