"""GtoPdb pharmacology supplementation and ChEMBL integration.

The GtoPdb API requires a user API key.  This module never stores that key and
does not write it to URLs, JSON, reports, or logs.  When the key is unavailable
the run produces explicit QC outputs and can still consume the official
download CSV at ``data/raw/gtopdb/gtopdb_interactions.csv``.

The API and CSV paths are normalized to the same schema.  Raw API JSON and raw
CSV rows are retained; no missing activity is converted to zero and no affinity
types, units, species, or beta subtypes are merged.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urljoin
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
DRUG_LIST = ROOT / "data" / "raw" / "drug_list.csv"
CHEMBL_MOLECULES = ROOT / "data" / "raw" / "pharmacology" / "chembl" / "molecules.csv"
CHEMBL_ACTIVITY = ROOT / "data" / "processed" / "pharmacology.csv"
DEFAULT_RAW = ROOT / "data" / "raw" / "gtopdb"
DEFAULT_RESULTS = ROOT / "results" / "tables"
DEFAULT_OUTPUTS = ROOT / "outputs"

NORM_COLUMNS = [
    "drug_id", "input_drug_name", "gtopdb_ligand_id", "gtopdb_ligand_name",
    "gtopdb_target_id", "target_name", "target_type", "species",
    "interaction_type", "action", "affinity_parameter", "affinity_value",
    "affinity_low", "affinity_high", "affinity_units", "pKi", "pKd",
    "pIC50", "pEC50", "alpha_subunit", "beta_subunit", "gamma_subunit",
    "other_subunits", "receptor_composition", "alpha_class", "is_gabaa_related",
    "reference_id", "PMID", "DOI", "reference_citation", "assay_id",
    "source", "retrieval_date", "gtopdb_version", "raw_record_json",
]

MASTER_COLUMNS = [
    "evidence_id", "drug_id", "input_drug_name", "molecule_chembl_id",
    "gtopdb_ligand_id", "receptor_composition", "target_id", "target_name",
    "target_type", "species", "metric", "value", "value_low", "value_high",
    "unit", "assay", "reference_id", "PMID", "DOI", "reference",
    "source_database", "source_record_id", "alpha_subtype", "beta_subtype",
    "gamma_subtype", "other_subunits", "duplicate_group_id",
    "duplicate_candidate",
]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def norm_text(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip()).casefold()


def first_value(obj: Mapping[str, Any], aliases: Sequence[str], default: Any = "") -> Any:
    """Get a field while tolerating API/official-download naming variants."""
    wanted = {re.sub(r"[^a-z0-9]", "", a.casefold()) for a in aliases}
    for key, value in obj.items():
        if re.sub(r"[^a-z0-9]", "", str(key).casefold()) in wanted:
            return value
    return default


def stringify(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def numeric(value: Any) -> Any:
    if value in (None, "", "NA", "N/A", "null"):
        return np.nan
    try:
        return float(str(value).replace(",", "").strip())
    except (TypeError, ValueError):
        return np.nan


def flatten_records(payload: Any) -> List[Dict[str, Any]]:
    """Extract likely interaction/candidate objects without assuming one schema."""
    if isinstance(payload, list):
        return [x for x in payload if isinstance(x, dict)]
    if not isinstance(payload, dict):
        return []
    preferred = ("interactions", "interaction", "results", "data", "items", "objects", "ligands")
    for key in preferred:
        value = payload.get(key)
        if isinstance(value, list):
            return [x for x in value if isinstance(x, dict)]
        if isinstance(value, dict):
            nested = flatten_records(value)
            if nested:
                return nested
    interaction_keys = {"target", "targetid", "targetname", "affinity", "action", "interactiontype", "ligandid"}
    if any(re.sub(r"[^a-z0-9]", "", str(k).casefold()) in interaction_keys for k in payload):
        return [payload]
    found: List[Dict[str, Any]] = []
    for value in payload.values():
        if isinstance(value, dict):
            found.extend(flatten_records(value))
        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    found.extend(flatten_records(item))
    return found


def parse_subunits(target_name: Any, subunit_value: Any = "") -> Dict[str, str]:
    text = f"{target_name or ''} {subunit_value or ''}"
    text = text.replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    found: Dict[str, List[str]] = {"alpha": [], "beta": [], "gamma": [], "other": []}
    pattern = re.compile(r"(alpha|beta|gamma|delta|epsilon|theta|pi)\s*[-_/ ]?\s*(\d+)", re.I)
    for match in pattern.finditer(text):
        family = match.group(1).lower()
        value = f"{family}{match.group(2)}"
        if value not in found[family if family in found else "other"]:
            found[family if family in found else "other"].append(value)
    alpha = "/".join(found["alpha"])
    beta = "/".join(found["beta"])
    gamma = "/".join(found["gamma"])
    other = "/".join(found["other"])
    pieces = [x for x in [alpha, beta, gamma, other] if x]
    return {
        "alpha_subunit": alpha,
        "beta_subunit": beta,
        "gamma_subunit": gamma,
        "other_subunits": other,
        "receptor_composition": "/".join(pieces) if pieces else "generic_gabaa",
        "alpha_class": "alpha1-containing" if "alpha1" in alpha else ("alpha2-containing" if "alpha2" in alpha else "alpha-unknown"),
    }


def is_gabaa(target_name: Any) -> bool:
    text = norm_text(target_name).replace("α", "alpha")
    if re.search(r"gaba\s*[- ]?b|gaba\s*[- ]?c", text):
        return False
    return bool(re.search(r"gaba\s*[- ]?a|gabaa|gamma-aminobutyric acid receptor", text))


def species_is_human(value: Any) -> bool:
    text = norm_text(value)
    return text in {"human", "homo sapiens", "h. sapiens"} or "homo sapiens" in text


def _reference_fields(record: Mapping[str, Any]) -> Dict[str, str]:
    ref = first_value(record, ["reference", "citation", "sourceReference", "publication"], "")
    if isinstance(ref, dict):
        pmid = first_value(ref, ["PMID", "pmid", "pubmedId", "pubmed_id"], "")
        doi = first_value(ref, ["DOI", "doi"], "")
        rid = first_value(ref, ["referenceId", "reference_id", "id"], "")
        citation = first_value(ref, ["citation", "title", "referenceCitation"], "")
    else:
        pmid = first_value(record, ["PMID", "pmid", "pubmedId", "pubmed_id"], "")
        doi = first_value(record, ["DOI", "doi"], "")
        rid = first_value(record, ["referenceId", "reference_id", "reference_id"], "")
        citation = ref
    text = " ".join([stringify(ref), stringify(citation)])
    if not pmid:
        m = re.search(r"\bPMID\s*[:#]?\s*(\d{5,9})\b", text, re.I)
        pmid = m.group(1) if m else ""
    if not doi:
        m = re.search(r"10\.\d{4,9}/[-._;()/:a-z0-9]+", text, re.I)
        doi = m.group(0) if m else ""
    return {"reference_id": stringify(rid), "PMID": stringify(pmid), "DOI": stringify(doi), "reference_citation": stringify(citation or ref)}


def normalize_interaction(record: Mapping[str, Any], drug_id: str, input_name: str, ligand_id: str = "", ligand_name: str = "", retrieval_date: str = "", version: str = "") -> Dict[str, Any]:
    target = first_value(record, ["target", "targetName", "target_name", "targetPrefName"], "")
    target_id = first_value(record, ["targetId", "target_id", "gtopdbTargetId", "target_chembl_id"], "")
    target_type = first_value(record, ["targetType", "target_type", "type"], "")
    if isinstance(target, dict):
        target_id = target_id or first_value(target, ["id", "targetId", "target_id"], "")
        target_type = target_type or first_value(target, ["type", "targetType", "target_type"], "")
        target = first_value(target, ["name", "targetName", "prefName", "label"], "")
    subunit_value = first_value(record, ["subunits", "subunit", "receptorComposition", "receptor_composition", "composition"], "")
    subs = parse_subunits(target, subunit_value)
    # Official downloads may provide separate subtype columns instead of
    # embedding the composition in target_name.
    direct_alpha = stringify(first_value(record, ["alpha_subunit", "alphaSubtype", "alpha_subtype"], ""))
    direct_beta = stringify(first_value(record, ["beta_subunit", "betaSubtype", "beta_subtype"], ""))
    direct_gamma = stringify(first_value(record, ["gamma_subunit", "gammaSubtype", "gamma_subtype"], ""))
    direct_other = stringify(first_value(record, ["other_subunits", "otherSubunits"], ""))
    if direct_alpha or direct_beta or direct_gamma or direct_other:
        subs["alpha_subunit"] = direct_alpha
        subs["beta_subunit"] = direct_beta
        subs["gamma_subunit"] = direct_gamma
        subs["other_subunits"] = direct_other
        pieces = [x for x in [direct_alpha, direct_beta, direct_gamma, direct_other] if x]
        subs["receptor_composition"] = "/".join(pieces) if pieces else "generic_gabaa"
        subs["alpha_class"] = "alpha1-containing" if "alpha1" in direct_alpha else ("alpha2-containing" if "alpha2" in direct_alpha else "alpha-unknown")
    reference = _reference_fields(record)
    interaction_type = first_value(record, ["interactionType", "interaction_type", "type"], "")
    action = first_value(record, ["action", "effect", "activity", "modeOfAction"], "")
    parameter = first_value(record, ["affinityParameter", "affinity_parameter", "parameter", "measure", "standardType"], "")
    value = first_value(record, ["affinityValue", "affinity_value", "value", "standardValue", "standard_value"], "")
    low = first_value(record, ["affinityLow", "affinity_low", "lower", "low"], "")
    high = first_value(record, ["affinityHigh", "affinity_high", "upper", "high"], "")
    units = first_value(record, ["affinityUnits", "affinity_units", "units", "unit", "standardUnits"], "")
    species = first_value(record, ["species", "organism", "targetSpecies"], "")
    row = {
        "drug_id": drug_id, "input_drug_name": input_name, "gtopdb_ligand_id": stringify(ligand_id), "gtopdb_ligand_name": stringify(ligand_name),
        "gtopdb_target_id": stringify(target_id), "target_name": stringify(target), "target_type": stringify(target_type), "species": stringify(species),
        "interaction_type": stringify(interaction_type), "action": stringify(action), "affinity_parameter": stringify(parameter), "affinity_value": numeric(value),
        "affinity_low": numeric(low), "affinity_high": numeric(high), "affinity_units": stringify(units),
        "pKi": numeric(first_value(record, ["pKi", "pki", "pKiValue", "pki_value"], "")), "pKd": numeric(first_value(record, ["pKd", "pkd", "pKdValue", "pkd_value"], "")),
        "pIC50": numeric(first_value(record, ["pIC50", "pic50", "pIC50Value", "pic50_value"], "")), "pEC50": numeric(first_value(record, ["pEC50", "pec50", "pEC50Value", "pec50_value"], "")),
        **subs, "is_gabaa_related": is_gabaa(target), **reference, "assay_id": stringify(first_value(record, ["assayId", "assay_id", "assay"], "")),
        "source": "GtoPdb", "retrieval_date": retrieval_date, "gtopdb_version": version, "raw_record_json": json.dumps(record, ensure_ascii=False, sort_keys=True),
    }
    return row


class GtoPdbClient:
    def __init__(self, base_url: str, api_key: str, raw_dir: Path, timeout: float = 30, delay: float = 0.1):
        self.base_url = base_url.rstrip("/") + "/"
        self.api_key = api_key
        self.raw_dir = raw_dir
        self.timeout = timeout
        self.delay = delay

    def get(self, path: str, params: Optional[Mapping[str, Any]] = None, raw_name: Optional[str] = None) -> Any:
        url = urljoin(self.base_url, path.lstrip("/"))
        if params:
            url += ("&" if "?" in url else "?") + urlencode({k: v for k, v in params.items() if v not in (None, "")})
        req = Request(url, headers={"Accept": "application/json", "GTP-API-Key": self.api_key, "User-Agent": "gaba-a-fingerprint-poc/0.2"})
        try:
            with urlopen(req, timeout=self.timeout) as response:
                payload = json.load(response)
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
            code = getattr(exc, "code", "")
            raise RuntimeError(f"GtoPdb request failed ({code}): {path}") from exc
        if raw_name:
            api_dir = self.raw_dir / "api"
            api_dir.mkdir(parents=True, exist_ok=True)
            (api_dir / f"{raw_name}.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2))
        if self.delay:
            time.sleep(self.delay)
        return payload


def read_drugs() -> pd.DataFrame:
    if not DRUG_LIST.exists():
        return pd.DataFrame(columns=["drug_id", "input_drug_name"])
    df = pd.read_csv(DRUG_LIST, dtype="string")
    return df[["drug_id", "input_drug_name"]].assign(drug_id=lambda x: x.drug_id.str.strip().str.lower(), input_drug_name=lambda x: x.input_drug_name.str.strip())


def load_chembl_molecules() -> pd.DataFrame:
    if not CHEMBL_MOLECULES.exists():
        return pd.DataFrame(columns=["drug_id", "molecule_chembl_id", "canonical_smiles", "standard_inchikey"])
    return pd.read_csv(CHEMBL_MOLECULES, dtype="string")


def ligand_candidate_rows(payload: Any) -> List[Dict[str, Any]]:
    rows = []
    for item in flatten_records(payload):
        ligand_id = first_value(item, ["ligandId", "ligand_id", "id", "ligand"], "")
        name = first_value(item, ["name", "ligandName", "ligand_name", "preferredName", "label"], "")
        if isinstance(ligand_id, dict):
            ligand_id = first_value(ligand_id, ["id", "ligandId"], "")
        if ligand_id or name:
            rows.append({"gtopdb_ligand_id": stringify(ligand_id), "gtopdb_ligand_name": stringify(name), "canonical_smiles": stringify(first_value(item, ["canonicalSmiles", "canonical_smiles", "smiles"], "")), "standard_inchikey": stringify(first_value(item, ["inchiKey", "inchikey", "standardInchiKey", "standard_inchikey"], "")), "raw": item})
    return rows


def choose_ligand(candidates: List[Dict[str, Any]], expected: Mapping[str, Any]) -> Tuple[Optional[Dict[str, Any]], str, str]:
    exp_key = norm_text(expected.get("standard_inchikey", "")).upper()
    exp_smiles = norm_text(expected.get("canonical_smiles", ""))
    exp_name = norm_text(expected.get("input_drug_name", ""))
    chemical = []
    for candidate in candidates:
        key_match = bool(exp_key and norm_text(candidate.get("standard_inchikey", "")).upper() == exp_key)
        smiles_match = bool(exp_smiles and norm_text(candidate.get("canonical_smiles", "")) == exp_smiles)
        if key_match or smiles_match:
            chemical.append(candidate)
    unique = {x.get("gtopdb_ligand_id"): x for x in chemical if x.get("gtopdb_ligand_id")}
    if len(unique) == 1:
        return next(iter(unique.values())), "exact_chemical_identifier_match", ""
    if len(unique) > 1:
        return None, "ambiguous_chemical_match", "multiple GtoPdb candidates match ChEMBL chemical identifiers"
    name_only = [x for x in candidates if norm_text(x.get("gtopdb_ligand_name")) == exp_name]
    if name_only:
        return None, "name_match_unconfirmed", "name matched but canonical SMILES/InChIKey did not confirm identity"
    return None, "unresolved", "no exact chemical identifier match"


def fallback_rows(path: Path, drugs: pd.DataFrame, molecules: pd.DataFrame, retrieval_date: str, version: str) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], str]:
    if not path.exists():
        return [], [], "missing"
    df = pd.read_csv(path, dtype="string")
    required = ["drug_id", "input_drug_name", "target_name"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        return [], [], "invalid_missing_columns:" + ",".join(missing)
    raw, normalized = [], []
    for rec in df.fillna("").to_dict("records"):
        raw.append(rec)
        did = rec.get("drug_id", "")
        row = normalize_interaction(rec, did, rec.get("input_drug_name", ""), rec.get("gtopdb_ligand_id", ""), rec.get("gtopdb_ligand_name", ""), retrieval_date, version)
        normalized.append(row)
    return raw, normalized, "loaded"


def add_chembl_master(master_rows: List[Dict[str, Any]], molecules: pd.DataFrame) -> None:
    if not CHEMBL_ACTIVITY.exists():
        return
    chem = pd.read_csv(CHEMBL_ACTIVITY, dtype="string")
    mol_by_drug = molecules.set_index("drug_id").to_dict("index") if not molecules.empty else {}
    for idx, rec in chem.fillna("").iterrows():
        composition = rec.get("receptor_subtype", "") or parse_subunits(rec.get("target_name", "")).get("receptor_composition", "generic_gabaa")
        subs = parse_subunits(composition)
        ref = rec.get("reference", "") or rec.get("document_chembl_id", "")
        master_rows.append({
            "evidence_id": f"chembl:{idx}", "drug_id": rec.get("drug_id", ""), "input_drug_name": rec.get("input_drug_name", ""),
            "molecule_chembl_id": rec.get("molecule_chembl_id", ""), "gtopdb_ligand_id": "", "receptor_composition": composition or "generic_gabaa",
            "target_id": rec.get("target_chembl_id", ""), "target_name": rec.get("target_name", ""), "target_type": rec.get("target_type", ""),
            "species": rec.get("organism", ""), "metric": rec.get("standard_type", ""), "value": numeric(rec.get("standard_value", "")),
            "value_low": np.nan, "value_high": np.nan, "unit": rec.get("standard_units", ""), "assay": rec.get("assay_chembl_id", ""),
            "reference_id": rec.get("document_chembl_id", ""), "PMID": "", "DOI": "", "reference": ref,
            "source_database": "ChEMBL", "source_record_id": rec.get("activity_id", idx), "alpha_subtype": subs.get("alpha_subunit", ""),
            "beta_subtype": subs.get("beta_subunit", ""), "gamma_subtype": subs.get("gamma_subunit", ""), "other_subunits": subs.get("other_subunits", ""),
            "duplicate_group_id": "", "duplicate_candidate": False,
        })


def add_gtopdb_master(master_rows: List[Dict[str, Any]], normalized: pd.DataFrame) -> None:
    for idx, rec in normalized.fillna("").iterrows():
        if not bool(rec.get("is_gabaa_related", False)):
            continue
        metric = rec.get("affinity_parameter", "")
        value = rec.get("affinity_value", np.nan)
        if pd.isna(value):
            for p in ["pKi", "pKd", "pIC50", "pEC50"]:
                if not pd.isna(rec.get(p, np.nan)):
                    metric, value = p, rec.get(p)
                    break
        master_rows.append({
            "evidence_id": f"gtopdb:{idx}", "drug_id": rec.get("drug_id", ""), "input_drug_name": rec.get("input_drug_name", ""),
            "molecule_chembl_id": "", "gtopdb_ligand_id": rec.get("gtopdb_ligand_id", ""), "receptor_composition": rec.get("receptor_composition", "generic_gabaa"),
            "target_id": rec.get("gtopdb_target_id", ""), "target_name": rec.get("target_name", ""), "target_type": rec.get("target_type", ""),
            "species": rec.get("species", ""), "metric": metric, "value": numeric(value), "value_low": rec.get("affinity_low", np.nan), "value_high": rec.get("affinity_high", np.nan),
            "unit": rec.get("affinity_units", ""), "assay": rec.get("assay_id", ""), "reference_id": rec.get("reference_id", ""),
            "PMID": rec.get("PMID", ""), "DOI": rec.get("DOI", ""), "reference": rec.get("reference_citation", ""), "source_database": "GtoPdb",
            "source_record_id": idx, "alpha_subtype": rec.get("alpha_subunit", ""), "beta_subtype": rec.get("beta_subunit", ""), "gamma_subtype": rec.get("gamma_subunit", ""), "other_subunits": rec.get("other_subunits", ""),
            "duplicate_group_id": "", "duplicate_candidate": False,
        })


def _reference_tokens(row: Mapping[str, Any]) -> set:
    tokens = set()
    for field in ["PMID", "DOI"]:
        value = norm_text(row.get(field, ""))
        if value:
            tokens.add(f"{field}:{value}")
    text = stringify(row.get("reference", ""))
    pmid = re.search(r"\bPMID\s*[:#]?\s*(\d{5,9})\b", text, re.I)
    doi = re.search(r"10\.\d{4,9}/[-._;()/:a-z0-9]+", text, re.I)
    if pmid:
        tokens.add(f"PMID:{pmid.group(1).casefold()}")
    if doi:
        tokens.add(f"DOI:{doi.group(0).casefold()}")
    return tokens


def mark_duplicates(master: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if master.empty:
        return master.assign(duplicate_group_id=pd.Series(dtype="string"), duplicate_candidate=pd.Series(dtype="bool")), pd.DataFrame(columns=list(master.columns) + ["duplicate_reason"])
    master = master.copy()
    qc_rows = []
    group_no = 0
    chem = master[master.source_database == "ChEMBL"]
    gtp = master[master.source_database == "GtoPdb"]
    for cidx, c in chem.iterrows():
        for gidx, g in gtp.iterrows():
            c_tokens, g_tokens = _reference_tokens(c), _reference_tokens(g)
            same_ref = bool(c_tokens & g_tokens)
            same_content = all([
                norm_text(c.drug_id) == norm_text(g.drug_id), norm_text(c.receptor_composition) == norm_text(g.receptor_composition),
                species_is_human(c.species) == species_is_human(g.species), norm_text(c.metric) == norm_text(g.metric),
                numeric(c.value) == numeric(g.value), norm_text(c.unit) == norm_text(g.unit),
            ])
            if same_ref and same_content:
                group_no += 1
                gid = f"DUP_{group_no:04d}"
                master.loc[cidx, "duplicate_group_id"] = gid
                master.loc[gidx, "duplicate_group_id"] = gid
                master.loc[cidx, "duplicate_candidate"] = True
                master.loc[gidx, "duplicate_candidate"] = True
                qc_rows.append({"duplicate_group_id": gid, "chembl_evidence_id": c.evidence_id, "gtopdb_evidence_id": g.evidence_id, "duplicate_reason": "same PMID/DOI and drug/composition/species/metric/value/unit", "PMID": g.PMID, "DOI": g.DOI})
            elif same_content and (norm_text(c.reference) and norm_text(c.reference) == norm_text(g.reference)):
                group_no += 1
                gid = f"DUP_{group_no:04d}"
                master.loc[cidx, "duplicate_group_id"] = gid
                master.loc[gidx, "duplicate_group_id"] = gid
                master.loc[cidx, "duplicate_candidate"] = True
                master.loc[gidx, "duplicate_candidate"] = True
                qc_rows.append({"duplicate_group_id": gid, "chembl_evidence_id": c.evidence_id, "gtopdb_evidence_id": g.evidence_id, "duplicate_reason": "same reference text and drug/composition/species/metric/value/unit", "PMID": g.PMID, "DOI": g.DOI})
    return master, pd.DataFrame(qc_rows, columns=["duplicate_group_id", "chembl_evidence_id", "gtopdb_evidence_id", "duplicate_reason", "PMID", "DOI"])


def comparable_pairs(master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if master.empty:
        return pd.DataFrame(columns=["drug_id", "alpha1_composition", "alpha2_composition", "metric", "unit", "same_beta_gamma", "same_reference", "tier"])
    data = master[(~master.duplicate_candidate.astype(bool)) & master.value.notna()].copy()
    data["human"] = data.species.map(species_is_human)
    for drug in sorted(data.drug_id.dropna().unique()):
        sub = data[data.drug_id == drug]
        a1 = sub[sub.alpha_subtype.astype(str).str.contains("alpha1", na=False)]
        a2 = sub[sub.alpha_subtype.astype(str).str.contains("alpha2", na=False)]
        for _, x in a1.iterrows():
            for _, y in a2.iterrows():
                same_metric = norm_text(x.metric) == norm_text(y.metric) and norm_text(x.unit) == norm_text(y.unit)
                same_bg = norm_text(x.beta_subtype) == norm_text(y.beta_subtype) and norm_text(x.gamma_subtype) == norm_text(y.gamma_subtype) and bool(x.beta_subtype) and bool(x.gamma_subtype) and bool(y.beta_subtype) and bool(y.gamma_subtype)
                same_ref = bool(_reference_tokens(x) & _reference_tokens(y)) or (norm_text(x.reference) and norm_text(x.reference) == norm_text(y.reference)) or (norm_text(x.assay) and norm_text(x.assay) == norm_text(y.assay))
                if bool(x.human and y.human) and same_bg and same_metric and same_ref:
                    tier = 1
                elif bool(x.human and y.human) and same_bg and same_metric:
                    tier = 2
                elif bool(x.human and y.human) and same_metric:
                    tier = 3
                else:
                    tier = 4
                rows.append({"drug_id": drug, "alpha1_composition": x.receptor_composition, "alpha2_composition": y.receptor_composition, "alpha1_source": x.source_database, "alpha2_source": y.source_database, "metric": x.metric, "unit": x.unit, "same_beta_gamma": same_bg, "same_reference_or_assay": bool(same_ref), "alpha1_species": x.species, "alpha2_species": y.species, "alpha1_value": x.value, "alpha2_value": y.value, "tier": tier, "alpha1_evidence_id": x.evidence_id, "alpha2_evidence_id": y.evidence_id})
    return pd.DataFrame(rows)


def build_presence(drugs: pd.DataFrame, normalized: pd.DataFrame, master: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, drug in drugs.iterrows():
        did = drug.drug_id
        g = normalized[(normalized.drug_id == did) & normalized.is_gabaa_related.astype(bool)] if not normalized.empty else normalized
        quant = g[g.affinity_value.notna()] if not g.empty else g
        a1 = quant[quant.alpha_subunit.astype(str).str.contains("alpha1", na=False)] if not quant.empty else quant
        a2 = quant[quant.alpha_subunit.astype(str).str.contains("alpha2", na=False)] if not quant.empty else quant
        pp = pairs[pairs.drug_id == did] if not pairs.empty else pairs
        rows.append({"drug_id": did, "input_drug_name": drug.input_drug_name, "gtopdb_alpha1_quantitative": bool(not a1.empty), "gtopdb_alpha2_quantitative": bool(not a2.empty), "gtopdb_both": bool(not a1.empty and not a2.empty), "gtopdb_alpha1_targets": int(a1.gtopdb_target_id.nunique()) if not a1.empty else 0, "gtopdb_alpha2_targets": int(a2.gtopdb_target_id.nunique()) if not a2.empty else 0, "gtopdb_compositions": "|".join(sorted(set(quant.receptor_composition))) if not quant.empty else "", "same_beta_gamma_pair_exists": bool((pp.same_beta_gamma == True).any()) if not pp.empty else False, "same_metric_pair_exists": bool(not pp.empty), "same_reference_pair_exists": bool((pp.same_reference_or_assay == True).any()) if not pp.empty else False, "best_tier": int(pp.tier.min()) if not pp.empty else (5 if not g.empty else np.nan)})
    return pd.DataFrame(rows)


def write_reports(fetch_meta: Dict[str, Any], presence: pd.DataFrame, tier_summary: pd.DataFrame, duplicate_qc: pd.DataFrame, fallback_status: str, normalized: pd.DataFrame, master: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    key_status = fetch_meta.get("api_key_status", "unknown")
    lines = [
        "# GTOPDB_FETCH_REPORT", "", "This is an independent GtoPdb supplementation phase. Existing ChEMBL and strict 3-drug outputs were not modified.", "",
        f"- API key status: **{key_status}**", f"- Fallback CSV status: **{fallback_status}**", f"- GtoPdb normalized GABA-A rows: **{len(normalized)}**", f"- GtoPdb version: **{fetch_meta.get('gtopdb_version') or 'not exposed by response'}**", "",
    ]
    if key_status != "available":
        lines += ["The GtoPdb REST API was not called because `GTP_API_KEY` is absent or authentication failed. Set the key only in the process environment:", "", "```bash", "export GTP_API_KEY='your-key'", "python fetch_gtopdb.py --config config/poc.yaml", "```", "", "Do not place the key in code, YAML, CSV, reports, or Git history. Alternatively place the official download at `data/raw/gtopdb/gtopdb_interactions.csv` and rerun.", ""]
    (out_dir / "GTOPDB_FETCH_REPORT.md").write_text("\n".join(lines))

    combined_both = int(tier_summary.combined_alpha1_alpha2.astype(bool).sum()) if not tier_summary.empty else 0
    t1 = int((tier_summary.best_tier == 1).sum()) if not tier_summary.empty else 0
    t2 = int((tier_summary.best_tier <= 2).sum()) if not tier_summary.empty else 0
    t3 = int((tier_summary.best_tier <= 3).sum()) if not tier_summary.empty else 0
    strict_three = {"diazepam", "triazolam", "zolpidem"}
    supplemented = tier_summary[(tier_summary.gtopdb_both == True) & (~tier_summary.drug_id.isin(strict_three))].drug_id.tolist() if not tier_summary.empty else []
    gtopdb_missing = tier_summary[(tier_summary.gtopdb_alpha1_quantitative == False) | (tier_summary.gtopdb_alpha2_quantitative == False)].drug_id.tolist() if not tier_summary.empty else []
    report = [
        "# PHARMACOLOGY_INTEGRATION_REPORT", "", "## Scope and rules", "", "ChEMBL raw activity is retained as-is. GtoPdb rows are separate evidence records, and duplicate candidates are flagged rather than silently removed. Ki, Kd, IC50 and EC50 remain distinct; units, species and beta subtypes are not merged.", "",
        f"- Drugs with alpha1/alpha2 quantitative data in GtoPdb alone: **{int(presence.gtopdb_both.astype(bool).sum()) if not presence.empty else 0}**",
        f"- Drugs with alpha1/alpha2 quantitative data in the combined master: **{combined_both}**",
        f"- Tier 1: **{t1}**", f"- Tier 2 or better: **{t2}**", f"- Tier 3 or better: **{t3}**", f"- ChEMBL-strict missing but GtoPdb-complemented: **{', '.join(supplemented) if supplemented else 'none in current run'}**", f"- GtoPdb without subtype-specific quantitative data: **{', '.join(gtopdb_missing) if gtopdb_missing else 'none'}**", f"- Duplicate candidate rows: **{len(duplicate_qc)}**", "",
        "The tier is a data comparability classification, not a potency or efficacy ranking.", "",
        "## Per-drug QC", "", tier_summary.to_string(index=False) if not tier_summary.empty else "No records available.", "",
        "## Current limitation", "", "If the API key was unavailable and the fallback CSV was absent, these counts describe the existing ChEMBL master only; they do not claim that GtoPdb has no corresponding data. Rerun after providing the key or official CSV.", "",
    ]
    (out_dir / "PHARMACOLOGY_INTEGRATION_REPORT.md").write_text("\n".join(report))


def run(cfg: Mapping[str, Any]) -> Dict[str, Any]:
    gcfg = cfg.get("gtopdb", {})
    raw_dir = ROOT / gcfg.get("raw_output_dir", "data/raw/gtopdb")
    result_dir = ROOT / gcfg.get("results_dir", "results/tables")
    output_dir = ROOT / gcfg.get("outputs_dir", "outputs")
    raw_dir.mkdir(parents=True, exist_ok=True); (raw_dir / "api").mkdir(exist_ok=True); result_dir.mkdir(parents=True, exist_ok=True); output_dir.mkdir(parents=True, exist_ok=True)
    retrieval_date = now_utc()
    api_key = os.environ.get("GTP_API_KEY", "")
    fetch_meta: Dict[str, Any] = {"retrieval_date": retrieval_date, "api_key_status": "absent" if not api_key else "available", "gtopdb_version": ""}
    drugs = read_drugs(); molecules = load_chembl_molecules(); molmap = molecules.set_index("drug_id").to_dict("index") if not molecules.empty else {}
    ligand_qc = []
    raw_records: List[Dict[str, Any]] = []
    normalized_rows: List[Dict[str, Any]] = []
    fallback_path = ROOT / gcfg.get("fallback_csv", "data/raw/gtopdb/gtopdb_interactions.csv")
    if api_key:
        # API mode takes precedence when a key is present; an offline CSV is
        # not mixed into the API evidence set.
        fallback_raw, fallback_norm, fallback_status = [], [], "not_used_api_key"
    else:
        fallback_raw, fallback_norm, fallback_status = fallback_rows(fallback_path, drugs, molecules, retrieval_date, "")
        raw_records.extend(fallback_raw); normalized_rows.extend(fallback_norm)
    if api_key:
        client = GtoPdbClient(gcfg.get("api_base", "https://www.guidetopharmacology.org/services/"), api_key, raw_dir, float(gcfg.get("request_timeout_seconds", 30)), float(gcfg.get("request_delay_seconds", 0.1)))
        for _, drug in drugs.iterrows():
            did, name = drug.drug_id, drug.input_drug_name
            expected = molmap.get(did, {})
            try:
                search = client.get("ligands", {"name": name}, f"ligand_search_{did}")
                if isinstance(search, dict):
                    fetch_meta["gtopdb_version"] = stringify(first_value(search, ["gtopdbVersion", "gtopdb_version", "version", "release"], fetch_meta.get("gtopdb_version", "")))
                candidates = ligand_candidate_rows(search)
                # Search responses often expose only an ID and display name.
                # Fetch the candidate detail before matching, while retaining
                # every detail response in the raw API directory.
                for candidate in candidates:
                    if candidate.get("gtopdb_ligand_id") and not (candidate.get("canonical_smiles") or candidate.get("standard_inchikey")):
                        try:
                            detail_payload = client.get(f"ligands/{candidate['gtopdb_ligand_id']}", raw_name=f"ligand_detail_{did}_{candidate['gtopdb_ligand_id']}")
                            detail_items = flatten_records(detail_payload)
                            detail = detail_items[0] if detail_items else (detail_payload if isinstance(detail_payload, dict) else {})
                            candidate["canonical_smiles"] = stringify(first_value(detail, ["canonicalSmiles", "canonical_smiles", "smiles"], candidate.get("canonical_smiles", "")))
                            candidate["standard_inchikey"] = stringify(first_value(detail, ["inchiKey", "inchikey", "standardInchiKey", "standard_inchikey"], candidate.get("standard_inchikey", "")))
                            candidate["gtopdb_ligand_name"] = candidate.get("gtopdb_ligand_name") or stringify(first_value(detail, ["name", "ligandName", "preferredName", "label"], ""))
                        except RuntimeError:
                            candidate["detail_lookup_status"] = "failed"
                selected, status, note = choose_ligand(candidates, {**expected, "input_drug_name": name})
                ligand_qc.append({"drug_id": did, "input_drug_name": name, "status": status, "gtopdb_ligand_id": selected.get("gtopdb_ligand_id", "") if selected else "", "gtopdb_ligand_name": selected.get("gtopdb_ligand_name", "") if selected else "", "candidate_count": len(candidates), "chembl_molecule_id": expected.get("molecule_chembl_id", ""), "chembl_inchikey": expected.get("standard_inchikey", ""), "chembl_canonical_smiles": expected.get("canonical_smiles", ""), "note": note})
                if not selected:
                    continue
                lid = selected["gtopdb_ligand_id"]
                payload = client.get(f"ligands/{lid}/interactions", {"species": "Human", "targetType": "LGIC"}, f"interactions_{did}_{lid}")
                if isinstance(payload, dict):
                    fetch_meta["gtopdb_version"] = stringify(first_value(payload, ["gtopdbVersion", "gtopdb_version", "version", "release"], fetch_meta.get("gtopdb_version", "")))
                items = flatten_records(payload)
                for item in items:
                    raw_records.append({"drug_id": did, "input_drug_name": name, "gtopdb_ligand_id": lid, "gtopdb_ligand_name": selected.get("gtopdb_ligand_name", ""), "raw_record_json": json.dumps(item, ensure_ascii=False, sort_keys=True)})
                    normalized_rows.append(normalize_interaction(item, did, name, lid, selected.get("gtopdb_ligand_name", ""), retrieval_date, fetch_meta.get("gtopdb_version", "")))
            except RuntimeError as exc:
                text = str(exc)
                if "401" in text or "403" in text:
                    fetch_meta["api_key_status"] = "rejected_401_or_403"
                ligand_qc.append({"drug_id": did, "input_drug_name": name, "status": "api_error", "gtopdb_ligand_id": "", "gtopdb_ligand_name": "", "candidate_count": 0, "chembl_molecule_id": expected.get("molecule_chembl_id", ""), "chembl_inchikey": expected.get("standard_inchikey", ""), "chembl_canonical_smiles": expected.get("canonical_smiles", ""), "note": text})
    elif fallback_status == "loaded":
        fetch_meta["api_key_status"] = "absent_fallback_csv_used"
        for _, drug in drugs.iterrows():
            expected = molmap.get(drug.drug_id, {})
            ligand_qc.append({"drug_id": drug.drug_id, "input_drug_name": drug.input_drug_name, "status": "fallback_csv", "gtopdb_ligand_id": "", "gtopdb_ligand_name": "", "candidate_count": 0, "chembl_molecule_id": expected.get("molecule_chembl_id", ""), "chembl_inchikey": expected.get("standard_inchikey", ""), "chembl_canonical_smiles": expected.get("canonical_smiles", ""), "note": "Official fallback CSV used; ligand identity is carried by CSV rows"})
    elif not api_key:
        for _, drug in drugs.iterrows():
            expected = molmap.get(drug.drug_id, {})
            ligand_qc.append({"drug_id": drug.drug_id, "input_drug_name": drug.input_drug_name, "status": "api_key_missing", "gtopdb_ligand_id": "", "gtopdb_ligand_name": "", "candidate_count": 0, "chembl_molecule_id": expected.get("molecule_chembl_id", ""), "chembl_inchikey": expected.get("standard_inchikey", ""), "chembl_canonical_smiles": expected.get("canonical_smiles", ""), "note": "GTP_API_KEY is absent; no ligand was auto-confirmed"})

    raw_df = pd.DataFrame(raw_records)
    if raw_df.empty:
        raw_df = pd.DataFrame(columns=["drug_id", "input_drug_name", "gtopdb_ligand_id", "gtopdb_ligand_name", "raw_record_json"])
    norm_df = pd.DataFrame(normalized_rows, columns=NORM_COLUMNS)
    if not norm_df.empty:
        norm_df = norm_df.reindex(columns=NORM_COLUMNS)
    _write_csv(raw_df, raw_dir / "gtopdb_interactions_raw.csv", ROOT / "data" / "processed" / "gtopdb_interactions_raw.csv")
    _write_csv(norm_df, raw_dir / "gtopdb_interactions_normalized.csv", ROOT / "data" / "processed" / "gtopdb_interactions_normalized.csv")
    ligand_qc_df = pd.DataFrame(ligand_qc)
    _write_csv(ligand_qc_df, result_dir / "gtopdb_ligand_match_qc.csv")
    target_counts = norm_df[norm_df.is_gabaa_related.astype(bool)].groupby(["drug_id", "gtopdb_target_id", "target_name", "species"], dropna=False).size().reset_index(name="interaction_count") if not norm_df.empty else pd.DataFrame(columns=["drug_id", "gtopdb_target_id", "target_name", "species", "interaction_count"])
    comp = norm_df[norm_df.is_gabaa_related.astype(bool)].groupby(["drug_id", "receptor_composition", "alpha_class", "species"], dropna=False).size().reset_index(name="interaction_count") if not norm_df.empty else pd.DataFrame(columns=["drug_id", "receptor_composition", "alpha_class", "species", "interaction_count"])
    master_rows: List[Dict[str, Any]] = []; add_chembl_master(master_rows, molecules); add_gtopdb_master(master_rows, norm_df)
    master = pd.DataFrame(master_rows, columns=MASTER_COLUMNS)
    master, duplicate_qc = mark_duplicates(master)
    _write_csv(target_counts, result_dir / "gtopdb_drug_target_counts.csv"); _write_csv(comp, result_dir / "gtopdb_receptor_composition.csv"); _write_csv(master, ROOT / "data" / "processed" / "pharmacology_master_raw.csv"); _write_csv(duplicate_qc, result_dir / "pharmacology_duplicate_qc.csv")
    pairs_all = comparable_pairs(master)
    pairs_gtopdb = pairs_all[(pairs_all.alpha1_source == "GtoPdb") | (pairs_all.alpha2_source == "GtoPdb")] if not pairs_all.empty else pairs_all
    _write_csv(pairs_gtopdb, result_dir / "gtopdb_comparable_metrics.csv")
    _write_csv(pairs_all, result_dir / "pharmacology_comparable_metrics.csv")
    g_presence = build_presence(drugs, norm_df, master, pairs_gtopdb)
    _write_csv(g_presence, result_dir / "gtopdb_alpha1_alpha2_presence.csv")
    tier_rows = []
    for _, drug in drugs.iterrows():
        did = drug.drug_id
        p = pairs_all[pairs_all.drug_id == did] if not pairs_all.empty else pairs_all
        g = norm_df[(norm_df.drug_id == did) & norm_df.is_gabaa_related.astype(bool)] if not norm_df.empty else norm_df
        q = g[g.affinity_value.notna()] if not g.empty else g
        a1 = q[q.alpha_subunit.astype(str).str.contains("alpha1", na=False)] if not q.empty else q
        a2 = q[q.alpha_subunit.astype(str).str.contains("alpha2", na=False)] if not q.empty else q
        cm = master[(master.drug_id == did) & (~master.duplicate_candidate.astype(bool))] if not master.empty else master
        ca1 = cm[cm.alpha_subtype.astype(str).str.contains("alpha1", na=False)] if not cm.empty else cm
        ca2 = cm[cm.alpha_subtype.astype(str).str.contains("alpha2", na=False)] if not cm.empty else cm
        tier_rows.append({"drug_id": did, "input_drug_name": drug.input_drug_name, "gtopdb_alpha1_quantitative": bool(not a1.empty), "gtopdb_alpha2_quantitative": bool(not a2.empty), "gtopdb_both": bool(not a1.empty and not a2.empty), "combined_alpha1_quantitative": bool(not ca1.empty), "combined_alpha2_quantitative": bool(not ca2.empty), "combined_alpha1_alpha2": bool(not ca1.empty and not ca2.empty), "best_tier": int(p.tier.min()) if not p.empty else (5 if not g.empty else np.nan), "tier1_pairs": int((p.tier == 1).sum()) if not p.empty else 0, "tier2_pairs": int((p.tier == 2).sum()) if not p.empty else 0, "tier3_pairs": int((p.tier == 3).sum()) if not p.empty else 0, "chEMBL_strict_comparable": did in {"diazepam", "triazolam", "zolpidem"}, "gtopdb_subtype_specific_missing": bool(a1.empty or a2.empty)})
    tier_summary = pd.DataFrame(tier_rows)
    _write_csv(tier_summary, result_dir / "pharmacology_tier_summary.csv")
    fetch_meta["gtopdb_interaction_count"] = len(norm_df); fetch_meta["fallback_status"] = fallback_status; fetch_meta["retrieval_date"] = retrieval_date
    (raw_dir / "gtopdb_fetch_metadata.json").write_text(json.dumps(fetch_meta, ensure_ascii=False, indent=2))
    write_reports(fetch_meta, g_presence, tier_summary, duplicate_qc, fallback_status, norm_df, master, output_dir)
    return {"api_key_status": fetch_meta["api_key_status"], "fallback_status": fallback_status, "gtopdb_rows": len(norm_df), "master_rows": len(master), "duplicate_candidates": len(duplicate_qc)}


def _write_csv(df: pd.DataFrame, *paths: Path) -> None:
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch and integrate GtoPdb GABA-A pharmacology")
    parser.add_argument("--config", default="config/poc.yaml")
    args = parser.parse_args()
    cfg = yaml.safe_load((ROOT / args.config).read_text())
    print(json.dumps(run(cfg), ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
