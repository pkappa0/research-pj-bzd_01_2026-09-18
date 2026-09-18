"""ID-driven, lossless ChEMBL acquisition in immutable run directories.

Uses only the Python standard library. API JSON bytes are gzip archived before
parsing; all derived rows point back to a response and JSON pointer. No pooling,
unit conversion, name-based identification, or missing-value imputation.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urljoin, urlsplit
from urllib.request import Request, urlopen

from . import run_manager

ROOT = Path(__file__).resolve().parents[1]
PROVENANCE = ['source_request_id', 'source_file', 'source_json_pointer', 'source_sha256', 'retrieved_at_utc']


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')


def write_csv(path, rows, required=()):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys([*required, *(k for r in rows for k in r)]))
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False, sort_keys=True)
                             if isinstance(v, (dict, list)) else v for k, v in row.items()})


def read_registry(path):
    with path.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError('Compound registry is empty')
    seen_drug, seen_molecule = set(), set()
    for row in rows:
        drug, mid = row.get('drug_id', ''), row.get('molecule_chembl_id', '')
        if not drug or not re.fullmatch(r'CHEMBL[0-9]+', mid):
            raise ValueError('Every compound requires drug_id and explicit molecule_chembl_id; names are never resolved')
        if drug in seen_drug or mid in seen_molecule:
            raise ValueError('Duplicate drug_id or ChEMBL ID in registry; resolve explicitly')
        seen_drug.add(drug)
        seen_molecule.add(mid)
        key = row.get('expected_standard_inchikey', '')
        if key and not re.fullmatch(r'[A-Z]{14}-[A-Z]{10}-[A-Z]', key):
            raise ValueError(f'Invalid expected InChIKey for {mid}')
    return rows


class ArchiveClient:
    def __init__(self, run_dir, config, replay=None):
        self.run_dir, self.config = run_dir, config
        self.base = config['api_base'].rstrip('/') + '/'
        self.entries = []
        self.replay = replay
        self.replay_index = {}
        if replay:
            for line in (replay / 'raw/requests.jsonl').read_text().splitlines():
                entry = json.loads(line)
                if entry.get('success'):
                    self.replay_index.setdefault(entry['url'], entry)

    def record(self, entry):
        self.entries.append(entry)
        with (self.run_dir / 'raw/requests.jsonl').open('a') as handle:
            handle.write(json.dumps(entry, ensure_ascii=False) + '\n')

    def get(self, endpoint, params=None):
        url = urljoin(self.base, endpoint)
        if params:
            url += ('&' if '?' in url else '?') + urlencode(sorted(params.items()))
        if urlsplit(url).scheme != 'https' or urlsplit(url).netloc != urlsplit(self.base).netloc:
            raise ValueError(f'Unexpected API/pagination host: {url}')
        request_id = f'request_{len(self.entries) + 1:05d}'
        entry = {'request_id': request_id, 'url': url, 'method': 'GET', 'started_at_utc': now(),
                 'success': False, 'attempts': []}
        payload = None
        try:
            if self.replay:
                old = self.replay_index[url]
                archive = self.replay / old['source_file']
                archive_bytes = archive.read_bytes()
                if digest(archive_bytes) != old['archive_sha256']:
                    raise ValueError(f'Archive hash mismatch: {archive}')
                payload = gzip.decompress(archive_bytes)
                if digest(payload) != old['source_sha256']:
                    raise ValueError(f'Response hash mismatch: {archive}')
                entry.update({k: old[k] for k in ['retrieved_at_utc', 'http_status', 'response_headers']})
                entry['replayed_from'] = {'run_id': self.replay.name, 'request_id': old['request_id']}
            else:
                for attempt in range(self.config['attempts']):
                    try:
                        request = Request(url, headers={'Accept': 'application/json',
                                                        'User-Agent': 'bzd-poc-acquisition/1.0'})
                        with urlopen(request, timeout=self.config['timeout_seconds']) as response:
                            payload = response.read()
                            entry.update(http_status=response.status, retrieved_at_utc=now(),
                                         response_headers={k: response.headers.get(k) for k in
                                                           ['Content-Type', 'ETag', 'Last-Modified', 'Date']})
                        # Retry malformed/truncated JSON, but never treat it as an empty page.
                        json.loads(payload)
                        break
                    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError, OSError) as exc:
                        entry['attempts'].append({'attempt': attempt + 1, 'error': str(exc)})
                        if isinstance(exc, HTTPError) and exc.code not in (429, 500, 502, 503, 504):
                            raise
                        if attempt + 1 == self.config['attempts']:
                            raise
                        time.sleep(min(2 ** attempt, 8))
                time.sleep(self.config.get('request_delay_seconds', 0))
            relative = f'raw/responses/{request_id}.json.gz'
            path = self.run_dir / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            archive_bytes = gzip.compress(payload, mtime=0)
            path.write_bytes(archive_bytes)
            entry.update(source_file=relative, source_sha256=digest(payload), archive_sha256=digest(archive_bytes))
            data = json.loads(payload)
            if not isinstance(data, dict):
                raise ValueError('API response is not an object')
            entry['success'] = True
            self.record(entry)
            return data, entry
        except Exception as exc:
            entry['error'] = f'{type(exc).__name__}: {exc}'
            self.record(entry)
            raise

    @staticmethod
    def sourced(record, entry, pointer):
        return dict(record, source_request_id=entry['request_id'], source_file=entry['source_file'],
                    source_json_pointer=pointer, source_sha256=entry['source_sha256'],
                    retrieved_at_utc=entry['retrieved_at_utc'])

    def detail(self, resource, identifier):
        data, entry = self.get(f'{resource}/{identifier}.json')
        return self.sourced(data, entry, '')

    def all(self, resource, params):
        params = dict(params, limit=self.config['page_size'], offset=0)
        endpoint, first, seen, rows, expected = f'{resource}.json', True, set(), [], None
        key = {'activity': 'activities', 'target': 'targets', 'assay': 'assays',
               'document': 'documents'}[resource]
        while True:
            data, entry = self.get(endpoint, params if first else None)
            first = False
            if entry['url'] in seen:
                raise ValueError('Pagination cycle')
            seen.add(entry['url'])
            if not isinstance(data.get(key), list) or 'page_meta' not in data:
                raise ValueError(f'Missing {key}/page_meta in response')
            count = data['page_meta'].get('total_count')
            if not isinstance(count, int) or (expected is not None and count != expected):
                raise ValueError('Invalid or changing page total_count')
            expected = count
            rows.extend(self.sourced(r, entry, f'/{key}/{i}') for i, r in enumerate(data[key]))
            endpoint = data['page_meta'].get('next')
            if not endpoint:
                break
            if not data[key] or len(rows) >= expected:
                raise ValueError('Non-progressing or inconsistent pagination')
        if len(rows) != expected:
            raise ValueError(f'Incomplete {resource}: expected {expected}, received {len(rows)}')
        return rows


def subunit_mentions(text):
    """Literal text hints only. Never promoted to measured receptor composition."""
    text = (text or '').lower()
    for greek, latin in [('α', 'alpha'), ('β', 'beta'), ('γ', 'gamma'), ('δ', 'delta')]:
        text = text.replace(greek, latin)
    numbered = {''.join(parts) for parts in re.findall(r'(alpha|beta|gamma)[-\s]?(\d+)', text)}
    unnumbered = set(re.findall(r'\b(?:delta|epsilon|theta|pi)\b', text))
    return sorted(numbered | unnumbered)


def target_scope(target):
    name = (target.get('pref_name') or '').lower()
    if 'benzodiazepine' in name and 'peripheral' in name and 'central' in name:
        return 'benzodiazepine_mixed_central_peripheral_review'
    if 'associated protein' in name or 'translocator' in name or 'peripheral' in name:
        return 'excluded_non_gabaa'
    if any(x in name for x in ['gaba-a', 'gaba a receptor', 'gaba(a)', 'gamma-aminobutyric acid receptor']):
        if 'rho' in name:
            return 'gabaa_rho_separate_review'
        return 'gabaa_related_site_unverified'
    if 'benzodiazepine' in name:
        return 'benzodiazepine_target_review'
    return 'excluded_unresolved'


def check_identity(registry_row, molecule):
    if molecule.get('molecule_chembl_id') != registry_row['molecule_chembl_id']:
        return 'blocked_id_mismatch'
    actual = (molecule.get('molecule_structures') or {}).get('standard_inchi_key')
    expected = registry_row.get('expected_standard_inchikey')
    if expected and actual != expected:
        return 'blocked_inchikey_mismatch_or_missing'
    return 'id_and_inchikey_verified' if expected else 'explicit_id_only_structure_review'


def batch_records(client, resource, identifiers):
    ids, rows = sorted(identifiers), []
    size = client.config['batch_size']
    for start in range(0, len(ids), size):
        rows.extend(client.all(resource, {f'{resource}_chembl_id__in': ','.join(ids[start:start + size])}))
    return rows


def numeric_issue(value):
    if value in (None, ''):
        return 'missing_standard_value'
    try:
        if not Decimal(str(value)).is_finite():
            return 'nonfinite_standard_value'
    except InvalidOperation:
        return 'nonnumeric_standard_value'
    return None


def activity_qc(activity, assay, target, document, duplicate_count):
    flags = []
    issue = numeric_issue(activity.get('standard_value'))
    if issue:
        flags.append(issue)
    for field in ['standard_units', 'standard_type', 'standard_relation', 'assay_chembl_id', 'document_chembl_id']:
        if activity.get(field) in (None, ''):
            flags.append('missing_' + field)
    if activity.get('standard_relation') not in (None, '', '='):
        flags.append('censored_or_nonexact_relation')
    if activity.get('data_validity_comment'):
        flags.append('source_data_validity_comment')
    if activity.get('potential_duplicate'):
        flags.append('source_potential_duplicate')
    if duplicate_count > 1:
        flags.append('repeated_activity_id')
    if not assay:
        flags.append('missing_assay_detail')
    elif assay.get('target_chembl_id') != activity.get('target_chembl_id'):
        flags.append('assay_activity_target_id_conflict')
    if not document:
        flags.append('missing_document_detail')
    if target.get('target_type') != 'PROTEIN COMPLEX':
        flags.append('not_specific_protein_complex')
    target_subunits = set(subunit_mentions(target.get('pref_name')))
    assay_subunits = set(subunit_mentions((assay or {}).get('description')))
    if not target_subunits:
        flags.append('target_subunits_unspecified')
    if assay_subunits and target_subunits and assay_subunits != target_subunits:
        flags.append('assay_target_subunit_text_discrepancy_review')
    for family in ['alpha', 'beta', 'gamma']:
        target_family = {x for x in target_subunits if x.startswith(family)}
        assay_family = {x for x in assay_subunits if x.startswith(family)}
        if target_family and assay_family and target_family != assay_family:
            flags.append('assay_target_subunit_conflict_review')
            break
    if not (assay or {}).get('assay_parameters'):
        flags.append('no_structured_assay_parameters')
    # Conditions in free text stay verbatim. A matched target does not establish BZD-site binding.
    return flags


def export_tables(run_dir, registry, molecules, targets, activities, assays, documents, version):
    processed, tables = run_dir / 'processed', run_dir / 'tables'
    reg = {r['molecule_chembl_id']: r for r in registry}
    molrows, compound_qc = [], []
    for m in molecules:
        mid = m['molecule_chembl_id']
        r = reg[mid]
        s, h = m.get('molecule_structures') or {}, m.get('molecule_hierarchy') or {}
        status = check_identity(r, m)
        molrows.append(dict(m, drug_id=r['drug_id'], identity_status=status,
                            canonical_smiles=s.get('canonical_smiles'), standard_inchikey=s.get('standard_inchi_key'),
                            parent_chembl_id=h.get('parent_chembl_id'), active_chembl_id=h.get('active_chembl_id'),
                            chembl_version=version))
        compound_qc.append({'drug_id': r['drug_id'], 'molecule_chembl_id': mid, 'identity_status': status,
                            'missing_smiles': not bool(s.get('canonical_smiles')),
                            'missing_inchikey': not bool(s.get('standard_inchi_key')),
                            'multicomponent_smiles': '.' in (s.get('canonical_smiles') or ''),
                            'different_parent_id': bool(h.get('parent_chembl_id') and h['parent_chembl_id'] != mid),
                            'chirality_as_reported': m.get('chirality'),
                            'activity_rows': sum(a.get('molecule_chembl_id') == mid for a in activities),
                            'ifp_available': r.get('ifp_available')})
    amap = {a['assay_chembl_id']: a for a in assays}
    tmap = {t['target_chembl_id']: t for t in targets}
    dmap = {d['document_chembl_id']: d for d in documents}
    counts = Counter(a.get('activity_id') for a in activities)
    joined, row_qc = [], []
    for a in activities:
        assay, target, doc = amap.get(a.get('assay_chembl_id')), tmap[a['target_chembl_id']], dmap.get(a.get('document_chembl_id'))
        row = dict(a, drug_id=reg[a['molecule_chembl_id']]['drug_id'], chembl_version=version,
                   target_scope=target_scope(target), target_type=target.get('target_type'),
                   target_record_pref_name=target.get('pref_name'),
                   target_component_descriptions=[c.get('component_description') for c in target.get('target_components', [])],
                   target_name_subunit_mentions=subunit_mentions(target.get('pref_name')),
                   assay_description_subunit_mentions=subunit_mentions((assay or {}).get('description')),
                   bzd_site_status='not_adjudicated')
        for prefix, record in [('assay_detail', assay), ('target_detail', target), ('document_detail', doc)]:
            if record:
                # Full target records (including xrefs) live once in targets.jsonl/CSV.
                selected = record if prefix != 'target_detail' else {k: record.get(k) for k in
                    ['target_chembl_id', 'pref_name', 'target_type', 'organism', 'species_group_flag', *PROVENANCE]}
                row.update({prefix + '_' + k: v for k, v in selected.items()})
        flags = activity_qc(a, assay, target, doc, counts[a.get('activity_id')])
        row['qc_flags'] = flags
        joined.append(row)
        row_qc.append({'activity_id': a.get('activity_id'), 'molecule_chembl_id': a.get('molecule_chembl_id'),
                       'source_request_id': a['source_request_id'], 'source_json_pointer': a['source_json_pointer'],
                       'qc_flags': flags})
    for name, rows, required in [
        ('compounds', molrows, ['drug_id', 'molecule_chembl_id', 'canonical_smiles', 'standard_inchikey']),
        ('targets', targets, ['target_chembl_id']), ('assays', assays, ['assay_chembl_id']),
        ('documents', documents, ['document_chembl_id']),
        ('activities', joined, ['activity_id', 'molecule_chembl_id', 'assay_chembl_id', 'target_chembl_id'])]:
        write_csv(processed / f'{name}.csv', rows, required)
        with (processed / f'{name}.jsonl').open('w') as handle:
            for row in rows:
                handle.write(json.dumps(row, ensure_ascii=False) + '\n')
    write_csv(tables / 'compound_qc.csv', compound_qc)
    write_csv(tables / 'activity_qc.csv', row_qc, ['activity_id', 'qc_flags'])
    review_fields = ['activity_id', 'drug_id', 'molecule_chembl_id', 'target_chembl_id',
                     'target_record_pref_name', 'assay_chembl_id', 'assay_detail_description',
                     'target_name_subunit_mentions', 'assay_description_subunit_mentions',
                     'document_chembl_id', 'qc_flags', *PROVENANCE]
    write_csv(tables / 'subunit_review.csv', [{k: r.get(k) for k in review_fields} for r in joined
              if 'assay_target_subunit_text_discrepancy_review' in r['qc_flags']], review_fields)
    breakdown = Counter((a['molecule_chembl_id'], a.get('target_chembl_id'), a.get('standard_type'),
                         a.get('standard_units'), a.get('standard_relation'), a.get('target_organism')) for a in activities)
    write_csv(tables / 'coverage.csv', [dict(zip(
        ['molecule_chembl_id', 'target_chembl_id', 'standard_type', 'standard_units', 'standard_relation', 'target_organism', 'rows'], (*key, n)))
        for key, n in sorted(breakdown.items(), key=lambda x: str(x[0]))], ['molecule_chembl_id', 'rows'])
    flags = Counter(f for row in joined for f in row['qc_flags'])
    write_csv(tables / 'qc_summary.csv', [{'flag': k, 'rows': n} for k, n in sorted(flags.items())], ['flag', 'rows'])
    return {'compounds': len(molecules), 'targets_in_scope': sum(not target_scope(t).startswith('excluded') for t in targets),
            'activities': len(activities), 'assays': len(assays), 'documents': len(documents),
            'identity_blocked': sum(r['identity_status'].startswith('blocked') for r in compound_qc),
            'missing_structures': sum(r['missing_smiles'] or r['missing_inchikey'] for r in compound_qc),
            'flags': dict(flags), 'compound_qc': compound_qc}


def verify_archive(run_dir):
    entries = [json.loads(line) for line in (run_dir / 'raw/requests.jsonl').read_text().splitlines()]
    for entry in entries:
        if not entry['success']:
            continue
        archive = (run_dir / entry['source_file']).read_bytes()
        if digest(archive) != entry['archive_sha256'] or digest(gzip.decompress(archive)) != entry['source_sha256']:
            raise ValueError(f'Corrupt raw response: {entry["source_file"]}')
    manifest_path = run_dir / 'run_manifest.json'
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        for relative, expected in manifest.get('artifact_sha256', {}).items():
            if digest((run_dir / relative).read_bytes()) != expected:
                raise ValueError(f'Corrupt artifact: {relative}')
    return sum(e['success'] for e in entries)


def run(config_path, replay=None):
    if replay:
        config_path = replay / 'config/chembl_acquisition.json'
    config = json.loads(config_path.read_text())
    registry_path = replay / 'config/compound_registry.csv' if replay else ROOT / config['registry']
    registry = read_registry(registry_path)
    for k in ['page_size', 'batch_size', 'timeout_seconds', 'attempts']:
        if not isinstance(config[k], (int, float)) or config[k] <= 0:
            raise ValueError(f'Invalid {k}')
    code_paths = [Path(__file__), ROOT / 'src/run_manager.py']
    run_dir, manifest = run_manager.create_run(
        short_phase_name='chembl_acquisition', analysis_phase='ID-driven compound and pharmacology acquisition',
        analysis_purpose='Expand the data-ready compound panel; preserve assay context and provenance for QC.',
        parent_run=config.get('parent_run'), drugs=[r['drug_id'] for r in registry], receptors={},
        docking_parameters={}, input_files=[], notes=['No name matching, pooling, imputation, or model training.'],
        source_run=replay.name if replay else None, code_paths=code_paths)
    for p, dest in [(config_path, 'config/chembl_acquisition.json'), (registry_path, 'config/compound_registry.csv'),
                    *[(p, f'config/code_snapshot/{p.name}') for p in code_paths]]:
        (run_dir / dest).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, run_dir / dest)
        manifest['input_files'].append({'source': str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p),
                                        'copy': dest, 'sha256': digest(p.read_bytes())})
    manifest['mode'] = 'offline_replay' if replay else 'live'
    manifest['schema_version'] = 1
    client = ArchiveClient(run_dir, config, replay)
    print(f'Run: {run_dir.name}', flush=True)
    try:
        status, _ = client.get('status.json')
        version = status.get('chembl_db_version')
        manifest['chembl_status'] = status
        if status.get('status') != 'UP' or not version:
            raise ValueError('ChEMBL status/version unavailable')
        targets = {}
        for term in config['target_search_terms']:
            for target in client.all('target', {'pref_name__icontains': term}):
                targets.setdefault(target['target_chembl_id'], target)
        for tid in config.get('target_ids', []):
            targets[tid] = client.detail('target', tid)
        target_ids = sorted(tid for tid, t in targets.items() if not target_scope(t).startswith('excluded'))
        if not target_ids:
            raise ValueError('No target candidates: refuse misleading zero-activity result')
        write_csv(run_dir / 'tables/target_scope.csv',
                  [dict(target_chembl_id=t['target_chembl_id'], pref_name=t.get('pref_name'), scope=target_scope(t))
                   for t in targets.values()])
        print(f'Targets in search scope: {len(target_ids)}', flush=True)
        molecules, activities = [], []
        for r in registry:
            mid = r['molecule_chembl_id']
            molecule = client.detail('molecule', mid)
            # Preserve requested ID separately if the service unexpectedly returns another ID.
            if molecule.get('molecule_chembl_id') != mid:
                raise ValueError(f'Molecule ID mismatch for {mid}')
            molecules.append(molecule)
            identity = check_identity(r, molecule)
            before = len(activities)
            if not identity.startswith('blocked'):
                for start in range(0, len(target_ids), config['batch_size']):
                    batch = target_ids[start:start + config['batch_size']]
                    records = client.all('activity', {'molecule_chembl_id': mid, 'target_chembl_id__in': ','.join(batch)})
                    if any(a.get('molecule_chembl_id') != mid or a.get('target_chembl_id') not in batch for a in records):
                        raise ValueError('API filter returned unexpected molecule or target')
                    activities.extend(records)
            print(f'{r["drug_id"]} {mid}: {identity}; {len(activities)-before} activities', flush=True)
        assays = batch_records(client, 'assay', {a['assay_chembl_id'] for a in activities if a.get('assay_chembl_id')})
        print(f'Assay details: {len(assays)}', flush=True)
        documents = batch_records(client, 'document', {a['document_chembl_id'] for a in activities if a.get('document_chembl_id')})
        final_status, _ = client.get('status.json')
        if final_status.get('chembl_db_version') != version:
            raise ValueError('ChEMBL release changed during run')
        summary = export_tables(run_dir, registry, molecules, list(targets.values()), activities, assays, documents, version)
        summary['verified_raw_responses'] = verify_archive(run_dir)
        hard_flags = ['missing_assay_detail', 'missing_document_detail', 'repeated_activity_id', 'assay_activity_target_id_conflict']
        status_name = 'PARTIAL_QC_FAILURE' if summary['identity_blocked'] or summary['missing_structures'] or any(
            summary['flags'].get(k) for k in hard_flags) else 'COMPLETE_WITH_REVIEW_FLAGS'
        write_json(run_dir / 'tables/acquisition_summary.json', summary)
        manifest['summary'] = summary
        report = ['# ChEMBL acquisition', '', f'Run: {run_dir.name}', f'ChEMBL: {version}',
                  f'Status: {status_name}', '',
                  f'Compounds: {len(molecules)}; activity rows: {len(activities)}; assays: {len(assays)}; documents: {len(documents)}.', '',
                  '| Compound | ChEMBL ID | Activities | Identity |', '|---|---|---:|---|']
        report += [f'| {r["drug_id"]} | {r["molecule_chembl_id"]} | {r["activity_rows"]} | {r["identity_status"]} |' for r in summary['compound_qc']]
        report += ['', 'Data acquisition completed only within the documented target-name query scope; this is not an exhaustive BZD-site dataset.',
                   'GABA-A/rho/benzodiazepine candidates are kept with scope labels; BZD-site binding is not automatically assigned.',
                   'Subunit text mentions are review hints; target components do not prove assay stoichiometry.',
                   'Missing source values stay null (JSONL) or blank (CSV); raw and standard endpoints, units, relations and duplicate flags remain separate.',
                   'Assay descriptions, structured parameters, target components and document identifiers are retained without pooling.',
                   'No activity is excluded for censoring, missing values, validity comments, or potential duplication.',
                   'Zero retrieved activities means no rows within this query scope, not biological inactivity.',
                   'No new docking/IFP/phenotype labels or ML models were produced. Review experimental contexts before joining.',
                   '', 'See tables/qc_summary.csv, tables/activity_qc.csv, tables/coverage.csv and raw/requests.jsonl.']
        (run_dir / 'FINAL_REPORT.md').write_text('\n'.join(report) + '\n')
        (run_dir / 'QC_REPORT.md').write_text('# Acquisition QC\n\n' + json.dumps(summary, indent=2) + '\n')
        manifest['artifact_sha256'] = {str(p.relative_to(run_dir)): digest(p.read_bytes())
                                       for p in sorted(run_dir.rglob('*')) if p.is_file() and p.name != 'run_manifest.json'}
        run_manager.finalize_run(run_dir, manifest, qc_status=status_name)
        return run_dir, status_name
    except Exception as exc:
        manifest['failure'] = f'{type(exc).__name__}: {exc}'
        (run_dir / 'QC_REPORT.md').write_text('# FAILED acquisition\n\n' + manifest['failure'] + '\nPartial raw responses preserved. No complete dataset claim.\n')
        (run_dir / 'FINAL_REPORT.md').write_text('# FAILED acquisition\n\nSee QC_REPORT.md and raw/requests.jsonl.\n')
        run_manager.finalize_run(run_dir, manifest, qc_status='FAILED')
        raise


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, default=ROOT / 'config/chembl_acquisition.json')
    ap.add_argument('--replay-run', type=Path, help='Create a new run without network using archived responses')
    ap.add_argument('--verify-run', type=Path, help='Verify raw and artifact hashes without modifying a run')
    args = ap.parse_args()
    if args.verify_run:
        print(f'Verified {verify_archive(args.verify_run)} raw responses')
        return
    run_dir, status = run(args.config.resolve(), args.replay_run.resolve() if args.replay_run else None)
    print(json.dumps({'run_id': run_dir.name, 'qc_status': status}))
    if status == 'PARTIAL_QC_FAILURE':
        sys.exit(2)


if __name__ == '__main__':
    main()
