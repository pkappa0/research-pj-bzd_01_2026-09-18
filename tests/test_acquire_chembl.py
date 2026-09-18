"""Synthetic API fixtures only; never used as research observations."""
import csv
import gzip
import json
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pytest

from src import acquire_chembl as ac


KEY = 'AAAAAAAAAAAAAA-BBBBBBBBBB-C'
CONFIG = {'api_base': 'https://example.org/api/', 'registry': 'registry.csv',
          'target_search_terms': ['GABA-A'], 'target_ids': [], 'parent_run': None,
          'page_size': 2, 'batch_size': 20, 'timeout_seconds': 1, 'attempts': 1,
          'request_delay_seconds': 0}


def make_registry(path, rows=None):
    ac.write_csv(path, rows or [{'drug_id': 'synthetic', 'molecule_chembl_id': 'CHEMBL1',
                               'expected_standard_inchikey': KEY}], ['drug_id', 'molecule_chembl_id'])


def test_registry_requires_id_not_name(tmp_path):
    path = tmp_path / 'registry.csv'
    make_registry(path, [{'drug_id': 'name', 'input_drug_name': 'Diazepam'}])
    with pytest.raises(ValueError, match='explicit'):
        ac.read_registry(path)


def test_registry_rejects_duplicate_id(tmp_path):
    path = tmp_path / 'registry.csv'
    make_registry(path, [{'drug_id': name, 'molecule_chembl_id': 'CHEMBL1'} for name in ['a', 'b']])
    with pytest.raises(ValueError, match='Duplicate'):
        ac.read_registry(path)


def test_identity_does_not_collapse_parent_or_stereoisomer():
    r = {'molecule_chembl_id': 'CHEMBL1', 'expected_standard_inchikey': KEY}
    m = {'molecule_chembl_id': 'CHEMBL2', 'molecule_hierarchy': {'parent_chembl_id': 'CHEMBL1'}}
    assert ac.check_identity(r, m) == 'blocked_id_mismatch'
    m['molecule_chembl_id'] = 'CHEMBL1'
    m['molecule_structures'] = {'standard_inchi_key': KEY.replace('BBBBBBBBBB', 'CCCCCCCCCC')}
    assert ac.check_identity(r, m).startswith('blocked_inchikey')
    m['molecule_structures']['standard_inchi_key'] = KEY
    assert ac.check_identity(r, m) == 'id_and_inchikey_verified'


def test_subunits_do_not_parse_gamma_acid_or_pi_inside_words():
    assert ac.subunit_mentions('Gamma-aminobutyric acid receptor; alpha-1/beta-3/gamma-2') == ['alpha1', 'beta3', 'gamma2']
    assert ac.subunit_mentions('Binding to alpha1beta2gamma2; hippocampus; zolpidem') == ['alpha1', 'beta2', 'gamma2']
    assert ac.subunit_mentions('α2β3γ2') == ['alpha2', 'beta3', 'gamma2']
    assert ac.subunit_mentions(None) == []


def test_scope_separates_peripheral_and_rho():
    assert ac.target_scope({'pref_name': 'Peripheral benzodiazepine receptor'}) == 'excluded_non_gabaa'
    assert ac.target_scope({'pref_name': 'Benzodiazepine receptors; peripheral & central'}) == 'benzodiazepine_mixed_central_peripheral_review'
    assert ac.target_scope({'pref_name': 'Gamma-aminobutyric acid receptor subunit rho-1'}) == 'gabaa_rho_separate_review'


def test_qc_preserves_censoring_missing_and_conflicts():
    a = {'standard_relation': '>', 'standard_value': None, 'standard_units': None,
         'target_chembl_id': 'CHEMBL2', 'potential_duplicate': 1}
    assay = {'target_chembl_id': 'CHEMBL2', 'description': 'alpha1beta2gamma2'}
    t = {'pref_name': 'GABA-A alpha1/beta3/gamma2', 'target_type': 'PROTEIN COMPLEX'}
    flags = ac.activity_qc(a, assay, t, {}, 2)
    assert 'missing_standard_value' in flags
    assert 'missing_standard_units' in flags
    assert 'censored_or_nonexact_relation' in flags
    assert 'assay_target_subunit_text_discrepancy_review' in flags
    assert 'repeated_activity_id' in flags
    assert a['standard_value'] is None and a['standard_relation'] == '>'


class Response:
    status = 200
    headers = {'Content-Type': 'application/json'}

    def __init__(self, data):
        self.data = json.dumps(data).encode()

    def read(self):
        return self.data

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


@pytest.fixture
def client(tmp_path):
    (tmp_path / 'raw').mkdir()
    return ac.ArchiveClient(tmp_path, CONFIG)


def page(rows, total, next_page=None, key='activities'):
    return {key: rows, 'page_meta': {'total_count': total, 'next': next_page}}


def test_pagination_raw_bytes_and_pointer(client, monkeypatch):
    responses = iter([page([{'activity_id': 1}], 2, '/api/activity.json?offset=1'),
                      page([{'activity_id': 2}], 2)])
    monkeypatch.setattr(ac, 'urlopen', lambda *a, **kw: Response(next(responses)))
    rows = client.all('activity', {})
    assert [r['activity_id'] for r in rows] == [1, 2]
    assert rows[1]['source_json_pointer'] == '/activities/0'
    raw = gzip.decompress((client.run_dir / rows[1]['source_file']).read_bytes())
    assert json.loads(raw)['activities'][0]['activity_id'] == 2
    assert ac.verify_archive(client.run_dir) == 2


@pytest.mark.parametrize('data', [page([], 1), page([{}], 2), {'error': 'not an empty dataset'}])
def test_incomplete_or_invalid_pages_fail(client, monkeypatch, data):
    monkeypatch.setattr(ac, 'urlopen', lambda *a, **kw: Response(data))
    with pytest.raises(ValueError):
        client.all('activity', {})


def test_cross_host_pagination_rejected(client):
    with pytest.raises(ValueError, match='host'):
        client.get('https://other.example/api/activity.json')


def test_replay_and_corruption(client, monkeypatch, tmp_path):
    monkeypatch.setattr(ac, 'urlopen', lambda *a, **kw: Response(page([], 0)))
    client.all('activity', {})
    other = tmp_path / 'other'
    (other / 'raw').mkdir(parents=True)
    monkeypatch.setattr(ac, 'urlopen', lambda *a, **kw: pytest.fail('offline replay used network'))
    replay = ac.ArchiveClient(other, CONFIG, client.run_dir)
    assert replay.all('activity', {}) == []
    archive = client.run_dir / client.entries[0]['source_file']
    archive.write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='Corrupt'):
        ac.verify_archive(client.run_dir)


@pytest.fixture
def integration_env(tmp_path, monkeypatch):
    monkeypatch.setattr(ac, 'ROOT', tmp_path)
    monkeypatch.setattr(ac.run_manager, 'ROOT', tmp_path)
    monkeypatch.setattr(ac.run_manager, 'RUNS', tmp_path / 'runs')
    (tmp_path / 'src').mkdir()
    (tmp_path / 'src/run_manager.py').write_text('# synthetic snapshot')
    make_registry(tmp_path / 'registry.csv')
    cfg = tmp_path / 'config.json'
    ac.write_json(cfg, CONFIG)
    target = {'target_chembl_id': 'CHEMBL2', 'pref_name': 'GABA-A alpha1/beta3/gamma2',
              'target_type': 'PROTEIN COMPLEX', 'target_components': []}
    def api(request, **kwargs):
        path = urlsplit(request.full_url).path
        if path.endswith('status.json'):
            return Response({'status': 'UP', 'chembl_db_version': 'SYNTHETIC_TEST'})
        if path.endswith('target.json'):
            return Response(page([target], 1, key='targets'))
        if path.endswith('molecule/CHEMBL1.json'):
            return Response({'molecule_chembl_id': 'CHEMBL1', 'molecule_structures':
                             {'standard_inchi_key': KEY, 'canonical_smiles': 'C'}})
        if path.endswith('activity.json'):
            return Response(page([{'activity_id': 1, 'molecule_chembl_id': 'CHEMBL1',
                                   'target_chembl_id': 'CHEMBL2', 'assay_chembl_id': 'CHEMBL3',
                                   'document_chembl_id': 'CHEMBL4', 'standard_relation': '<',
                                   'standard_value': None}], 1))
        if path.endswith('assay.json'):
            return Response(page([{'assay_chembl_id': 'CHEMBL3', 'target_chembl_id': 'CHEMBL2',
                                   'description': 'Synthetic fixture only', 'assay_parameters': []}], 1, key='assays'))
        if path.endswith('document.json'):
            return Response(page([{'document_chembl_id': 'CHEMBL4', 'doi': None}], 1, key='documents'))
        pytest.fail(request.full_url)
    monkeypatch.setattr(ac, 'urlopen', api)
    return cfg


def test_end_to_end_replay_does_not_mutate_source(integration_env, monkeypatch):
    first, status = ac.run(integration_env)
    assert status == 'COMPLETE_WITH_REVIEW_FLAGS'
    before = {str(p.relative_to(first)): p.read_bytes() for p in first.rglob('*') if p.is_file()}
    monkeypatch.setattr(ac, 'urlopen', lambda *a, **kw: pytest.fail('replay used network'))
    second, status = ac.run(integration_env, first)
    assert second != first
    assert {str(p.relative_to(first)): p.read_bytes() for p in first.rglob('*') if p.is_file()} == before
    assert (first / 'processed/activities.jsonl').read_bytes() == (second / 'processed/activities.jsonl').read_bytes()
    activity = json.loads((second / 'processed/activities.jsonl').read_text())
    assert activity['standard_value'] is None
    assert activity['standard_relation'] == '<'
    assert activity['document_detail_doi'] is None


def test_failed_run_is_explicit(integration_env, monkeypatch):
    monkeypatch.setattr(ac, 'urlopen', lambda *a, **kw: Response({'error': 'API unavailable'}))
    with pytest.raises(ValueError, match='status/version'):
        ac.run(integration_env)
    manifests = list(ac.run_manager.RUNS.glob('*/run_manifest.json'))
    assert len(manifests) == 1
    assert json.loads(manifests[0].read_text())['qc_status'] == 'FAILED'


def test_transient_http_error_retries_and_records_attempt(client, monkeypatch):
    client.config = dict(CONFIG, attempts=2)
    calls = []
    def fetch(*args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise ac.HTTPError('https://example.org/api/activity.json', 503, 'Unavailable', {}, None)
        return Response(page([], 0))
    monkeypatch.setattr(ac, 'urlopen', fetch)
    monkeypatch.setattr(ac.time, 'sleep', lambda _: None)
    assert client.all('activity', {}) == []
    assert len(calls) == 2
    assert len(client.entries[0]['attempts']) == 1


def test_processed_artifact_tampering_detected(integration_env):
    run_dir, _ = ac.run(integration_env)
    with (run_dir / 'processed/activities.csv').open('a') as handle:
        handle.write('tampered\n')
    with pytest.raises(ValueError, match='Corrupt artifact'):
        ac.verify_archive(run_dir)
