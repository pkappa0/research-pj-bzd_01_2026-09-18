"""Evidence-conflict, numerical integrity, comparison and immutable-run regressions."""
from copy import deepcopy
import csv
import json
from pathlib import Path
import shutil

import pytest

from src import bzd_pharmacology_filter as bf


SOURCE = bf.ROOT / 'runs/20260918_1117_chembl_acquisition'


def example(alpha=1, beta=3, molecule='CHEMBL12', activity=1):
    composition = f'alpha{alpha}beta{beta}gamma2'
    desc = f'Displacement of [3H]flumazenil from human GABA-A {composition} expressed in HEK293 cells'
    target = dict(target_chembl_id=f'CHEMBL{alpha}{beta}', pref_name=f'GABA-A {composition}',
                  organism='Homo sapiens', target_type='PROTEIN COMPLEX', target_components=[
                      dict(component_description=t, relationship='PROTEIN SUBUNIT')
                      for t in [f'alpha{alpha}', f'beta{beta}', 'gamma2']])
    row = dict(activity_id=activity, drug_id='synthetic', molecule_chembl_id=molecule,
               target_chembl_id=target['target_chembl_id'], target_organism='Homo sapiens',
               assay_chembl_id=f'CHEMBL10{alpha}{beta}', document_chembl_id='CHEMBL999',
               assay_detail_target_chembl_id=target['target_chembl_id'], assay_detail_document_chembl_id='CHEMBL999',
               assay_detail_assay_chembl_id=f'CHEMBL10{alpha}{beta}', assay_type='B', assay_detail_assay_type='B',
               assay_description=desc, assay_detail_description=desc, assay_detail_assay_organism='Homo sapiens',
               standard_type='Ki', standard_units='nM', standard_value='12.00', standard_relation='=',
               type='Ki', units='nM', value='12.00', relation='=', potential_duplicate=0)
    return row, target


def description(row, text):
    row.update(assay_description=text, assay_detail_description=text)


def classified(**kwargs):
    return bf.classify_activity(*example(**kwargs))


@pytest.mark.parametrize('text,expected', [('α1β3γ2', 'alpha1_beta3_gamma2'),
    ('alpha-2-beta-2-gamma-2', 'alpha2_beta2_gamma2'), ('GABRA1 GABRB3 GABRG2', 'alpha1_beta3_gamma2'),
    ('Gamma-aminobutyric acid receptor, hippocampus, zolpidem', None)])
def test_explicit_subunit_tokens(text, expected):
    assert bf.extract_subunits(text)['composition'] == expected


@pytest.mark.parametrize('old,new,family', [('alpha1', 'alpha4', 'alpha'), ('alpha1', 'alpha5', 'alpha'), ('beta3', 'beta2', 'beta')])
def test_general_subunit_conflict_not_activity_specific(old, new, family):
    row, target = example(activity=777777)
    description(row, row['assay_description'].replace(old, new))
    result = bf.classify_activity(row, target)
    assert result['classification'] == 'REVIEW'
    assert f'target_assay_{family}_conflict' in result['review_flags']
    assert result['receptor_composition'] is None


def test_component_conflict_is_independent_evidence():
    row, target = example()
    target['target_components'][1]['component_description'] = 'beta2'
    result = bf.classify_activity(row, target)
    assert result['classification'] == 'REVIEW'
    assert 'target_name_component_beta_conflict' in result['review_flags']


def test_group_members_not_coassembled():
    row, target = example()
    target.update(pref_name='GABA-A receptor', target_type='PROTEIN COMPLEX GROUP', target_components=[
        dict(component_description='alpha1 alpha2 alpha4 beta2 beta3 gamma2', relationship='GROUP MEMBER')])
    result = bf.classify_activity(row, target)
    assert result['classification'] == 'SECONDARY'
    assert not result['target_components_are_coassembly_evidence']
    assert not result['review_flags']


def test_single_subunit_not_promoted_by_complete_assay():
    row, target = example()
    target.update(pref_name='GABA-A alpha1', target_type='SINGLE PROTEIN', target_components=[])
    assert bf.classify_activity(row, target)['classification'] == 'SECONDARY'


def test_mixed_alpha_not_assigned_to_two_values():
    row, target = example()
    description(row, row['assay_description'].replace('alpha1', 'alpha1/alpha2'))
    target['pref_name'] = 'GABA-A alpha1/alpha2 beta3 gamma2'
    target['target_components'] = []
    result = bf.classify_activity(row, target)
    assert result['subtype'] == 'alpha1_alpha2'
    assert result['classification'] == 'SECONDARY'
    assert not result['point_comparison_usable']


def test_species_missing_never_filled_from_prose_or_target():
    row, target = example()
    row['assay_detail_assay_organism'] = None
    result = bf.classify_activity(row, target)
    assert result['classification'] == 'PRIMARY'
    assert result['species'] == 'unknown' and result['target_species'] == 'human'
    assert not result['point_comparison_usable']


@pytest.mark.parametrize('organism', ['Rattus norvegicus', 'Xenopus laevis'])
def test_species_conflict_no_automatic_host_reinterpretation(organism):
    row, target = example()
    row['assay_detail_assay_organism'] = organism
    result = bf.classify_activity(row, target)
    assert result['classification'] == 'REVIEW' and result['species'] == 'unknown'
    assert result['assay_detail_assay_organism'] == organism


def test_host_cell_text_does_not_change_receptor_species():
    row, target = example()
    description(row, row['assay_description'].replace('HEK293 cells', 'mouse LTK cells'))
    result = bf.classify_activity(row, target)
    assert result['species'] == 'human' and result['classification'] == 'PRIMARY'


@pytest.mark.parametrize('typ,desc,assay,expected', [
    ('IC50', 'Displacement of [3H]flumazenil from human GABA-A alpha1beta3gamma2', 'B', 'binding'),
    ('IC50', 'Inhibition of GABA current in human GABA-A alpha1beta3gamma2', 'F', 'functional_potency'),
    ('EC50', 'Potentiation of GABA current in human GABA-A alpha1beta3gamma2', 'F', 'functional_potency'),
    ('Emax', 'Maximal potentiation of GABA current in human GABA-A alpha1beta3gamma2', 'F', 'functional_efficacy')])
def test_endpoint_uses_description(typ, desc, assay, expected):
    row, target = example()
    description(row, desc)
    row.update(standard_type=typ, assay_type=assay, assay_detail_assay_type=assay,
               standard_units='%' if typ == 'Emax' else 'nM')
    result = bf.classify_activity(row, target)
    assert result['endpoint_class'] == expected and result['classification'] == 'PRIMARY'


def test_functional_assay_with_b_annotation_requires_review():
    row, target = example()
    description(row, 'Potentiation of GABA current in human GABA-A alpha1beta3gamma2')
    row['standard_type'] = 'EC50'
    result = bf.classify_activity(row, target)
    assert result['endpoint_class'] == 'functional_potency'
    assert result['classification'] == 'REVIEW'
    assert 'assay_type_description_conflict' in result['review_flags']


def test_binding_ec50_is_ambiguous():
    row, target = example()
    row['standard_type'] = 'EC50'
    assert bf.classify_activity(row, target)['classification'] == 'REVIEW'


@pytest.mark.parametrize('relation', ['>', '<', '>=', '<=', '~'])
def test_censoring_retained_not_point_observation(relation):
    row, target = example()
    row.update(relation=relation, standard_relation=relation)
    result = bf.classify_activity(row, target)
    assert result['classification'] == 'PRIMARY'
    assert result['relation'] == result['standard_relation'] == relation
    assert result['standard_value'] == '12.00'
    assert not result['point_comparison_usable']


@pytest.mark.parametrize('field', ['standard_value', 'standard_units', 'standard_relation'])
def test_missing_not_zero_or_exact(field):
    row, target = example()
    row[field] = None
    result = bf.classify_activity(row, target)
    assert result[field] is None and result['classification'] == 'REVIEW'


def test_potential_duplicate_not_counted_independent():
    row, target = example()
    row['potential_duplicate'] = 1
    result = bf.classify_activity(row, target)
    assert result['classification'] == 'PRIMARY' and not result['point_comparison_usable']


def test_variant_not_wild_type():
    row, target = example()
    row['assay_variant_mutation'] = 'H101R'
    assert bf.classify_activity(row, target)['classification'] == 'REVIEW'


def test_nonbzd_readout_retained_with_exclude_reason():
    row, target = example()
    description(row, row['assay_description'].replace('flumazenil', 'TBPS'))
    result = bf.classify_activity(row, target)
    assert result['classification'] == 'EXCLUDE'
    assert result['exclude_reason'] == 'non_bzd_probe_or_channel_site_readout'


def test_names_never_control_classification():
    row, target = example(molecule='CHEMBL911')
    original = bf.classify_activity(row, target)
    row.update(drug_id='zolpidem', molecule_pref_name='Nonbenzodiazepine Z-drug')
    result = bf.classify_activity(row, target)
    assert result['classification'] == original['classification'] == 'PRIMARY'
    assert result['molecule_chembl_id'] == 'CHEMBL911'


def test_comparison_keys_independent_of_compound_but_not_assay():
    a, b = classified(activity=1), classified(activity=2, molecule='CHEMBL9999')
    groups, pairs = bf.build_comparisons([a, b])
    assert len(groups) == 1 and groups[0]['compound_count'] == 2
    assert groups[0]['activity_ids'] == [1, 2] and 'standard_value' not in groups[0]
    b['assay_chembl_id'] = 'CHEMBL_DIFFERENT'
    assert len(bf.build_comparisons([a, b])[0]) == 2


@pytest.mark.parametrize('field,value', [('document_chembl_id', 'CHEMBLOTHER'), ('standard_units', 'uM'),
    ('endpoint_class', 'functional_potency'), ('species', 'rat'), ('standard_relation', '>')])
def test_grouping_does_not_mix_context(field, value):
    a, b = classified(activity=1), classified(activity=2)
    b[field] = value
    assert len(bf.build_comparisons([a, b])[0]) == 2


def test_pairing_masks_only_alpha_not_beta():
    a, b = classified(activity=1, alpha=1, beta=3), classified(activity=2, alpha=2, beta=3)
    groups, pairs = bf.build_comparisons([a, b])
    assert len(groups) == 2 and len(pairs) == 1
    c = classified(activity=3, alpha=2, beta=2)
    assert not bf.build_comparisons([a, c])[1]
    b['document_chembl_id'] = 'CHEMBLOTHER'
    assert not bf.build_comparisons([a, b])[1]


@pytest.fixture
def real_data():
    source = bf.read_jsonl(SOURCE / 'processed/activities.jsonl')
    targets = {r['target_chembl_id']: r for r in bf.read_jsonl(SOURCE / 'processed/targets.jsonl')}
    with (SOURCE / 'config/compound_registry.csv').open(newline='') as handle:
        registry = list(csv.DictReader(handle))
    rows = [bf.classify_activity(a, targets[a['target_chembl_id']]) for a in source]
    return source, rows, registry


def test_real_434_complete_preserved_and_known_conflicts(real_data):
    source, rows, registry = real_data
    groups, pairs = bf.build_comparisons(rows)
    assert all(bf.validate_rows(source, rows, registry, groups, pairs, 434).values())
    lookup = {r['activity_id']: r for r in rows}
    for activity in [604850, 595754, 1876214, 1701379, 1701378]:
        assert lookup[activity]['classification'] == 'REVIEW'
        assert any(f.startswith('target_assay_') and f.endswith('_conflict') for f in lookup[activity]['review_flags'])
    for activity in [12662653, 12662654]:
        assert lookup[activity]['classification'] == 'REVIEW'
    coverage, candidates = bf.coverage_tables(rows, registry, pairs)
    assert len(coverage) == len(candidates) == 10
    assert len({r['molecule_chembl_id'] for r in coverage}) == 10
    assert {r['drug_id'] for r in pairs} == {'diazepam', 'alprazolam'}
    assert not any(r['additional_docking_coverage_rank'] for r in candidates if r['ifp_available'] == 'False')


def test_qc_rejects_source_tampering_and_duplicate_ids(real_data):
    source, rows, registry = real_data
    groups, pairs = bf.build_comparisons(rows)
    rows[0]['standard_relation'] = '>'
    with pytest.raises(ValueError, match='original_field'):
        bf.validate_rows(source, rows, registry, groups, pairs, 434)
    source[1]['activity_id'] = source[0]['activity_id']
    with pytest.raises(ValueError, match='unique_activity_ids'):
        bf.validate_rows(source, rows, registry, groups, pairs, 434)


def test_run_roundtrip_and_prior_files_immutable(tmp_path, monkeypatch):
    # Real acquired inputs, isolated output root. No API calls and no prior-run writes.
    original_root = bf.ROOT
    (tmp_path / 'runs').mkdir()
    (tmp_path / 'runs' / SOURCE.name).symlink_to(SOURCE, target_is_directory=True)
    pointer = tmp_path / 'runs/LATEST_RUN.txt'
    pointer.write_text('prior_run\n')
    original_source_hash = bf.file_inventory(SOURCE, SOURCE.rglob('*'))
    for name in ['src/bzd_pharmacology_filter.py', 'src/acquire_chembl.py', 'src/run_manager.py', 'tests/test_bzd_pharmacology_filter.py']:
        dest = tmp_path / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original_root / name, dest)
    config = tmp_path / 'config.json'
    shutil.copyfile(original_root / 'config/bzd_pharmacology_filter.json', config)
    monkeypatch.setattr(bf, 'ROOT', tmp_path)
    monkeypatch.setattr(bf.run_manager, 'ROOT', tmp_path)
    monkeypatch.setattr(bf.run_manager, 'RUNS', tmp_path / 'runs')
    monkeypatch.setattr(bf, 'tracked_inventory', lambda: {'runs/LATEST_RUN.txt': bf.sha(pointer.read_bytes())})
    first, qc = bf.run(config)
    assert qc['status'] == 'PASS'
    assert pointer.read_text() == 'prior_run\n'
    assert bf.file_inventory(SOURCE, SOURCE.rglob('*')) == original_source_hash
    manifest = json.loads((first / 'run_manifest.json').read_text())
    assert all(bf.sha((first / name).read_bytes()) == digest for name, digest in manifest['artifact_sha256'].items())
    assert len(bf.read_jsonl(first / 'processed/bzd_activity_classification.jsonl')) == 434
    first_hashes = bf.file_inventory(first, first.rglob('*'))
    second, _ = bf.run(config)
    assert second != first
    assert bf.file_inventory(first, first.rglob('*')) == first_hashes


def test_native_does_not_become_defined_recombinant_complex():
    row, target = example()
    description(row, row['assay_description'].replace('expressed in HEK293 cells', 'in native brain membrane'))
    result = bf.classify_activity(row, target)
    assert result['classification'] == 'SECONDARY'
    assert 'native_receptor_composition_not_unique' in result['classification_reasons']


def test_fixed_response_and_emax_not_same_group():
    row, target = example()
    row.update(standard_type='Efficacy', standard_units='%', assay_type='F', assay_detail_assay_type='F')
    description(row, 'Potentiation of GABA current in human GABA-A alpha1beta3gamma2 at 1 uM')
    fixed = bf.classify_activity(row, target)
    description(row, 'Maximal potentiation of GABA current in human GABA-A alpha1beta3gamma2')
    maximal = bf.classify_activity(row, target)
    maximal['activity_id'] = 2
    assert fixed['endpoint_detail'] == 'fixed_concentration_response'
    assert maximal['endpoint_detail'] == 'maximal_response'
    assert len(bf.build_comparisons([fixed, maximal])[0]) == 2


@pytest.mark.parametrize('value', ['NaN', 'Infinity', 'not-numeric', '0', '-1'])
def test_nonfinite_or_invalid_concentration_requires_review(value):
    row, target = example()
    row['standard_value'] = value
    result = bf.classify_activity(row, target)
    assert result['classification'] == 'REVIEW'
    assert result['standard_value'] == value and not result['point_comparison_usable']
