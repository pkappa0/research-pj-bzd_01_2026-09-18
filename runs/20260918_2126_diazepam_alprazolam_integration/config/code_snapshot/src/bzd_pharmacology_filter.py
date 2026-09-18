"""Conservative, auditable BZD-site pharmacology triage; no source data edits.

This is a deterministic candidate classifier, not an experimental adjudication.
Every source activity remains one row; conflicts take precedence over inclusion.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

from . import run_manager
from .acquire_chembl import verify_archive

ROOT = Path(__file__).resolve().parents[1]
CATEGORIES = ('PRIMARY', 'SECONDARY', 'REVIEW', 'EXCLUDE')
ENDPOINTS = ('binding', 'functional_potency', 'functional_efficacy', 'other')
SPECIES = {'Homo sapiens': 'human', 'Rattus norvegicus': 'rat', 'Mus musculus': 'mouse', 'Bos taurus': 'bovine'}
FAMILIES = ('alpha', 'beta', 'gamma', 'rho', 'delta', 'epsilon', 'theta', 'pi')
SUBUNIT_RE = re.compile(r'(alpha|beta|gamma|rho)[ -]?(\d+)([sl](?![a-z]))?|\b(delta|epsilon|theta|pi)\b')
BINDING_RE = re.compile(r'bind|affinity|displac|radioligand|radiolabeled|membrane filtration')
FUNCTION_RE = re.compile(r'current|potentiat|chloride influx|cl-? influx|voltage[- ]clamp|patch[- ]clamp|electrophysiolog|allosteric modulat|gaba.?(?:induced|elicited).*response')
BZD_PROBE_RE = re.compile(r'flumazenil|flunitraz(?:epam|apam)|\bro\W*15\W*1788\b|\bro\W*15\W*4513\b|\bcgs\W*8216\b|\bl\W*655708\b|\[(?:3h|11c)\][ -]*(?:diazepam|zolpidem)')
BZD_SITE_RE = re.compile(r'benzodiazepine|\bbzd\b|\bbz\b|\bbzr\b|\bcbr\b')
OFFSITE_RE = re.compile(r'\btbps\b|picrotoxin|muscimol|gabazine|bicuculline')
CONCENTRATION_TYPES = {'ki', 'kd', 'ic50', 'ec50', 'ac50', 'potency'}
CONCENTRATION_UNITS = {'M', 'mM', 'uM', 'µM', 'μM', 'nM', 'pM', 'fM'}
RESPONSE_TYPES = {'efficacy', 'emax', 'max activation', 'activity', 'inhibition', 'control', 'gaba current',
                  'change in cl- current', 'cl - current change', 'fc'}


def norm(value):
    text = str(value or '').lower().replace('−', '-').replace('–', '-').replace('—', '-')
    for greek, latin in [('α', 'alpha'), ('β', 'beta'), ('γ', 'gamma'), ('ρ', 'rho'), ('δ', 'delta')]:
        text = text.replace(greek, latin)
    return re.sub(r'\s+', ' ', text).strip()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def encode(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def key_for(prefix, value):
    return prefix + '_' + sha(encode(value).encode())[:20]


def read_jsonl(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def write_csv(path, rows, required=()):
    fields = list(dict.fromkeys([*required, *(k for r in rows for k in r)]))
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: encode(v) if isinstance(v, (dict, list)) else v for k, v in row.items()})


def extract_subunits(text):
    """Explicit tokens only; do not infer coassembly or missing subunits."""
    sequence, spans = [], []
    text = norm(text)
    for m in SUBUNIT_RE.finditer(text):
        token = (m[1] + m[2] + (m[3] or '')) if m[1] else m[4]
        sequence.append(token)
        spans.append(m.group())
    # Explicit gene symbols are annotations, not compound-name heuristics.
    for m in re.finditer(r'\bgabr([abg])(\d+)\b', text):
        sequence.append({'a': 'alpha', 'b': 'beta', 'g': 'gamma'}[m[1]] + m[2])
        spans.append(m.group())
    parts = {family: sorted({s for s in sequence if s.startswith(family)}) for family in FAMILIES}
    tokens = [s for family in FAMILIES for s in parts[family]]
    return {'parts': parts, 'sequence': sequence, 'evidence_tokens': spans,
            'composition': '_'.join(tokens) if tokens else None,
            'complete_single_abg': (len(parts['alpha']) == len(parts['beta']) == len(parts['gamma']) == 1
                                    and not any(parts[k] for k in ('rho', 'delta', 'epsilon', 'theta', 'pi')))}


def family_conflicts(left, right):
    return [family for family in FAMILIES if left['parts'][family] and right['parts'][family]
            and left['parts'][family] != right['parts'][family]]


def subtype_label(parts):
    alphas = set(parts['alpha'])
    if alphas == {'alpha1'}:
        return 'alpha1'
    if alphas == {'alpha2'}:
        return 'alpha2'
    if alphas == {'alpha1', 'alpha2'}:
        return 'alpha1_alpha2'
    return 'other_explicit_subtype' if alphas else 'unknown'


def species_label(organism):
    return SPECIES.get(organism, 'other' if organism else 'unknown')


def prose_species(text):
    # Do not assign host-cell species to the receptor (e.g. mouse LTK or Xenopus).
    receptor_text = re.split(r'\bexpressed in\b', norm(text), maxsplit=1)[0]
    result = []
    for pattern, organism in [(r'\bhuman\b', 'Homo sapiens'), (r'\brats?\b', 'Rattus norvegicus'),
                              (r'\b(?:mouse|mice|murine)\b', 'Mus musculus'), (r'\b(?:bovine|calf)\b', 'Bos taurus'),
                              (r'\bguinea[- ]pig\b', 'Cavia porcellus')]:
        if re.search(pattern, receptor_text):
            result.append(organism)
    return result


def endpoint_info(row):
    text, typ = norm(row.get('assay_detail_description')), norm(row.get('standard_type'))
    binding, functional = bool(BINDING_RE.search(text)), bool(FUNCTION_RE.search(text))
    flags, endpoint, detail, assay_class = [], 'other', 'unresolved', 'unknown'
    ratio = typ in {'ratio', 'gaba ratio', 'gaba shift', 'gs', 'gr', 'tbps shift', '[35s] tbps shift', 'shift'}
    if ratio:
        detail, assay_class = 'ratio_or_shift_not_affinity', 'binding' if binding or 'gaba' in text else 'unknown'
    elif typ in {'id50', 'ed50'}:
        detail, assay_class = 'whole_animal_dose_not_receptor_concentration', 'in_vivo_ex_vivo'
    elif typ in {'ki', 'kd'}:
        endpoint, detail, assay_class = 'binding', typ, 'binding'
        if functional:
            flags.append('affinity_endpoint_with_functional_description')
        if not binding:
            flags.append('affinity_endpoint_without_binding_description')
    elif typ in {'ic50', 'ac50', 'ec50', 'potency', 'log ic50'}:
        if binding and functional:
            flags.append('binding_functional_description_ambiguous')
        elif binding and typ != 'ec50':
            endpoint, detail, assay_class = 'binding', 'binding_' + typ, 'binding'
        elif functional:
            endpoint, detail, assay_class = 'functional_potency', typ, 'functional'
        else:
            flags.append('concentration_endpoint_context_ambiguous')
        if typ == 'log ic50':
            flags.append('log_endpoint_scale_requires_review')
    elif functional and typ in RESPONSE_TYPES:
        endpoint, assay_class = 'functional_efficacy', 'functional'
        detail = ('maximal_response' if re.search(r'maximal|maximum|\bemax\b|max activation', text + ' ' + typ)
                  else 'fixed_concentration_response' if re.search(r'(?:at|concentration of) \d', text)
                  else 'response_not_established_emax')
        if binding:
            flags.append('binding_functional_description_ambiguous')
    elif binding and typ in {'displacement', 'inhibition', 'activity', 'binding affinity', 'control'}:
        endpoint, detail, assay_class = 'binding', 'binding_' + typ, 'binding'
    else:
        flags.append('endpoint_interpretation_ambiguous')
    raw_class = row.get('assay_type')
    if assay_class == 'binding' and raw_class == 'F':
        flags.append('assay_type_description_conflict')
    if assay_class == 'functional' and raw_class in {'B', 'A', 'P', 'T'}:
        flags.append('assay_type_description_conflict')
    if row.get('assay_detail_assay_type') != raw_class:
        flags.append('assay_activity_type_conflict')
    # Concentration-valued IC50 is not the same thing as percent inhibition at one dose.
    if typ in CONCENTRATION_TYPES and re.search(r'percent inhibition.*\bat \d', text):
        flags.append('concentration_endpoint_vs_fixed_percent_description')
    return endpoint, detail, assay_class, flags


def classify_activity(source, target):
    row = dict(source)  # Every original column/value survives, including null and raw relations.
    desc = source.get('assay_detail_description')
    text = norm(desc)
    target_name = target.get('pref_name')
    target_type = target.get('target_type')
    tn, ad = extract_subunits(target_name), extract_subunits(desc)
    components = target.get('target_components') or []
    component_text = ' '.join(c.get('component_description') or '' for c in components)
    tc = extract_subunits(component_text)
    direct_components = target_type in {'PROTEIN COMPLEX', 'SINGLE PROTEIN'} and all(
        c.get('relationship') in {'PROTEIN SUBUNIT', 'SINGLE PROTEIN'} for c in components)
    review, notes = [], []
    for field, expected in [('target_chembl_id', target.get('target_chembl_id')),
                            ('assay_detail_target_chembl_id', source.get('target_chembl_id')),
                            ('assay_detail_assay_chembl_id', source.get('assay_chembl_id')),
                            ('assay_detail_document_chembl_id', source.get('document_chembl_id'))]:
        if source.get(field) != expected:
            review.append(field + '_inconsistent')
    if source.get('target_organism') != target.get('organism'):
        review.append('activity_target_metadata_organism_conflict')
    conflicts = family_conflicts(tn, ad)
    review += ['target_assay_' + f + '_conflict' for f in conflicts]
    if direct_components:
        review += ['target_name_component_' + f + '_conflict' for f in family_conflicts(tn, tc)]
        review += ['target_component_assay_' + f + '_conflict' for f in family_conflicts(tc, ad)]
    else:
        notes.append('component_group_members_not_treated_as_coassembled_subunits')
    target_species = source.get('target_organism')
    assay_species = source.get('assay_detail_assay_organism')
    explicit_species = prose_species(desc)
    species_conflict = bool(target_species and assay_species and target_species != assay_species)
    if species_conflict:
        review.append('target_assay_organism_conflict')
    if explicit_species and any(s != target_species for s in explicit_species if target_species):
        review.append('target_description_organism_conflict')
        species_conflict = True
    if explicit_species and any(s != assay_species for s in explicit_species if assay_species):
        review.append('assay_metadata_description_organism_conflict')
        species_conflict = True
    if not assay_species:
        notes.append('assay_organism_missing_not_filled')
    species = 'unknown' if species_conflict else species_label(assay_species)
    if re.search(r'unknown origin|not specified|not tested|not determined', text):
        review.append('assay_origin_or_measurement_uncertain')
    if not desc:
        review.append('missing_assay_description')
    if source.get('assay_description') != desc:
        review.append('activity_assay_description_conflict')
    if re.search(r'\bmutant\b|\bmutat\w*|\bchimera\w*', text) or source.get('assay_variant_mutation') or source.get('assay_detail_variant_sequence'):
        review.append('variant_construct_requires_review')
    if len(ad['sequence']) != len(set(ad['sequence'])):
        review.append('repeated_subunit_sequence_or_concatemer_review')
    if re.search(r'\btimes ki\b', text):
        review.append('compound_relative_test_concentration')
    if source.get('data_validity_comment'):
        review.append('chembl_data_validity_comment')
    endpoint, endpoint_detail, assay_class, endpoint_flags = endpoint_info(source)
    review.extend(endpoint_flags)
    numeric = None
    try:
        if source.get('standard_value') not in (None, ''):
            numeric = Decimal(str(source['standard_value']))
            if not numeric.is_finite():
                review.append('nonfinite_standard_value')
        else:
            review.append('missing_standard_value')
    except InvalidOperation:
        review.append('nonnumeric_standard_value')
    typ, unit, relation = norm(source.get('standard_type')), source.get('standard_units'), source.get('standard_relation')
    if typ in CONCENTRATION_TYPES and unit not in CONCENTRATION_UNITS:
        review.append('missing_or_incompatible_concentration_unit')
    if typ in CONCENTRATION_TYPES and numeric is not None and numeric.is_finite() and numeric <= 0:
        review.append('nonpositive_concentration')
    if endpoint in {'binding', 'functional_potency', 'functional_efficacy'} and not unit:
        review.append('missing_endpoint_unit')
    if relation not in {'=', '<', '>', '<=', '>=', '~'}:
        review.append('missing_or_unrecognised_relation')
    if relation and relation != '=':
        notes.append('censored_or_nonexact_value_retained')
    if source.get('potential_duplicate'):
        notes.append('chembl_potential_duplicate_not_independent_replication')
    if not source.get('document_chembl_id') or not source.get('assay_chembl_id'):
        review.append('missing_assay_or_document_id')
    gtext = norm(target_name) + ' ' + text
    gabaa = bool(re.search(r'gaba|aminobutyric', gtext)) and not re.search(r'gaba[- ]?b\b|translocator', norm(target_name))
    bzd_text = bool(BZD_SITE_RE.search(text) or BZD_PROBE_RE.search(text))
    offsite = bool(OFFSITE_RE.search(text))
    ratio = endpoint_detail == 'ratio_or_shift_not_affinity'
    evidence = ('non_bzd_probe_readout' if offsite else 'explicit_site_or_probe' if bzd_text else
                'explicit_abg_receptor_candidate' if ad['complete_single_abg'] and ad['parts']['gamma'] == ['gamma2'] else
                'gaba_shift_indirect' if ratio and 'gaba' in text else 'unconfirmed')
    # No merged/union composition: an accepted composition must be fully stated in assay prose.
    composition = None
    composition_source = None
    if not conflicts and ad['composition']:
        composition, composition_source = ad['composition'], 'assay_description_explicit_tokens'
    elif conflicts:
        notes.append('composition_conflict_not_resolved')
    subtype = subtype_label(ad['parts']) if not conflicts else 'unknown'
    if not ad['parts']['alpha']:
        subtype = 'unknown'
    observed_alpha = sorted(set(tn['parts']['alpha'] + ad['parts']['alpha']))
    complete_agreement = bool(tn['complete_single_abg'] and ad['complete_single_abg'] and tn['composition'] == ad['composition'])
    native = bool(re.search(r'\bnative\b|cerebellar granule|brain membrane|cerebral cortex|brain homogenate', text)) and 'recombinant' not in text
    in_vivo = bool(re.search(r'\bin vivo\b|\bex vivo\b|intraperitoneal|convulsion|after (?:ip|po) administration', text)) or typ in {'id50', 'ed50'}
    if evidence == 'unconfirmed' and not offsite and not in_vivo:
        review.append('bzd_site_relevance_unconfirmed')
    # Evidence-specific rules never depend on compound name or hard-coded activity IDs.
    review = sorted(set(review))
    reasons = []
    if review:
        category, reasons = 'REVIEW', review.copy()
    elif not gabaa:
        category, reasons = 'EXCLUDE', ['outside_gabaa_scope']
    elif in_vivo:
        category, reasons = 'EXCLUDE', ['whole_animal_or_ex_vivo_not_receptor_pharmacology']
    elif offsite:
        category, reasons = 'EXCLUDE', ['non_bzd_probe_or_channel_site_readout']
    elif (complete_agreement and target_type == 'PROTEIN COMPLEX' and subtype in {'alpha1', 'alpha2'}
          and ad['parts']['gamma'] == ['gamma2'] and not native and endpoint != 'other'):
        category, reasons = 'PRIMARY', ['explicit_alpha1_or_alpha2_abg2_in_both_sources', 'no_detected_consistency_conflict']
    else:
        category = 'SECONDARY'
        if target_type != 'PROTEIN COMPLEX':
            reasons.append('generic_or_single_subunit_target')
        if not complete_agreement:
            reasons.append('full_composition_not_independently_explicit_in_both_sources')
        if ad['parts']['gamma'] != ['gamma2']:
            reasons.append('gamma2_not_explicit_in_assay')
        if subtype == 'other_explicit_subtype':
            reasons.append('other_subtype_context_outside_alpha1_alpha2_primary')
        if subtype == 'alpha1_alpha2' or len(ad['parts']['alpha']) > 1:
            reasons.append('mixed_alpha_composition_not_allocated_to_single_subtypes')
        if native:
            reasons.append('native_receptor_composition_not_unique')
        if endpoint == 'other':
            reasons.append('indirect_ratio_or_other_endpoint')
        if not reasons:
            reasons.append('insufficient_direct_comparison_context')
    point_issues = []
    if category != 'PRIMARY': point_issues.append('not_primary')
    if relation != '=': point_issues.append('not_exact_relation')
    if numeric is None or not numeric.is_finite(): point_issues.append('no_finite_value')
    if not unit: point_issues.append('no_reported_unit')
    if species == 'unknown': point_issues.append('species_not_resolved_from_assay_metadata')
    if source.get('potential_duplicate'): point_issues.append('potential_duplicate')
    row.update({
        'classification': category, 'subtype': subtype, 'endpoint_class': endpoint,
        'endpoint_detail': endpoint_detail, 'derived_assay_class': assay_class,
        'species': species, 'species_basis': 'assay_organism' if species != 'unknown' else 'unresolved_no_imputation',
        'target_species': species_label(target_species), 'assay_species': species_label(assay_species),
        'description_organisms_explicit': explicit_species,
        'species_priority': 1 if species == 'human' else 2 if species in {'rat', 'mouse'} else 4 if species == 'unknown' else 3,
        'receptor_composition': composition, 'composition_source': composition_source,
        'target_name_composition': tn['composition'], 'target_component_composition': tc['composition'],
        'target_components_are_coassembly_evidence': direct_components,
        'assay_description_composition': ad['composition'], 'assay_subunit_sequence': ad['sequence'],
        'target_name_subtype': subtype_label(tn['parts']), 'assay_description_subtype': subtype_label(ad['parts']),
        'observed_alpha_mentions_for_review_only': observed_alpha,
        'bzd_relevance_evidence': evidence, 'classification_reasons': reasons,
        'inclusion_reason': ';'.join(reasons) if category in {'PRIMARY', 'SECONDARY'} else None,
        'exclude_reason': ';'.join(reasons) if category == 'EXCLUDE' else None,
        'manual_review_required': category == 'REVIEW', 'review_flags': review, 'quality_notes': sorted(set(notes)),
        'point_comparison_usable': not point_issues, 'point_comparison_limitations': point_issues,
        'target_components_used': [{k: c.get(k) for k in ['component_id', 'accession', 'component_description', 'relationship']} for c in components],
        'comparison_group_id': None, 'comparison_key_json': None,
    })
    for label, parsed in [('target_name', tn), ('target_components', tc), ('assay_description', ad)]:
        for family in ('alpha', 'beta', 'gamma'):
            row[label + '_' + family + '_subunits'] = parsed['parts'][family]
    return row


def primary_group_key(row):
    """Compound-independent, but strictly assay- and document-specific."""
    return {k: row.get(k) for k in [
        'species', 'assay_detail_assay_organism', 'receptor_composition', 'derived_assay_class',
        'endpoint_class', 'endpoint_detail', 'standard_type', 'standard_units', 'standard_relation',
        'type', 'units', 'target_chembl_id', 'assay_chembl_id', 'document_chembl_id',
        'assay_detail_description', 'assay_detail_assay_parameters', 'assay_detail_variant_sequence']}


def subtype_mask(text):
    return re.sub(r'alpha[ -]?[12](?!\d)', 'alphaX', norm(text))


def pairing_key(row):
    """Candidate cross-subtype context; only alpha1/2 is masked, never beta/gamma."""
    if row['subtype'] not in {'alpha1', 'alpha2'}:
        return None
    description = subtype_mask(row.get('assay_detail_description'))
    context = {k: row.get(k) for k in [
        'species', 'assay_detail_assay_organism', 'derived_assay_class', 'endpoint_class', 'endpoint_detail',
        'standard_type', 'standard_units', 'standard_relation', 'type', 'units', 'document_chembl_id',
        'assay_type', 'assay_detail_assay_cell_type', 'assay_detail_assay_tissue', 'assay_detail_assay_strain',
        'assay_detail_assay_subcellular_fraction', 'assay_detail_assay_parameters',
        'assay_detail_bao_format', 'assay_detail_variant_sequence', 'assay_detail_relationship_type']}
    context.update(receptor_background=subtype_mask(row['receptor_composition']),
                   subtype_masked_description=description)
    return context


def build_comparisons(rows):
    groups, paired = defaultdict(list), defaultdict(list)
    for row in rows:
        if row['classification'] != 'PRIMARY':
            continue
        key = primary_group_key(row)
        row['comparison_key_json'] = key
        row['comparison_group_id'] = key_for('CG', key)
        groups[row['comparison_group_id']].append(row)
        if row['point_comparison_usable']:
            paired[(row['molecule_chembl_id'], key_for('PAIR', pairing_key(row)))].append(row)
    group_rows, pair_rows = [], []
    for gid, members in sorted(groups.items()):
        first = members[0]
        group_rows.append(dict(first['comparison_key_json'], comparison_group_id=gid,
            activity_count=len(members), compound_count=len({r['molecule_chembl_id'] for r in members}),
            activity_ids=[r['activity_id'] for r in members], compound_ids=sorted({r['molecule_chembl_id'] for r in members}),
            drug_ids=sorted({r['drug_id'] for r in members}), point_usable_count=sum(r['point_comparison_usable'] for r in members),
            numeric_values_not_aggregated=True))
    for (mid, pid), members in sorted(paired.items()):
        alpha1 = [r for r in members if r['subtype'] == 'alpha1']
        alpha2 = [r for r in members if r['subtype'] == 'alpha2']
        if not alpha1 or not alpha2:
            continue
        pair_rows.append({'pair_candidate_id': pid, 'drug_id': members[0]['drug_id'], 'molecule_chembl_id': mid,
            'recorded_context_key': pairing_key(members[0]), 'alpha1_activity_ids': [r['activity_id'] for r in alpha1],
            'alpha2_activity_ids': [r['activity_id'] for r in alpha2],
            'alpha1_composition': alpha1[0]['receptor_composition'], 'alpha2_composition': alpha2[0]['receptor_composition'],
            'species': members[0]['species'], 'endpoint_class': members[0]['endpoint_class'],
            'standard_type': members[0]['standard_type'], 'standard_units': members[0]['standard_units'],
            'document_chembl_id': members[0]['document_chembl_id'],
            'alpha1_assay_ids': sorted({r['assay_chembl_id'] for r in alpha1}),
            'alpha2_assay_ids': sorted({r['assay_chembl_id'] for r in alpha2}),
            'status': 'recorded_context_candidate_unreported_conditions_not_verified'})
    return group_rows, pair_rows


def coverage_tables(rows, registry, pair_rows):
    coverage, candidates = [], []
    for compound in registry:
        mid = compound['molecule_chembl_id']
        members = [r for r in rows if r['molecule_chembl_id'] == mid]
        primary = [r for r in members if r['classification'] == 'PRIMARY']
        usable = [r for r in primary if r['point_comparison_usable']]
        review = [r for r in members if r['classification'] == 'REVIEW']
        cell = {k: compound[k] for k in ['drug_id', 'molecule_chembl_id', 'ifp_available']}
        for sub in ('alpha1', 'alpha2'):
            for ep in ENDPOINTS[:3]:
                subset = [r for r in primary if r['subtype'] == sub and r['endpoint_class'] == ep]
                numeric = [r for r in subset if r['point_comparison_usable']]
                uncertain = [r for r in review if sub in r['observed_alpha_mentions_for_review_only'] and r['endpoint_class'] == ep]
                prefix = sub + '_' + ep
                cell.update({prefix + '_primary_count': len(subset), prefix + '_usable_activity_count': len(numeric),
                    prefix + '_comparison_group_count': len({r['comparison_group_id'] for r in numeric}),
                    prefix + '_human_data': any(r['species'] == 'human' for r in numeric),
                    prefix + '_manual_review_count': len(uncertain), prefix + '_manual_review_present': bool(uncertain)})
        pairs = [r for r in pair_rows if r['molecule_chembl_id'] == mid]
        candidate = dict({k: compound[k] for k in ['drug_id', 'molecule_chembl_id', 'ifp_available']},
            alpha1_primary_data=any(r['subtype'] == 'alpha1' for r in primary),
            alpha2_primary_data=any(r['subtype'] == 'alpha2' for r in primary),
            alpha1_usable_count=sum(r['subtype'] == 'alpha1' for r in usable),
            alpha2_usable_count=sum(r['subtype'] == 'alpha2' for r in usable),
            binding_data=any(r['endpoint_class'] == 'binding' for r in usable),
            functional_data=any(r['endpoint_class'].startswith('functional_') for r in usable),
            human_data=any(r['species'] == 'human' for r in usable),
            primary_activity_count=len(primary), comparable_activity_count=len(usable),
            review_activity_count=len(review), review_flag_count=sum(len(r['review_flags']) for r in review),
            recorded_context_pair_count=len(pairs),
            both_subtypes_present=all(any(r['subtype'] == sub for r in primary) for sub in ('alpha1', 'alpha2')),
            both_subtypes_point_usable=all(any(r['subtype'] == sub for r in usable) for sub in ('alpha1', 'alpha2')),
            reviewed_alpha1_mentions=sum('alpha1' in r['observed_alpha_mentions_for_review_only'] for r in review),
            reviewed_alpha2_mentions=sum('alpha2' in r['observed_alpha_mentions_for_review_only'] for r in review),
            category_counts=dict(Counter(r['classification'] for r in members)))
        candidate['docking_triage'] = ('existing_ifp' if str(compound['ifp_available']).lower() == 'true' else
            'data_supported_candidate_pending_methods_check' if pairs else
            'both_subtypes_but_no_matched_recorded_context' if candidate['both_subtypes_point_usable'] else
            'manual_review_before_docking_priority' if candidate['reviewed_alpha1_mentions'] and candidate['reviewed_alpha2_mentions'] else
            'insufficient_subtype_data')
        coverage.append(dict(cell, both_subtypes_present=candidate['both_subtypes_present'],
                             both_subtypes_point_usable=candidate['both_subtypes_point_usable'],
                             recorded_context_pair_count=len(pairs), total_review_activity_count=len(review)))
        candidates.append(candidate)
    # Rank only records with qualifying pairs. No artificial winner among unsupported compounds.
    eligible = [r for r in candidates if r['docking_triage'] == 'data_supported_candidate_pending_methods_check']
    eligible.sort(key=lambda r: (-int(r['human_data']), -r['recorded_context_pair_count'],
                                 -r['comparable_activity_count'], r['review_activity_count'], r['molecule_chembl_id']))
    rank = {r['molecule_chembl_id']: i + 1 for i, r in enumerate(eligible)}
    for r in candidates:
        r['additional_docking_coverage_rank'] = rank.get(r['molecule_chembl_id'])
    return coverage, candidates


def validate_rows(source, rows, registry, groups, pairs, expected):
    """Fail closed on lost source values, row duplication, or incompatible grouping."""
    ids = [r['activity_id'] for r in source]
    checks = {
        'expected_source_count': len(source) == expected,
        'activity_id_one_to_one_in_order': ids == [r['activity_id'] for r in rows],
        'unique_activity_ids': len(ids) == len(set(ids)) == len(rows),
        'all_rows_exactly_one_category': all(r['classification'] in CATEGORIES for r in rows),
        'category_total_equals_source': sum(Counter(r['classification'] for r in rows).values()) == expected,
        'every_original_field_preserved_including_null_relation_and_ids': all(
            all(k in dst and dst[k] == value for k, value in src.items()) for src, dst in zip(source, rows)),
        'no_review_in_primary': all(not r['review_flags'] and not r['manual_review_required']
                                    for r in rows if r['classification'] == 'PRIMARY'),
        'registry_ids_unique_and_exactly_cover_source': len(registry) == len({r['molecule_chembl_id'] for r in registry})
            and {r['molecule_chembl_id'] for r in rows} == {r['molecule_chembl_id'] for r in registry},
        'species_not_imputed': all(r['species'] == 'unknown' or
            r['species'] == species_label(r.get('assay_detail_assay_organism')) for r in rows)
            and all(r['species'] == 'unknown' for r in rows if not r.get('assay_detail_assay_organism')),
        'composition_not_imputed': all(r['receptor_composition'] is None or
            r['receptor_composition'] == extract_subunits(r.get('assay_detail_description'))['composition'] for r in rows),
        'censored_values_never_point_usable': all(not r['point_comparison_usable'] for r in rows if r.get('standard_relation') != '='),
    }
    indexed = {r['activity_id']: r for r in rows}
    checks['group_keys_keep_composition_endpoint_species_unit_document_assay_separate'] = all(
        all(primary_group_key(indexed[i]) == primary_group_key(indexed[g['activity_ids'][0]]) for i in g['activity_ids'])
        for g in groups)
    group_members = [i for g in groups for i in g['activity_ids']]
    checks['groups_cover_primary_once'] = sorted(group_members) == sorted(r['activity_id'] for r in rows if r['classification'] == 'PRIMARY')
    checks['pairs_match_recorded_context_and_beta_gamma'] = all(
        all(indexed[i]['point_comparison_usable'] and pairing_key(indexed[i]) == p['recorded_context_key']
            and indexed[i]['molecule_chembl_id'] == p['molecule_chembl_id']
            for i in p['alpha1_activity_ids'] + p['alpha2_activity_ids']) for p in pairs)
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError('QC failed: ' + ', '.join(failed))
    return checks


def file_inventory(root, paths):
    return {str(p.relative_to(root)): sha(p.read_bytes()) for p in sorted(paths) if p.is_file()}


def tracked_inventory():
    names = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().split('\0')
    return file_inventory(ROOT, [ROOT / name for name in names if name])


def tests_summary(path):
    if path is None:
        return {'status': 'NOT_PROVIDED'}
    tree = ET.parse(path).getroot()
    suites = [tree] if tree.tag == 'testsuite' else list(tree.iter('testsuite'))
    result = {key: sum(int(s.get(key, 0)) for s in suites) for key in ['tests', 'failures', 'errors', 'skipped']}
    if not result['tests'] or result['failures'] or result['errors']:
        raise ValueError('Provided test results are empty or failing')
    return dict(result, status='PASSED', source_sha256=sha(path.read_bytes()))


def make_report(rows, groups, pairs, candidates, qc, source_run, run_id):
    primary = [r for r in rows if r['classification'] == 'PRIMARY']
    def counts(items, field):
        return encode(dict(sorted(Counter(r[field] for r in items).items())))
    flags = Counter(f for r in rows for f in r['review_flags'])
    names = lambda test: ', '.join(c['drug_id'] for c in candidates if test(c)) or 'なし'
    lines = [f'# BZD-site pharmacology filter — {run_id}', '',
        f'入力: `{source_run}/processed/activities.jsonl`（ChEMBL 37、434 activity、10剤）。',
        '取得済みアーカイブだけを使用。追加API取得・ML・Docking・phenotype予測は実施していない。', '',
        '## 結果', '', '| category | activity数 |', '|---|---:|']
    total = Counter(r['classification'] for r in rows)
    lines += [f'| {cat} | {total[cat]} |' for cat in CATEGORIES]
    lines += [f'| 合計 | {sum(total.values())} |', '',
        '434行すべてを1カテゴリへ分類し、元activity IDと1対1対応。PRIMARYは直接比較の**候補**であり、独立した検証済み教師ラベルではない。',
        'PRIMARYには数値を保持したcensored値、potential_duplicate、assay organism欠損も含む。厳しいpoint比較条件を満たす件数を別に示す。', '',
        f'- PRIMARY subtype: `{counts(primary, "subtype")}`。alpha1_alpha2 / other_explicit_subtypeはPRIMARY 0件。',
        f'- 全件endpoint: `{counts(rows, "endpoint_class")}`',
        f'- PRIMARY endpoint: `{counts(primary, "endpoint_class")}`',
        f'- 全件species: `{counts(rows, "species")}`',
        f'- PRIMARY species: `{counts(primary, "species")}`',
        f'- PRIMARY point比較候補: {sum(r["point_comparison_usable"] for r in primary)}件、comparison group: {len(groups)}件。', '',
        '## Compound coverage', '',
        'usableはPRIMARYかつ有限値・relation `=`・報告単位あり・assay species既知・potential_duplicateなし。',
        '各α1/α2 × binding/potency/efficacyの件数、group数、human有無、review有無は `compound_pharmacology_coverage.csv` に保存。', '',
        '| compound | ChEMBL ID | IFP済 | PRIMARY | α1 usable | α2 usable | binding | functional | human | REVIEW行 | flags延べ数 | 記録条件pair数 |',
        '|---|---|---|---:|---:|---:|---|---|---|---:|---:|---:|']
    for c in candidates:
        lines.append('| ' + ' | '.join(str(c[k]) for k in [
            'drug_id', 'molecule_chembl_id', 'ifp_available', 'primary_activity_count', 'alpha1_usable_count', 'alpha2_usable_count',
            'binding_data', 'functional_data', 'human_data', 'review_activity_count', 'review_flag_count', 'recorded_context_pair_count']) + ' |')
    lines += ['', f'α1/α2双方にPRIMARYあり: **{names(lambda c: c["both_subtypes_present"])}**。',
        f'α1/α2双方にpoint比較候補あり: **{names(lambda c: c["both_subtypes_point_usable"])}**。',
        f'同一compound・文献・species・β/γ・endpoint・単位・記録条件のα1/α2 pair候補あり: **{names(lambda c: c["recorded_context_pair_count"] > 0)}**。', '',
        'pair候補はassay IDをまたぐため、α1/α2表記だけをマスクしたdescriptionと、保存されたcell/tissue/strain/parameter等の完全一致を要求する。',
        '未記載条件の一致を証明するものではない。descriptionの表記差で見逃す可能性があり、pairなしは生物学的比較不能の証明ではない。',
        'comparison group自体はassay ID・文献・species・full composition・type・unit・relationごと。compound名はkeyに含めず、平均・代表値・比は計算しない。', '',
        '## 既存IFP 4剤と追加6剤', '',
        'diazepam・alprazolam・triazolam・zolpidemは既存IFP群として保持。triazolamのPRIMARYはassay organism欠損のためpoint比較候補に数えない。',
        'zolpidemのα1/α2 usable行が存在しても、厳密な記録条件pairが成立するとは限らない。',
        'alprazolamのpair候補はβ2背景を含み、既存IFPのβ3背景と同一構成とは扱わない。pharmacology内のpairとIFP構造との一致は別途確認が必要。',
        'lorazepam・clonazepam・midazolam・temazepam・zopiclone・zaleplonは未IFP群。骨格・compound名による除外は行っていない。',
        f'現在のルールで追加Dockingを検討するための記録条件pairを持つ未解析compound: **{names(lambda c: c["docking_triage"] == "data_supported_candidate_pending_methods_check")}**。',
        '該当しない薬剤に順位や科学的best drugを付与しない。`additional_docking_coverage_rank` の欠損は評価不十分を意味する。',
        'zaleplonは原著確認の候補。current modulationのEC50なのにassay_type BとされるためREVIEW。またα1β2γ2とα2β3γ2はβ背景が異なり、注釈を確認してもそのままmatched pairにはできない。',
        'zopicloneのnative receptor/current efficacyも構成とassay分類の確認が必要。残る4剤はsubtype-completeな比較データの追加探索が先となる。', '',
        '## Major review issues', '', '| review flag | activity数（重複計上あり） |', '|---|---:|']
    lines += [f'| {name} | {count} |' for name, count in flags.most_common()]
    lines += ['', 'α1対α4、α1対α5、β1対β2の矛盾は一般的なfamily別集合比較で検出。activity IDに依存した分類分岐はない。',
        'target/assay species不一致、host organism由来の可能性、unknown originの記述も自動修正しない。',
        f'assay organism欠損は全件中{sum(not r.get("assay_detail_assay_organism") for r in rows)}件。target/proseから補完せずunknownとする。矛盾時も解析用speciesはunknown、各元speciesは保持。',
        'B/F assay分類は粗い注釈なので、不一致は実験が誤りという断定ではなく原著確認要求。REVIEW優先のため、off-siteと判断できる行でも別の矛盾があればREVIEWに残る。',
        'component metadataのGROUP MEMBERは受容体の共集合構成と解釈しない。single-subunit targetやgeneric/nativeはPRIMARYに昇格させない。',
        'functional_efficacyは広いresponse分類。固定濃度応答と真のEmaxはendpoint_detailで区別し、同じgroupに入れない。pChEMBLへの統一やKi/IC50/EC50/Emaxの数値統合は行わない。', '',
        '## QC・provenance・再現', '',
        f'QC: **{qc["status"]}**。{len(qc["checks"])}チェックを通過。入力raw gzipと元manifestのartifact hashも検証。',
        f'tests: `{encode(qc["tests"])}`。詳細は `logs/pytest.xml`。',
        '`reports/QC_REPORT.json` はID一意性、元全fieldの完全保持、null・relation・censoring・species・composition・group整合性を記録。',
        '`raw/lineage/` は入力processed JSONL、registry、元manifest、raw request ledgerのbyte-identical snapshot。',
        '元APIレスポンスはsource runの `raw/responses/` に保持し、コピーし直さない。元source_file / JSON pointerは **source run基準**。',
        '`processed/bzd_activity_classification.jsonl` はnullと元数値文字列を保持する正本。CSVの空セルは欠損（0ではない）、list/dictはJSON表現。',
        '各activityへsource run、元processed行番号・ファイルhash、rule versionを追記。元target/component/assay/document情報とURL・取得時刻・raw hashは失わない。',
        '実行時コード、config、tests、依存関数を `config/code_snapshot/` に保存。manifestのgit_commitは実行前HEADであり、実行コードの厳密な識別はsnapshot SHA-256による。',
        '既存tracked filesは実行前後のhash一致を確認。旧run、IFP、phenotype、`runs/LATEST_RUN.txt` を更新しない。', '',
        '再実行: `MPLBACKEND=Agg .venv/bin/python -m pytest -q --junitxml=/tmp/bzd-tests.xml` の後、',
        '` .venv/bin/python -m src.bzd_pharmacology_filter --config config/bzd_pharmacology_filter.json --test-results /tmp/bzd-tests.xml`。',
        '再実行は別runを生成する。ルールは保守的な機械抽出であり、negationや複雑なconstructを完全には解釈できない。PRIMARYも原著とIFP構造のspecies/subunit/条件を照合してから使用する。', '',
        '## 判定設計の参照', '',
        '- [ChEMBL data FAQ](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions): assay分類、target注釈、potential_duplicateを別々の情報として扱う。',
        '- [BZD-site構造研究](https://www.nature.com/articles/s41586-018-0255-3)、[GABAA構造薬理](https://www.nature.com/articles/s41586-018-0832-5): α/γ界面と明示構成を重視。',
        '- [TBPS binding study](https://pubmed.ncbi.nlm.nih.gov/3035434/): channel-site readoutをBZD-site affinityと混同しない。', '']
    return '\n'.join(lines)


def run(config_path, test_results=None):
    config = json.loads(config_path.read_text())
    source_dir = ROOT / config['source_run']
    raw_count = verify_archive(source_dir)
    before = tracked_inventory()
    inputs = ['processed/activities.jsonl', 'processed/targets.jsonl', 'processed/assays.jsonl',
              'processed/documents.jsonl', 'processed/compounds.jsonl', 'config/compound_registry.csv',
              'run_manifest.json', 'raw/requests.jsonl']
    input_hashes = file_inventory(source_dir, [source_dir / name for name in inputs])
    if len(input_hashes) != len(inputs):
        raise ValueError('Missing required source inputs')
    activities = read_jsonl(source_dir / 'processed/activities.jsonl')
    targets_list = read_jsonl(source_dir / 'processed/targets.jsonl')
    targets = {t['target_chembl_id']: t for t in targets_list}
    if len(targets) != len(targets_list):
        raise ValueError('Duplicate target IDs')
    with (source_dir / 'config/compound_registry.csv').open(newline='') as handle:
        registry = list(csv.DictReader(handle))
    rows = []
    for index, a in enumerate(activities, 1):
        row = classify_activity(a, targets[a['target_chembl_id']])
        row.update(classification_source_run=config['source_run'], classification_rule_version=config['rule_version'],
                   source_processed_line=index, source_processed_sha256=input_hashes['processed/activities.jsonl'])
        rows.append(row)
    groups, pairs = build_comparisons(rows)
    coverage, candidates = coverage_tables(rows, registry, pairs)
    checks = validate_rows(activities, rows, registry, groups, pairs, config['expected_activity_count'])
    test_summary = tests_summary(test_results)
    code_paths = [ROOT / name for name in ['src/bzd_pharmacology_filter.py', 'src/acquire_chembl.py', 'src/run_manager.py',
                                          'tests/test_bzd_pharmacology_filter.py']]
    run_dir, manifest = run_manager.create_run(short_phase_name='bzd_pharmacology_filter',
        analysis_phase='BZD-site subtype pharmacology candidate triage',
        analysis_purpose='Classify acquired activities without imputation, averaging, or source modification',
        parent_run=source_dir.name, source_run=config['source_run'],
        drugs=[c['drug_id'] for c in registry], receptors={}, docking_parameters={},
        input_files=[str((source_dir / name).relative_to(ROOT)) for name in inputs],
        notes=['Offline classification only', 'LATEST_RUN.txt intentionally unchanged', 'Not validated ML labels'], code_paths=code_paths)
    try:
        for name in inputs:
            dest = run_dir / 'raw/lineage' / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_dir / name, dest)
        shutil.copyfile(config_path, run_dir / 'config/bzd_pharmacology_filter.json')
        for path in code_paths:
            dest = run_dir / 'config/code_snapshot' / path.relative_to(ROOT)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, dest)
        tables = {
            'bzd_activity_classification': rows,
            'primary_bzd_pharmacology': [r for r in rows if r['classification'] == 'PRIMARY'],
            'bzd_manual_review': [r for r in rows if r['classification'] == 'REVIEW'],
            'compound_pharmacology_coverage': coverage, 'comparison_groups': groups,
            'docking_candidate_coverage': candidates, 'alpha1_alpha2_pair_candidates': pairs}
        for name, data in tables.items():
            write_csv(run_dir / 'tables' / (name + '.csv'), data, required=('activity_id',) if name.startswith(('bzd_', 'primary_')) else
                      ('pair_candidate_id', 'molecule_chembl_id') if name == 'alpha1_alpha2_pair_candidates' else
                      ('comparison_group_id',) if name == 'comparison_groups' else ('molecule_chembl_id',))
        (run_dir / 'processed/bzd_activity_classification.jsonl').write_text(''.join(encode(r) + '\n' for r in rows))
        if test_results:
            shutil.copyfile(test_results, run_dir / 'logs/pytest.xml')
        checks['input_snapshots_byte_identical'] = all(sha((run_dir / 'raw/lineage' / name).read_bytes()) == digest
                                                      for name, digest in input_hashes.items())
        checks['preexisting_tracked_files_unchanged'] = tracked_inventory() == before
        checks['jsonl_source_values_roundtrip'] = read_jsonl(run_dir / 'processed/bzd_activity_classification.jsonl') == rows
        with (run_dir / 'tables/bzd_activity_classification.csv').open(newline='') as handle:
            csv_rows = list(csv.DictReader(handle))
        checks['csv_id_and_relation_roundtrip'] = all(
            r['activity_id'] == str(a['activity_id']) and all(r[k] == ('' if a.get(k) is None else str(a[k]))
                for k in ['standard_value', 'standard_relation', 'value', 'relation', 'standard_units', 'molecule_chembl_id'])
            for r, a in zip(csv_rows, activities)) and len(csv_rows) == len(activities)
        if not all(checks.values()):
            raise ValueError('Artifact QC failed')
        qc = dict(status='PASS', checks=checks, tests=test_summary, source_activity_count=len(activities),
                  category_counts=dict(Counter(r['classification'] for r in rows)), source_raw_responses_verified=raw_count,
                  source_input_sha256=input_hashes, preexisting_tracked_file_count=len(before))
        write_json(run_dir / 'reports/QC_REPORT.json', qc)
        write_json(run_dir / 'reports/preexisting_files_sha256.json', before)
        (run_dir / 'reports/FINAL_REPORT.md').write_text(make_report(rows, groups, pairs, candidates, qc, config['source_run'], run_dir.name))
        manifest.update(completed_at=run_manager.iso_now(), qc_status='COMPLETE_WITH_MANUAL_REVIEW',
                        rule_version=config['rule_version'], category_counts=qc['category_counts'],
                        source_input_sha256=input_hashes,
                        output_files=sorted(str(p.relative_to(run_dir)) for p in run_dir.rglob('*') if p.is_file() and p != run_dir / 'run_manifest.json'))
        manifest['artifact_sha256'] = file_inventory(run_dir, [p for p in run_dir.rglob('*') if p.is_file() and p != run_dir / 'run_manifest.json'])
        run_manager.write_manifest(run_dir, manifest)  # Never call finalize_run: it changes LATEST_RUN.txt.
    except Exception as exc:
        manifest.update(completed_at=run_manager.iso_now(), qc_status='FAILED', error=str(exc))
        run_manager.write_manifest(run_dir, manifest)
        raise
    return run_dir, qc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'config/bzd_pharmacology_filter.json')
    parser.add_argument('--test-results', type=Path)
    args = parser.parse_args()
    run_dir, qc = run(args.config, args.test_results)
    print(encode({'run_directory': str(run_dir.relative_to(ROOT)), 'qc': qc['status'], 'category_counts': qc['category_counts']}))


if __name__ == '__main__':
    main()
