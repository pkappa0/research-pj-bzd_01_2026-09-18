"""Immutable two-drug descriptive integration; no fitted model or selectivity labels."""
from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal, InvalidOperation
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import run_manager
from .bzd_pharmacology_filter import (encode, sha, read_jsonl, write_json, write_csv,
                                      file_inventory, tests_summary, pairing_key)

ROOT = Path(__file__).resolve().parents[1]
DRUGS = ('diazepam', 'alprazolam')
SUBTYPES = ('alpha1', 'alpha2')
RECEPTORS = {'alpha1': 'alpha1_beta3_gamma2', 'alpha2': 'alpha2_beta3_gamma2_local_9CTJ'}
# Common-position labels use alpha1 numbering. Actual alpha2 residue numbers are retained separately.
SITES = {
    'gamma2_TYR58': ('BZD_GAMMA2_058', 'TYR', ('hydrophobic_interaction', 'pi_stack')),
    'gamma2_PHE77': ('BZD_GAMMA2_077', 'PHE', ('hydrophobic_interaction', 'pi_stack')),
    'alpha_HIS102': ('BZD_SITE_102', 'HIS', ('hydrophobic_interaction', 'pi_stack')),
    'alpha_LYS156': ('BZD_SITE_156', 'LYS', ('hydrogen_bond', 'hydrophobic_interaction')),
    'alpha_SER205': ('BZD_SITE_205', 'SER', ('hydrogen_bond',)),
    'gamma2_ASN60': ('BZD_GAMMA2_060', 'ASN', ('hydrogen_bond', 'hydrophobic_interaction', 'halogen_bond')),
}
TYPE_SHORT = {'hydrophobic_interaction': 'hydrophobic', 'pi_stack': 'pi', 'hydrogen_bond': 'H-bond', 'halogen_bond': 'halogen'}


def read_csv(path):
    import csv
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def finite(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None


def decimal_value(value):
    try:
        x = Decimal(str(value))
    except InvalidOperation:
        return None
    return x if x.is_finite() else None


def arithmetic(a, b):
    """No unit conversion or pooling. A percent-change ratio is deliberately withheld."""
    same = all(a.get(k) == b.get(k) for k in ['standard_type', 'standard_units', 'endpoint_class'])
    exact = a.get('standard_relation') == b.get('standard_relation') == '='
    x, y = decimal_value(a.get('standard_value')), decimal_value(b.get('standard_value'))
    comparable = same and exact and x is not None and y is not None and bool(a.get('standard_units'))
    difference = str(y - x) if comparable else None
    ratio = None
    reason = 'different_endpoint_unit_nonexact_or_missing_value'
    if comparable:
        if a['standard_type'] in {'Ki', 'Kd'} and a['endpoint_class'] == 'binding' and x > 0 and y > 0:
            ratio, reason = str(x / y), 'positive_exact_binding_constant_same_recorded_context'
        elif a['endpoint_class'] == 'functional_efficacy':
            reason = 'percent_change_normalization_and_dose_not_established'
        else:
            reason = 'not_an_established_positive_ratio_scale_for_this_analysis'
    return {'alpha2_minus_alpha1': difference,
            'difference_unit': ('percentage_points' if a.get('standard_units') == '%' else a.get('standard_units')) if comparable else None,
            'alpha1_over_alpha2_ratio': ratio, 'ratio_eligible': ratio is not None, 'ratio_reason': reason}


def extract_pairs(pair_records, activity_rows, compound_ids):
    activities = {r['activity_id']: r for r in activity_rows}
    if len(activities) != len(activity_rows):
        raise ValueError('Duplicate source activity ID')
    pairs, selected, seen = [], [], set()
    for p in pair_records:
        drug = p['drug_id']
        if drug not in compound_ids:
            continue
        if p['molecule_chembl_id'] != compound_ids[drug]:
            raise ValueError('Compound identity mismatch')
        # Do not choose a representative, average, or invent a Cartesian pair.
        ids = [json.loads(p[sub + '_activity_ids']) for sub in SUBTYPES]
        if any(len(v) != 1 for v in ids):
            raise ValueError('Pair has multiple observations; explicit observation pairing required')
        a, b = [activities[v[0]] for v in ids]
        unique_id = p['molecule_chembl_id'] + '__' + p['pair_candidate_id']
        if unique_id in seen:
            raise ValueError('Duplicate compound/pair context')
        seen.add(unique_id)
        key = json.loads(p['recorded_context_key'])
        for sub, row in zip(SUBTYPES, (a, b)):
            if (row['molecule_chembl_id'] != compound_ids[drug] or row['drug_id'] != drug or row['subtype'] != sub
                    or row['classification'] != 'PRIMARY' or row['review_flags'] or not row['point_comparison_usable']
                    or pairing_key(row) != key or row['standard_relation'] != '='):
                raise ValueError('Pair no longer satisfies recorded context or PRIMARY rules')
            if row['receptor_composition'] != p[sub + '_composition'] or [row['assay_chembl_id']] != json.loads(p[sub + '_assay_ids']):
                raise ValueError('Pair source identity/composition mismatch')
            selected.append(dict(row, integration_pair_id=unique_id))
        derived = arithmetic(a, b)
        row = dict(p, integration_pair_id=unique_id, **derived)
        row['recorded_context_key'] = key
        for sub, source in zip(SUBTYPES, (a, b)):
            for field in ['activity_id', 'assay_chembl_id', 'document_chembl_id', 'target_chembl_id', 'species',
                          'receptor_composition', 'endpoint_class', 'endpoint_detail', 'standard_type', 'standard_value',
                          'standard_units', 'standard_relation', 'value', 'relation', 'units', 'assay_type',
                          'assay_detail_description', 'assay_detail_assay_parameters', 'activity_properties',
                          'document_detail_title', 'document_detail_doi', 'document_detail_pubmed_id',
                          'target_organism', 'assay_detail_assay_organism', 'source_file', 'source_json_pointer', 'source_sha256']:
                row[sub + '_' + field] = source.get(field)
        row['value_interpretation'] = ('Smaller Ki means higher reported binding affinity; no categorical selectivity label'
                                       if a['standard_type'] == 'Ki' else
                                       'Reported electrophysiological percent change, not established Emax; no efficacy ratio')
        row['uncertainty'] = 'No matched replicate/SE information in selected records; nonzero difference is descriptive, not significance'
        row['structure_context_limit'] = ('beta2_pharmacology_vs_beta3_structures_and_provisional_alpha2_local_construct'
                                         if 'beta2' in a['receptor_composition'] else
                                         'beta3_family_matches_but_alpha2_is_local_construct_not_matched_full_pentamer')
        pairs.append(row)
    return pairs, selected


def pdb_residues(path):
    return {(line[21].strip(), int(line[22:26]), line[17:20].strip()) for line in path.read_text().splitlines()
            if line.startswith('ATOM  ')}


def fingerprint_rows(assigned, docking, mapping, prior_features, present_residues):
    """Extract requested features and recover proven absences from the complete contact archive."""
    output, raw_rows, geometry = [], [], []
    for sub in SUBTYPES:
        for drug in DRUGS:
            poses = [r for r in docking if r['drug_id'] == drug and r['receptor_id'] == RECEPTORS[sub] and r['status'] == 'docked']
            pose_ids = {str(r['pose_id']) for r in poses}
            if len(pose_ids) != len(poses):
                raise ValueError('Duplicate pose identifier')
            all_contacts = [r for r in assigned if r['drug_id'] == drug and r['receptor'] == sub]
            archived_poses = {str(r['pose_id']) for r in all_contacts}
            coverage_complete = bool(pose_ids) and archived_poses == pose_ids
            for site, (common, residue, kinds) in SITES.items():
                maps = [r for r in mapping if r['common_position'] == common]
                m = maps[0] if len(maps) == 1 else {}
                chain, number, actual = m.get(sub + '_chain'), finite(m.get(sub + '_residue_number')), m.get(sub + '_residue')
                mapping_valid = bool(m and str(m.get('qc_required')).lower() == 'false' and actual == residue
                                     and number is not None and (chain, int(number), actual) in present_residues[sub])
                for kind in kinds:
                    contacts = [r for r in all_contacts if r['common_position'] == common and r['interaction_type'] == kind]
                    for r in contacts:
                        if r['residue_chain'] != chain or finite(r['residue_number']) != number or r['raw_residue_name'] != actual:
                            raise ValueError('Contact mapping disagrees with source chain/residue')
                        if str(r['mapping_qc_required']).lower() != 'false':
                            raise ValueError('Unresolved source mapping')
                    observed = {str(r['pose_id']) for r in contacts}
                    known = coverage_complete and mapping_valid and observed <= pose_ids
                    n, count = len(pose_ids), len(observed)
                    prior = [r for r in prior_features if r['drug_id'] == drug and r['receptor'] == sub
                             and r['common_position'] == common and r['interaction_type'] == kind]
                    if len(prior) > 1:
                        raise ValueError('Duplicate prior feature')
                    freq = count / n if known else None
                    if prior and freq is not None and not math.isclose(float(prior[0]['frequency']), freq, abs_tol=1e-12):
                        raise ValueError('Recomputed frequency disagrees with existing corrected IFP')
                    if prior and known and (int(prior[0]['n_interacting_poses']) != count or int(prior[0]['interaction_count']) != len(contacts) or int(prior[0]['n_poses']) != n):
                        raise ValueError('Recomputed contact/pose counts disagree with existing corrected IFP')
                    attrs = [json.loads(r['raw_attributes']) for r in contacts]
                    types = {t: {str(r['pose_id']) for r, a in zip(contacts, attrs) if a.get('type') == t} for t in ('P', 'T')}
                    row = dict(drug_id=drug, receptor_subtype=sub, receptor_id=RECEPTORS[sub],
                               structure_pdb_id='6HUP' if sub == 'alpha1' else '9CTJ', structure_species='human',
                               structure_scope='full_pentamer_alpha1_beta3_gamma2L' if sub == 'alpha1' else 'provisional_local_CDE_beta3_alpha2_gamma2',
                               feature_id=site + '__' + kind, site_label=site, common_position=common,
                               residue_name=actual, residue_chain=chain, actual_residue_number=int(number) if number is not None else None,
                               interaction_type=kind, n_poses=n, n_interacting_poses=count if known else None,
                               interaction_count=len(contacts) if known else None, frequency=freq,
                               data_status=('observed' if count else 'verified_no_contact_in_archived_poses') if known else 'missing_or_unresolved_evidence',
                               feature_origin='existing_corrected_ifp' if prior else 'requested_feature_audited_against_complete_raw_contacts',
                               source_row_ids=[int(r['source_row_id']) for r in contacts], pose_ids=sorted(observed),
                               geometry_weighting='raw_contact_rows_not_independent_poses',
                               type_P_frequency=len(types['P']) / n if known and kind == 'pi_stack' else None,
                               type_T_frequency=len(types['T']) / n if known and kind == 'pi_stack' else None,
                               type_P_n_interacting_poses=len(types['P']) if known and kind == 'pi_stack' else None,
                               type_T_n_interacting_poses=len(types['T']) if known and kind == 'pi_stack' else None)
                    for measure, aliases in {'centdist': ['centdist'], 'angle': ['angle'], 'offset': ['offset'],
                                             'distance': ['dist_d-a', 'dist_d_a'] if kind == 'hydrogen_bond' else ['dist']}.items():
                        values = [finite(next((a[k] for k in aliases if k in a), None)) for a in attrs]
                        values = [x for x in values if x is not None]
                        row.update({measure + '_n_observed': len(values), measure + '_mean': statistics.mean(values) if values else None,
                                    measure + '_median': statistics.median(values) if values else None,
                                    measure + '_min': min(values) if values else None, measure + '_max': max(values) if values else None})
                    if prior and known:
                        for metric in ['centdist_mean', 'angle_mean', 'offset_mean', 'type_P_frequency', 'type_T_frequency']:
                            old, new = finite(prior[0][metric]), row[metric]
                            if (old is None) != (new is None) or (old is not None and not math.isclose(old, new, abs_tol=1e-12)):
                                raise ValueError('Geometry/type frequency disagrees with existing corrected IFP')
                    row['geometry_units'] = 'centdist/offset/distance: angstrom; angle: degrees'
                    output.append(row)
                    for raw, attr in zip(contacts, attrs):
                        raw_rows.append(dict(raw, integration_feature_id=row['feature_id']))
                        if kind == 'pi_stack':
                            geometry.append(dict(drug_id=drug, receptor_subtype=sub, site_label=site,
                                source_row_id=int(raw['source_row_id']), pose_id=raw['pose_id'], stacking_type=attr.get('type'),
                                centdist_angstrom=finite(attr.get('centdist')), angle_degrees=finite(attr.get('angle')),
                                offset_angstrom=finite(attr.get('offset')), raw_attributes=raw['raw_attributes']))
    return output, raw_rows, geometry


def structural_deltas(features):
    rows = []
    for drug in DRUGS:
        a = {r['feature_id']: r for r in features if r['drug_id'] == drug and r['receptor_subtype'] == 'alpha1'}
        b = {r['feature_id']: r for r in features if r['drug_id'] == drug and r['receptor_subtype'] == 'alpha2'}
        for fid, x in a.items():
            y = b[fid]
            row = dict(drug_id=drug, feature_id=fid, site_label=x['site_label'], interaction_type=x['interaction_type'],
                       alpha1_frequency=x['frequency'], alpha2_frequency=y['frequency'],
                       delta_alpha2_minus_alpha1=y['frequency'] - x['frequency'] if None not in (x['frequency'], y['frequency']) else None,
                       alpha1_n_interacting_poses=x['n_interacting_poses'], alpha2_n_interacting_poses=y['n_interacting_poses'],
                       alpha1_n_poses=x['n_poses'], alpha2_n_poses=y['n_poses'],
                       caveat='different receptor templates; alpha2 local construct; descriptive only')
            rows.append(row)
    return rows


def integrate(pairs, activities, features):
    by_pair = {p['integration_pair_id']: p for p in pairs}
    rows = []
    for a in activities:
        row = dict(a)
        p = by_pair[a['integration_pair_id']]
        row.update(matched_comparison_context=p['recorded_context_key'], pair_alpha2_minus_alpha1=p['alpha2_minus_alpha1'],
                   pair_alpha1_over_alpha2_ratio=p['alpha1_over_alpha2_ratio'], structure_context_limit=p['structure_context_limit'])
        for f in features:
            if f['drug_id'] != a['drug_id'] or f['receptor_subtype'] != a['subtype']:
                continue
            for k in ['frequency', 'n_interacting_poses', 'n_poses', 'data_status', 'actual_residue_number', 'residue_chain',
                      'centdist_mean', 'angle_mean', 'offset_mean', 'type_P_frequency', 'type_T_frequency']:
                row['ifp__' + f['feature_id'] + '__' + k] = f[k]
            row['ifp_receptor_id'] = f['receptor_id']
            row['ifp_structure_scope'] = f['structure_scope']
        rows.append(row)
    return rows


def heatmap_data(features):
    lookup = {(r['drug_id'], r['receptor_subtype'], r['feature_id']): r for r in features}
    labels, matrix = [], []
    for site, (_, _, kinds) in SITES.items():
        for kind in kinds:
            for field, suffix in [('frequency', TYPE_SHORT[kind])] + ([('type_P_frequency', 'pi / P'), ('type_T_frequency', 'pi / T')] if kind == 'pi_stack' else []):
                labels.append(site.replace('gamma2_', 'γ2 ').replace('alpha_', 'α ') + ' · ' + suffix)
                matrix.append([lookup[(drug, sub, site + '__' + kind)][field] for sub in SUBTYPES for drug in DRUGS])
    return labels, np.array([[np.nan if x is None else x for x in row] for row in matrix])


def draw_heatmap(ax, features, small=False):
    labels, data = heatmap_data(features)
    cmap = plt.get_cmap('YlGnBu').copy(); cmap.set_bad('#d9d9d9')
    im = ax.imshow(np.ma.masked_invalid(data), cmap=cmap, vmin=0, vmax=1, aspect='auto')
    ax.set_yticks(range(len(labels)), labels, fontsize=9 if small else 11)
    ax.set_xticks(range(4), ['Diazepam\nα1', 'Alprazolam\nα1', 'Diazepam\nα2*', 'Alprazolam\nα2*'], fontsize=10)
    ax.axvline(1.5, color='white', lw=4)
    for i in range(data.shape[0]):
        for j in range(4):
            v = data[i, j]
            ax.text(j, i, 'NA' if np.isnan(v) else f'{v:.2f}', ha='center', va='center', fontsize=9,
                    color='white' if not np.isnan(v) and v > .6 else '#17252b')
    return im


def plot_pair(ax, pair, idx):
    x, y = float(pair['alpha1_standard_value']), float(pair['alpha2_standard_value'])
    color = '#236aa4' if pair['drug_id'] == 'diazepam' else '#b05038'
    ax.plot([0, 1], [x, y], '-o', color=color, lw=2, markersize=7)
    ax.set_xticks([0, 1], ['α1', 'α2']); ax.set_xlim(-.4, 1.4)
    ax.set_ylim(min(0, x, y), max(x, y) * 1.42)
    ax.set_ylabel(f'{pair["standard_type"]} ({pair["standard_units"]})')
    beta = 'β2γ2' if 'beta2' in pair['alpha1_composition'] else 'β3γ2'
    ax.set_title(f'{idx}. {pair["drug_id"].capitalize()} · {beta}\n{pair["document_chembl_id"]}', fontsize=10, loc='left')
    for k, value in enumerate([x, y]):
        ax.annotate(f'{value:g}', (k, value), xytext=(0, 8), textcoords='offset points', ha='center', fontsize=11, color=color)
    diff = float(pair['alpha2_minus_alpha1'])
    ratio = f' | Ki α1/α2={float(pair["alpha1_over_alpha2_ratio"]):.3f}' if pair['ratio_eligible'] else ' | ratio withheld'
    difference_unit = ' pp' if pair['standard_units'] == '%' else ' ' + pair['standard_units']
    ax.text(.5, .96, f'Δ(α2−α1)={diff:+g}{difference_unit}{ratio}', transform=ax.transAxes, ha='center', va='top', fontsize=8)
    ax.grid(axis='y', alpha=.15); ax.spines[['top', 'right']].set_visible(False)


def figures(out, features, pairs, geometry):
    out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'axes.labelsize': 10, 'figure.facecolor': 'white'})
    paths = []
    def save(fig, name):
        for ext in ('png', 'pdf'):
            path = out / (name + '.' + ext)
            fig.savefig(path, dpi=170, bbox_inches='tight'); paths.append(path)
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 10), layout='constrained')
    im = draw_heatmap(ax, features); fig.colorbar(im, ax=ax, shrink=.75, label='Contacting poses / 9 archived poses')
    ax.set_title('Residue interaction fingerprint\nSeparate α1 and α2 blocks; interaction types retained', pad=18)
    fig.supxlabel('*α2: provisional 9CTJ local construct. Common labels use α1 numbering.\n0 = verified absence in archived poses; NA = missing evidence. P/T are PLIP stacking types.', fontsize=9)
    save(fig, 'residue_interaction_fingerprint_heatmap')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    for idx, (ax, p) in enumerate(zip(axes.flat, pairs), 1): plot_pair(ax, p, idx)
    fig.suptitle('Recorded-context α1/α2 pharmacology · human · relation =', fontsize=15)
    fig.supxlabel('Each panel is one compound × study × endpoint; no pooling. Smaller Ki = higher affinity.\nEfficacy is reported percent change, not established Emax. No uncertainty estimates available.', fontsize=10)
    save(fig, 'matched_pharmacology_comparison')
    fig = plt.figure(figsize=(15, 13), layout='constrained')
    gs = fig.add_gridspec(4, 2, width_ratios=[1.25, 1])
    left = fig.add_subplot(gs[:, 0]); im = draw_heatmap(left, features, small=True)
    left.set_title('Structural pose frequencies (0–1)', fontsize=13)
    fig.colorbar(im, ax=left, orientation='horizontal', shrink=.75, pad=.03, label='Contacting poses / 9')
    for idx, p in enumerate(pairs, 1): plot_pair(fig.add_subplot(gs[idx - 1, 1]), p, idx)
    fig.suptitle('Two-drug structure–pharmacology evidence map\nNo shared numerical scale, pooled estimate or fitted relationship', fontsize=16)
    fig.supxlabel('Structural α2* = provisional local construct; β2 pharmacology is a composition mismatch.\nLower Ki = higher affinity; response is not established Emax. Separate pharmacology units/scales.\nRecorded matches do not establish identical unreported conditions.', fontsize=10)
    save(fig, 'structure_pharmacology_summary')
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout='constrained')
    for ax, site in zip(axes, ['gamma2_TYR58', 'gamma2_PHE77', 'alpha_HIS102']):
        for drug in DRUGS:
            for sub in SUBTYPES:
                points = [r for r in geometry if r['drug_id'] == drug and r['receptor_subtype'] == sub and r['site_label'] == site]
                for typ in ('P', 'T'):
                    group = [r for r in points if r['stacking_type'] == typ and r['centdist_angstrom'] is not None and r['angle_degrees'] is not None]
                    if group:
                        ax.scatter([r['centdist_angstrom'] for r in group], [r['angle_degrees'] for r in group],
                                   color='#236aa4' if drug == 'diazepam' else '#b05038', marker='o' if typ == 'P' else '^',
                                   facecolors='none' if sub == 'alpha2' else None, s=60,
                                   label=f'{drug[:3]} {sub} {typ}')
        ax.set_title(site.replace('_', ' ')); ax.set_xlabel('Ring-centroid distance (Å)'); ax.set_ylabel('Ring angle (degrees)')
        ax.set_ylim(0, 95); ax.grid(alpha=.15); ax.legend(fontsize=7)
    fig.suptitle('Observed π-stacking geometry · individual raw contacts')
    fig.supxlabel('P: parallel; T: T-shaped. Hollow markers: α2 local construct. Missing contacts have no geometry.\nRepeated contacts within one pose are retained; these are not independent experimental replicates.', fontsize=9)
    save(fig, 'pi_stacking_geometry')
    return paths


def report_text(run_id, pairs, features, deltas, qc):
    lines = [f'# Diazepam / alprazolam integration — {run_id}', '',
        '## 目的と結論', '',
        '既存IFPとrecorded-context α1/α2薬理pairを並べた記述的PoC。新しいDocking、ML、相関検定、phenotype予測は実施していない。',
        '**残基別・interaction type別の構造差は存在するが、薬理差を一貫して説明することは未検証。単純な「contactが多い側ほどKiが低い」という一変数の説明は、このデータ全体には合わない。**', '',
        '## 1. Matched pharmacologyをpairごとに表示', '',
        'すべてhuman、relation `=`。Δはα2−α1、ratioはKi(α1)/Ki(α2)。元の数値文字列・ID・文献・組成・単位・relationを保持。study間の平均は計算しない。',
        'Kiは小さいほど報告されたbinding affinityが高い。ratio > 1はこのpairでα2のKiが低いことを示すが、selectivityカテゴリは付けない。', '',
        '| drug / document | endpoint・背景 | α1 value | α2 value | Δ α2−α1 | Ki α1/α2 | α1 activity / assay | α2 activity / assay |',
        '|---|---|---:|---:|---:|---:|---|---|']
    for p in pairs:
        ratio = f'{float(p["alpha1_over_alpha2_ratio"]):.6g}' if p['ratio_eligible'] else 'NA（未計算）'
        beta = 'β2γ2' if 'beta2' in p['alpha1_composition'] else 'β3γ2'
        lines.append(f'| {p["drug_id"]} / {p["document_chembl_id"]} | {p["standard_type"]} {p["standard_units"]} / {beta} | '
                     f'{p["alpha1_standard_value"]} | {p["alpha2_standard_value"]} | {p["alpha2_minus_alpha1"]} | {ratio} | '
                     f'{p["alpha1_activity_id"]} / {p["alpha1_assay_chembl_id"]} | {p["alpha2_activity_id"]} / {p["alpha2_assay_chembl_id"]} |')
    lines += ['', 'Efficacyの156%対89%は、記録された電気生理応答のpercent changeであり差は−67 percentage points。濃度・GABA条件・正規化基準の詳細が不十分なので、Emaxとは呼ばずratioも計算しない。',
        '各pairは前runの候補認定をそのまま使用。原著Methods/表を本runで再監査したものではなく、未記載条件・実験誤差・反復数の一致は未確認。数値差がゼロでないことと、実験的に有意な差があることを区別する。', '',
        '薬理文献（取得済みChEMBL document metadataから）：']
    seen = set()
    for p in pairs:
        doc = p['document_chembl_id']
        if doc in seen: continue
        seen.add(doc)
        doi = p['alpha1_document_detail_doi']
        lines.append(f'- {doc}: [{p["alpha1_document_detail_title"]}](https://doi.org/{doi})。PMID {p["alpha1_document_detail_pubmed_id"]}。')
    lines += ['', '## 2. 構造IFPと残基対応', '',
        '使用するのは `20260917_1706_structure_mouse_bridge` のreceptor-specific再対応済みIFP。以前のchain-D共有keyによる割当を再利用せず、元PLIP行と修正mappingを照合した。',
        'α1は6HUPのfull pentamer（human α1β3γ2L）。α2は9CTJのC/D/E鎖（β3/α2/γ2）だけの局所construct。',
        '9CTJ全体はβ2–α1–β3–α2–γ2のmixed native assemblyで、純粋なα2β3γ2 full pentamerではない。両側は厳密なmatched構造対ではない。',
        'β2薬理（両剤）とβ3構造の不一致、γ2 splice formの詳細不一致/未記載、6HUPのdiazepam-bound templateと9CTJ由来構造の状態差・box移送も交絡する。',
        'humanという種名が一致するだけでは実験条件やconstructの一致にはならない。', '',
        '| 共通ラベル | α1 実残基 | α2 実残基 |', '|---|---|---|',
        '| γ2 TYR58 | C:TYR58 | E:TYR58 |', '| γ2 PHE77 | C:PHE77 | E:PHE77 |',
        '| α HIS102 | D:HIS102 | D:HIS101 |', '| α LYS156 | D:LYS156 | D:LYS155 |',
        '| α SER205 | D:SER205 | D:SER204 |', '| γ2 ASN60 | C:ASN60 | E:ASN60 |', '',
        '36 archived poses（2剤×2系×9）を使用。頻度はunique interacting pose数 / 9。poseは同一seedの候補配座であり、9回の独立実験・平衡占有率・結合自由エネルギーではない。',
        '既存feature表に行がない要求featureも、mapping、prepared receptorの残基存在、全poseのcontact archiveを検証できた場合に限り「0 contacts」を計数。欠損値の0置換ではない。証拠不足ならNA。',
        'geometryは観測contact行のみの条件付き記述。π接触のないα2 TYR58/PHE77のcentdist/angle/offsetはNAのまま。', '',
        '| drug | site / interaction | α1 count/9 | α2 count/9 | Δ frequency α2−α1 |',
        '|---|---|---:|---:|---:|']
    for d in deltas:
        delta = 'NA' if d['delta_alpha2_minus_alpha1'] is None else f'{d["delta_alpha2_minus_alpha1"]:+.3f}'
        lines.append(f'| {d["drug_id"]} | {d["site_label"]} / {TYPE_SHORT[d["interaction_type"]]} | '
                     f'{d["alpha1_n_interacting_poses"]}/9 | {d["alpha2_n_interacting_poses"]}/9 | {delta} |')
    lines += ['', '### π-stacking geometry / P・T型', '',
        '[PLIPの定義](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md)でPはparallel、TはT-shaped。P/T frequencyもunique pose数 / 9。1 poseに複数contactがある場合、contact数とpose数は一致しない。',
        'centdist・offsetはÅ、angleはdegree。以下のmeanはraw contact行の記述統計で、精度や好ましさのscoreではない。P/T混合のangle平均だけで構造を代表させず、raw contactと型を併記する。', '',
        '| drug | block | site | raw contacts / poses | P / T frequency | centdist mean Å | angle mean ° | offset mean Å |',
        '|---|---|---|---:|---|---:|---:|---:|']
    fmt = lambda x: 'NA' if x is None else f'{x:.3f}'
    for f in features:
        if f['interaction_type'] != 'pi_stack': continue
        lines.append(f'| {f["drug_id"]} | {f["receptor_subtype"]} | {f["site_label"]} | {f["interaction_count"]} / {f["n_interacting_poses"]} | '
                     f'{fmt(f["type_P_frequency"])} / {fmt(f["type_T_frequency"])} | {fmt(f["centdist_mean"])} | {fmt(f["angle_mean"])} | {fmt(f["offset_mean"])} |')
    lines += ['', '## 3. Integrated comparison', '',
        '`integrated_structure_pharmacology.csv` は8 activity行。各行に元薬理fieldをすべて残し、対応するdrug×subtypeの6残基featureを横に付与する。',
        '`structure_pharmacology_differences.csv` は各pair×12featureの48行。薬理Δと構造Δを並べるための表で、構造featureをstudyごとに独立な観測として数えない。',
        'diazepamの同一IFPは3つの薬理pairに再利用される。この重複を標本数とする相関・回帰・有意差検定は行わない。', '',
        '### Supported observation', '',
        '- diazepamのKiはβ2 studyで14→20 nM（α2側が高い）、β3 studyで31→22 nM（α2側が低い）。study・背景によって方向が逆で、単一のsubtype affinity差へ統合できない。',
        '- alprazolamのβ2 Kiは0.8→0.6 nM。記録された数値差は−0.2 nM。誤差情報がないため統計的/生物学的有意性は判断不能。',
        '- 両剤でγ2 TYR58 πは4/9→0/9、PHE77 πはdiazepam 2/9→0/9、alprazolam 4/9→0/9。HIS102 πは1/9→2/9と1/9→3/9。',
        '- LYS156 hydrophobicはdiazepam 0/9→2/9、alprazolam 1/9→5/9。SER205 H-bondは2/9→1/9と2/9→0/9。ASN60はH-bond、hydrophobic、halogenを別々に保持。',
        '- α1 TYR58 π frequencyは両剤4/9で同じだがcentdist meanは3.9675と4.1775 Å。PHE77のhydrophobicは両剤7/9で同じだがπ頻度とP/T内訳は異なる。frequencyを1種類へ縮約すると情報が失われる。', '',
        '### Suggestive pattern', '',
        '- α1側のTYR58/PHE77 π接触の多さは、diazepam β2 pairの低いα1 Kiや大きいα1応答と定性的には並ぶ。しかし同じ構造パターンはdiazepam β3 Kiとalprazolam β2 Kiの方向には合わない。説明に都合のよいpairだけを選べない。',
        '- α2側のHIS102 πとLYS156 hydrophobicの増加は、diazepam β3およびalprazolam β2の低いα2 Kiと並ぶ可能性がある。しかしdiazepam β2では逆向きであり、異なる残基を事後的に選ぶことによる説明は未検証。',
        '- これらはmechanism候補の記述に留まる。全featureを保存する設計の情報保持には根拠があるが、geometryを足せば薬理の説明/予測が改善するという実証はない。', '',
        '### Not supported', '',
        '- 特定残基が薬理差を引き起こすこと、contact frequencyからaffinity・selectivityを予測できること。',
        '- study差をβ2/β3そのものの効果と断定すること（tracer、cell、構造、測定条件も異なる）。',
        '- 2剤だけから一般化すること、poseを独立標本として有意差・信頼区間を作ること。',
        '- IFPからphenotype・鎮静・副作用を説明すること、Kiとfunctional responseを同じラベルへ統合すること。', '',
        '## 4. 今回の6つの判断', '',
        '1. **薬理の数値差は支持**。4 pairすべてに非ゼロ差。ただし有意性・再現性・生物学的重要性は判断不能。diazepam Kiの方向はstudy依存。',
        '2. **定性的に並ぶfeatureはあるが、全pairに一貫した対応は支持されない**。TYR58/PHE77とHIS102/LYS156で逆向きの候補を作れてしまい、独立検証なしには説明力を選べない。',
        '3. **frequencyだけで十分とは支持できない**。同じfrequencyでもgeometry/P・Tが異なる事実がある。frequency-onlyモデルの性能を評価したわけではない。',
        '4. **geometryとinteraction typeを保持することは情報損失回避の点で支持**。薬理説明への増分寄与は判断不能。geometryがないcontactはNAを保つ。',
        '5. **検証目的の限定的拡張は条件付きで合理的**。原データを保持した設計で反証テストを増やせる。一方、予測力や創薬上の価値は現在の2剤から判断不能。追加薬の大量Dockingよりmatched構造と薬理条件の確保を優先する。',
        '6. **最優先の反証可能な仮説 H1**：両剤で観測したγ2 TYR58 πの「α1 > α2」が局所construct/単一seedだけの産物ではなく、揃えたhuman α1β3γ2 / α2β3γ2 full-pentamer条件でも再現する。',
        '   次段階では同一調製・box規則・pose採択規則・事前指定した独立5 seedsで比較し、各剤でΔ(α2−α1)<0が4/5 seeds以上かつseed中央値<0を暫定的な再現基準として事前登録する。満たさなければこのfeatureの頑健性仮説を支持せず、construct/pose依存性を再検討する。これは生物学的因果の検定ではなく構造特徴の再現性の検証であり、今回追加計算は実施していない。',
        '   H1が再現した場合も、同一study/β背景で両剤のKiと誤差を取得して薬理対応を別に検証する。「γ2接触が多いほど常にKiが低い」は本PoC全体では支持されないため、予測ルールとして採用しない。', '',
        '## 5. 図', '',
        '![IFP heatmap](figures/residue_interaction_fingerprint_heatmap.png)',
        '![Matched pharmacology](figures/matched_pharmacology_comparison.png)',
        '![Integrated summary](figures/structure_pharmacology_summary.png)',
        '![Geometry](figures/pi_stacking_geometry.png)', '',
        '各図はPNG/PDFを保存。薬理とfrequencyは別axis・別unit。欠損geometryはplotしない。薬理に根拠のないerror barを付けない。', '',
        '## 6. QC・provenance', '',
        f'QC **{qc["status"]}**、{len(qc["checks"])}項目PASS。tests: `{encode(qc["tests"])}`。',
        '4 pair / 8 activity / 48 IFP feature / 24 subtype feature差。元activityを平均・削除・単位変換していない。',
        '入力のSHA-256、元activity JSONL、PLIP raw属性、mapping、pose表、PDB/prepared receptor、解析コード・config・testsのsnapshotを保存。',
        'IFPの既存featureのfrequencyおよびcontact/pose数を再計数で照合。全poseのcontact archiveの存在、PDB/prepared receptorの残基存在を検証するが、PLIP再実行による再検出は本runでは行っていない。',
        'コードとmanifestは実行前Git HEADに加えてsnapshot hashで識別。`LATEST_RUN.txt`と既存tracked filesの実行前後hashは不変。',
        '再実行: `MPLBACKEND=Agg .venv/bin/python -m pytest -q --junitxml=/tmp/integration-tests.xml`、続いて',
        '` .venv/bin/python -m src.diazepam_alprazolam_integration --test-results /tmp/integration-tests.xml`。毎回新runとなる。', '',
        '構造注釈は同梱PDB header・既存調製コードと[RCSB 6HUP](https://www.rcsb.org/structure/6HUP)、[RCSB 9CTJ](https://www.rcsb.org/structure/9CTJ)（2026-09-18閲覧）を照合。原著薬理Methodsの追加確認は未完了として扱う。', '']
    return '\n'.join(lines)


def prepare(config):
    pharm = ROOT / config['pharmacology_run']; struct = ROOT / config['structure_run']
    manifest = json.loads((pharm / 'run_manifest.json').read_text())
    for name, digest in manifest['artifact_sha256'].items():
        if sha((pharm / name).read_bytes()) != digest:
            raise ValueError('Pharmacology source artifact hash mismatch')
    all_activity = read_jsonl(pharm / 'processed/bzd_activity_classification.jsonl')
    pair_records = read_csv(pharm / 'tables/alpha1_alpha2_pair_candidates.csv')
    pairs, activities = extract_pairs(pair_records, all_activity, config['compound_ids'])
    assigned = read_csv(struct / 'data/processed/reassigned_interactions.csv')
    docking = read_csv(struct / 'data/raw/lineage/docking_results_4drug.csv')
    mapping = read_csv(struct / 'results/tables/corrected_residue_mapping.csv')
    prior = read_csv(struct / 'results/tables/candidate_residue_features.csv')
    present = {sub: pdb_residues(ROOT / 'data/raw/structures/strict3' / (pdb + '.pdb')) &
                   pdb_residues(ROOT / 'data/raw/strict3/receptors' / (RECEPTORS[sub] + '.pdbqt'))
               for sub, pdb in [('alpha1', '6HUP'), ('alpha2', '9CTJ')]}
    features, raw_contacts, geometry = fingerprint_rows(assigned, docking, mapping, prior, present)
    for feature in features:
        feature['molecule_chembl_id'] = config['compound_ids'][feature['drug_id']]
    deltas = structural_deltas(features)
    integrated = integrate(pairs, activities, features)
    differences = [dict(d, integration_pair_id=p['integration_pair_id'], document_chembl_id=p['document_chembl_id'],
                        pharmacology_endpoint=p['standard_type'], pharmacology_unit=p['standard_units'],
                        pharmacology_alpha1=p['alpha1_standard_value'], pharmacology_alpha2=p['alpha2_standard_value'],
                        pharmacology_delta_alpha2_minus_alpha1=p['alpha2_minus_alpha1'],
                        pharmacology_ratio_alpha1_over_alpha2=p['alpha1_over_alpha2_ratio'],
                        structure_context_limit=p['structure_context_limit']) for p in pairs for d in deltas if p['drug_id'] == d['drug_id']]
    checks = {
        'exactly_four_source_pair_records_eight_distinct_activities': len(pairs) == 4 and len(activities) == 8 and len({a['activity_id'] for a in activities}) == 8,
        'two_id_verified_compounds': {a['molecule_chembl_id'] for a in activities} == set(config['compound_ids'].values()),
        'all_activity_fields_unchanged': all(all(a[k] == v for k, v in next(r for r in all_activity if r['activity_id'] == a['activity_id']).items()) for a in activities),
        'no_pooled_or_cartesian_pairs': len(pairs) == len([p for p in pair_records if p['drug_id'] in DRUGS]),
        'three_positive_exact_binding_ratios_one_withheld_efficacy_ratio': sum(p['ratio_eligible'] for p in pairs) == 3 and all(p['standard_type'] == 'Ki' for p in pairs if p['ratio_eligible']),
        'all_numeric_differences_exact_decimal': all(p['alpha2_minus_alpha1'] == str(Decimal(p['alpha2_standard_value']) - Decimal(p['alpha1_standard_value'])) for p in pairs),
        'requested_48_features_unique': len(features) == len({(r['drug_id'], r['receptor_subtype'], r['feature_id']) for r in features}) == 48,
        'all_36_pose_records_available_with_contact_archives_and_mapped_residues': all(r['n_poses'] == 9 and r['frequency'] is not None for r in features),
        'frequency_unique_poses_over_denominator': all(r['frequency'] == r['n_interacting_poses'] / r['n_poses'] and 0 <= r['frequency'] <= 1 for r in features),
        'no_geometry_zero_imputation_for_absent_contacts': all(r['centdist_mean'] is None and r['angle_mean'] is None and r['offset_mean'] is None for r in features if r['interaction_count'] == 0),
        'all_raw_geometry_records_traceable': len({r['source_row_id'] for r in geometry}) == len(geometry),
        'integrated_eight_activity_rows_no_source_field_loss': len(integrated) == 8 and all(all(row[k] == v for k, v in a.items()) for row, a in zip(integrated, activities)),
        'subtype_residue_numbers_not_conflated': all(r['actual_residue_number'] == (102 if r['receptor_subtype'] == 'alpha1' else 101) for r in features if r['site_label'] == 'alpha_HIS102'),
        'background_mismatch_and_local_construct_flags_present': all('local_construct' in r['structure_context_limit'] for r in integrated),
    }
    if not all(checks.values()):
        raise ValueError('QC failed: ' + ', '.join(k for k, v in checks.items() if not v))
    return pairs, activities, features, raw_contacts, geometry, deltas, integrated, differences, checks


def tracked_inventory():
    names = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().split('\0')
    return file_inventory(ROOT, [ROOT / name for name in names if name])


def run(config_path, test_results=None):
    config = json.loads(config_path.read_text()); before = tracked_inventory()
    pairs, activities, features, raw_contacts, geometry, deltas, integrated, differences, checks = prepare(config)
    p, s = config['pharmacology_run'], config['structure_run']
    source_names = [p + '/run_manifest.json', p + '/tables/alpha1_alpha2_pair_candidates.csv',
        p + '/processed/bzd_activity_classification.jsonl', s + '/run_manifest.json',
        s + '/data/processed/reassigned_interactions.csv', s + '/data/raw/lineage/docking_results_4drug.csv',
        s + '/results/tables/corrected_residue_mapping.csv', s + '/results/tables/candidate_residue_features.csv',
        'data/raw/structures/strict3/6HUP.pdb', 'data/raw/structures/strict3/9CTJ.pdb',
        'data/raw/strict3/receptors/alpha1_beta3_gamma2.pdbqt', 'data/raw/strict3/receptors/alpha2_beta3_gamma2_local_9CTJ.pdbqt',
        'src/strict_poc.py', 'src/structure_mouse_bridge.py']
    code_names = ['src/diazepam_alprazolam_integration.py', 'src/bzd_pharmacology_filter.py', 'src/acquire_chembl.py',
                  'src/run_manager.py', 'tests/test_diazepam_alprazolam_integration.py', 'pytest.ini']
    source_hashes = file_inventory(ROOT, [ROOT / name for name in source_names])
    test_info = tests_summary(test_results)
    run_dir, manifest = run_manager.create_run(short_phase_name='diazepam_alprazolam_integration',
        analysis_phase='Two-drug matched pharmacology and residue IFP descriptive integration',
        analysis_purpose='Expose concordant and discordant patterns without fitted model or causal inference',
        parent_run=Path(p).name, drugs=list(DRUGS), receptors=RECEPTORS, docking_parameters={'new_docking': False},
        input_files=source_names, notes=['No pooling, prediction or phenotype analysis', 'Alpha2 provisional local construct', 'LATEST_RUN.txt unchanged'],
        source_run=p, source_files=source_names, code_paths=[ROOT / name for name in code_names])
    try:
        for name in source_names:
            dest = run_dir / 'raw/lineage' / name; dest.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(ROOT / name, dest)
        for name in code_names:
            dest = run_dir / 'config/code_snapshot' / name; dest.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(ROOT / name, dest)
        shutil.copyfile(config_path, run_dir / 'config/diazepam_alprazolam_integration.json')
        tables = {'matched_pharmacology_pairs': pairs, 'diazepam_alprazolam_ifp': features,
                  'integrated_structure_pharmacology': integrated, 'structural_subtype_differences': deltas,
                  'structure_pharmacology_differences': differences, 'pi_stacking_geometry_contacts': geometry,
                  'selected_raw_plip_contacts': raw_contacts}
        for name, records in tables.items(): write_csv(run_dir / (name + '.csv'), records)
        for name, records in [('selected_pharmacology', activities), ('integrated_structure_pharmacology', integrated), ('diazepam_alprazolam_ifp', features)]:
            (run_dir / 'processed' / (name + '.jsonl')).write_text(''.join(encode(r) + '\n' for r in records))
        if test_results: shutil.copyfile(test_results, run_dir / 'logs/pytest.xml')
        figure_paths = figures(run_dir / 'figures', features, pairs, geometry)
        checks['eight_figure_artifacts_nonempty'] = len(figure_paths) == 8 and all(x.stat().st_size > 1000 for x in figure_paths)
        checks['input_snapshots_byte_identical'] = all(sha((run_dir / 'raw/lineage' / n).read_bytes()) == v for n, v in source_hashes.items())
        checks['all_preexisting_tracked_files_unchanged'] = before == tracked_inventory()
        checks['jsonl_null_and_numeric_string_roundtrip'] = read_jsonl(run_dir / 'processed/integrated_structure_pharmacology.jsonl') == integrated
        roundtrip = read_csv(run_dir / 'matched_pharmacology_pairs.csv')
        checks['csv_original_values_relations_and_ids_preserved'] = all(all(r[sub + '_' + k] == str(p[sub + '_' + k]) for sub in SUBTYPES
            for k in ['activity_id', 'assay_chembl_id', 'standard_value', 'standard_relation', 'standard_units']) for r, p in zip(roundtrip, pairs)) and len(roundtrip) == len(pairs)
        if not all(checks.values()): raise ValueError('Artifact QC failed')
        qc = dict(status='PASS', checks=checks, tests=test_info, pair_count=len(pairs), activity_count=len(activities),
                  feature_count=len(features), geometry_contact_count=len(geometry), source_input_sha256=source_hashes,
                  preexisting_tracked_file_count=len(before), limitations=['No experimental uncertainty estimates', 'Alpha2 local construct',
                    'Beta2 pharmacology versus beta3 structural background', 'Raw PLIP contacts audited but detection not rerun'])
        write_json(run_dir / 'QC_REPORT.json', qc)
        write_json(run_dir / 'reports/preexisting_files_sha256.json', before)
        (run_dir / 'FINAL_REPORT.md').write_text(report_text(run_dir.name, pairs, features, deltas, qc))
        manifest.update(completed_at=run_manager.iso_now(), qc_status='PASS_WITH_INTERPRETATION_LIMITS', source_input_sha256=source_hashes,
                        artifact_sha256=file_inventory(run_dir, [x for x in run_dir.rglob('*') if x.is_file() and x != run_dir / 'run_manifest.json']))
        manifest['output_files'] = list(manifest['artifact_sha256'])
        run_manager.write_manifest(run_dir, manifest)
    except Exception as exc:
        manifest.update(qc_status='FAILED', error=str(exc), completed_at=run_manager.iso_now()); run_manager.write_manifest(run_dir, manifest); raise
    return run_dir, qc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'config/diazepam_alprazolam_integration.json')
    parser.add_argument('--test-results', type=Path)
    args = parser.parse_args()
    out, qc = run(args.config, args.test_results)
    print(encode({'run_directory': str(out.relative_to(ROOT)), 'qc': qc['status'], 'pairs': qc['pair_count'], 'features': qc['feature_count']}))


if __name__ == '__main__': main()
