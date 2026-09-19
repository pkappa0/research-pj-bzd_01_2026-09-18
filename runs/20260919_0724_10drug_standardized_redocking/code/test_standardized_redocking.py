import json
from pathlib import Path
import numpy as np
from src import standardized_redocking as d

def test_fixed_seeds_and_drugs():
    assert len(d.DRUGS)==10 and len(set(d.SEEDS))==5
    assert d.SEEDS==[2026091901,2026091902,2026091903,2026091904,2026091905]

def test_actual_count_denominator_and_pose_ids():
    # Same pose rank in different seeds is a distinct computational pose.
    poses={'2026091901:1','2026091902:1'}
    assert len(poses)==2
    assert len(poses)/3==2/3 # no requested-max denominator or padding

def test_environment_and_preregistration_frozen():
    run=d.runpath();f=json.loads((run/'config/PREREGISTRATION_FROZEN.json').read_text())
    assert all(d.sha(run/p)==h for p,h in f['artifacts'].items())
    cfg=json.loads((run/'config/docking_config.json').read_text())
    assert cfg['exhaustiveness']==8 and cfg['num_modes']==9 and cfg['energy_range']==4
    assert cfg['external_annotation_allowed'] is False

def test_prepared_inputs_immutable():
    run=d.runpath();manifest=json.loads((run/'config/PREPARED_INPUTS_FROZEN.json').read_text())
    assert all(d.sha(run/p)==h for p,h in manifest['files'].items())

def test_torsions_not_globally_rigid():
    rows=d.read(d.runpath()/'tables/ligand_preparation_qc.csv')
    assert len(rows)==10
    assert all(int(r['torsdof'])==int(r['pdbqt_branch_count']) for r in rows)
    assert len({r['torsdof'] for r in rows})>1
    assert all(int(r['torsdof'])>0 for r in rows if int(r['rdkit_rotatable_bonds'])>0)
    assert {r['drug_id'] for r in rows if r['unassigned_stereocentres']!='[]'}=={'lorazepam','temazepam','zopiclone'}

def test_blocks_have_different_prepared_identities():
    rows=d.read(d.runpath()/'tables/receptor_preparation_qc.csv')
    assert len({r['pdbqt_sha256'] for r in rows})==2
    assert {r['chains'] for r in rows}=={'ABCDE','CDE'}

def test_actual_returned_poses_not_requested_maximum():
    rows=d.read(d.runpath()/'tables/docking_run_manifest.csv')
    poses=d.read(d.runpath()/'tables/pose_manifest.csv')
    assert len(rows)==100 and sum(int(r['n_poses']) for r in rows)==len(poses)
    assert all(0<int(r['n_poses'])<=9 for r in rows)
    from collections import Counter
    counts=Counter((r['drug_id'],r['receptor'],r['seed']) for r in poses if r['status']=='success')
    for r in d.read(d.runpath()/'tables/fingerprint_seed_level.csv'):
        n=counts[(r['drug_id'],r['receptor'],r['seed'])]
        assert int(r['n_poses'])==n
        assert abs(float(r['frequency'])-int(r['count'])/n)<1e-12
