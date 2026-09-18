from src import ten_drug_preflight as p

def test_environment_gate_fails_closed():
    assert not p.gate({},True)
    evidence={k:'documented' for k in ['vina_version','rdkit_version','plip_version','plip_effective_config','openbabel_version']}
    assert not p.gate(evidence,False)
    assert p.gate(evidence,True)
    evidence['plip_effective_config']=None
    assert not p.gate(evidence,True)

def test_legacy_four_drug_frequencies():
    rows=p.check_legacy()
    assert len(rows)==100
    assert all(r['match'] and r['n_poses']==9 for r in rows)
    assert len({r['drug_id'] for r in rows})==4

def test_preregistration_is_structural_only():
    import json
    c=json.loads((p.ROOT/'config/ten_drug_unsupervised_preregistered.json').read_text())
    assert len(c['drugs'])==10 and len(c['new_drugs'])==6
    assert c['primary']['linkage']=='average' and c['primary']['distance']=='euclidean'
    assert c['docking']['versions'] is None
    assert c['supervised_training'] is False
    assert 'before any activity/phenotype annotation load' in c['external_annotation_gate']
