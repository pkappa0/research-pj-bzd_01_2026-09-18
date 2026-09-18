from pathlib import Path
import yaml
from src.pipeline import run

def test_synthetic_pipeline(tmp_path):
    root=Path(__file__).parent/'fixtures'; out=tmp_path
    cfg={'inputs': {'interactions': str(root/'interactions.csv'), 'pharmacology': str(root/'pharmacology.csv'), 'phenotype': str(root/'phenotype.csv'), 'residue_mapping': str(root/'residue_mapping.csv')}, 'outputs': {'processed': str(out/'processed'), 'figures': str(out/'figures'), 'tables': str(out/'tables')}}
    result=run(cfg)
    assert result['order'] and len(result['order'])==3
    assert (out/'processed/fingerprint_binary.csv').exists()
    assert (out/'figures/integrated_heatmap.png').exists()
    assert (out/'tables/qc_summary.csv').exists()

def test_missing_inputs_do_not_stop(tmp_path):
    cfg={'inputs': {'interactions': str(tmp_path/'missing-interactions.csv'), 'pharmacology': str(tmp_path/'missing-pharmacology.csv'), 'phenotype': str(tmp_path/'missing-phenotype.csv'), 'residue_mapping': str(tmp_path/'missing-mapping.csv')}, 'outputs': {'processed': str(tmp_path/'processed'), 'figures': str(tmp_path/'figures'), 'tables': str(tmp_path/'tables')}}
    run(cfg)
    text=(tmp_path/'tables/qc_summary.csv').read_text()
    assert 'missing' in text
