"""Load PLIP/Vina-derived interaction records without inventing missing data."""
from pathlib import Path
import pandas as pd

REQUIRED = ['drug_id', 'receptor', 'pose_id', 'original_residue', 'interaction_type']
OPTIONAL = ['drug_name_raw', 'pdb_id', 'chain', 'residue_name', 'residue_number', 'source_file', 'score', 'confidence']

def resolve_input(path, expected_name=None):
    p = Path(path)
    if p.is_file(): return p
    if p.is_dir():
        candidate = p / expected_name if expected_name else None
        if candidate and candidate.exists(): return candidate
        csvs = sorted(c for c in p.glob('*.csv') if not c.name.endswith('.template.csv'))
        return csvs[0] if len(csvs) == 1 else None
    return None

def load_interactions(path):
    p = resolve_input(path, 'interactions.csv')
    if p is None: return pd.DataFrame(columns=REQUIRED), {'path': str(path), 'status': 'missing', 'missing_columns': ''}
    df = pd.read_csv(p)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing: return pd.DataFrame(columns=REQUIRED), {'path': str(p), 'status': 'invalid', 'missing_columns': '|'.join(missing)}
    df = df.copy(); df['drug_id'] = df['drug_id'].astype('string').str.strip().str.lower()
    for c in ['receptor', 'pose_id', 'original_residue', 'interaction_type']: df[c] = df[c].astype('string').str.strip()
    return df.drop_duplicates(), {'path': str(p), 'status': 'loaded', 'missing_columns': ''}
