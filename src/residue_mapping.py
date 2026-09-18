import pandas as pd
from .parse_interactions import resolve_input

COLS = ['receptor', 'original_residue', 'common_position', 'mapping_confidence', 'mapping_note']
def load_mapping(path):
    p = resolve_input(path, 'residue_mapping.csv')
    if p is None: return pd.DataFrame(columns=COLS), {'path': str(path), 'status': 'missing', 'missing_columns': ''}
    df = pd.read_csv(p); missing = [c for c in COLS[:3] if c not in df.columns]
    if missing: return pd.DataFrame(columns=COLS), {'path': str(p), 'status': 'invalid', 'missing_columns': '|'.join(missing)}
    for c in COLS[3:]:
        if c not in df: df[c] = pd.NA
    return df[COLS].drop_duplicates(), {'path': str(p), 'status': 'loaded', 'missing_columns': ''}
