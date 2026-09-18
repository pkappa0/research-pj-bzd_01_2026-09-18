import pandas as pd
from .parse_interactions import resolve_input
def process_pharmacology(path):
    p=resolve_input(path,'pharmacology.csv')
    if p is None: return pd.DataFrame(columns=['drug_id'])
    d=pd.read_csv(p)
    if 'drug_id' not in d: return pd.DataFrame(columns=['drug_id'])
    d=d.copy(); d['drug_id']=d.drug_id.astype('string').str.strip().str.lower(); return d
