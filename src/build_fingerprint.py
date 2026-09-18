import pandas as pd

def build_fingerprints(df, mapping):
    qc = {'n_interaction_rows': int(len(df)), 'n_drugs': int(df.drug_id.nunique()) if not df.empty else 0,
          'n_receptors': int(df.receptor.nunique()) if not df.empty else 0,
          'n_poses': int(df[['drug_id','receptor','pose_id']].drop_duplicates().shape[0]) if not df.empty else 0,
          'unmapped_interaction_rows': 0, 'dropped_zero_variance_features': 0, 'pose_counts_by_drug': ''}
    if df.empty: return pd.DataFrame(columns=['drug_id']), pd.DataFrame(columns=['drug_id']), qc
    x = df.merge(mapping[['receptor','original_residue','common_position']], how='left', on=['receptor','original_residue'])
    x['common_position'] = x['common_position'].astype('string'); missing = x['common_position'].isna() | (x['common_position'].str.strip() == '')
    qc['unmapped_interaction_rows'] = int(missing.sum())
    x.loc[missing, 'common_position'] = 'UNMAPPED_' + x.loc[missing, 'receptor'].astype(str) + '_' + x.loc[missing, 'original_residue'].astype(str)
    x['feature'] = x['common_position'].astype(str) + '_' + x['interaction_type'].str.upper().str.replace(r'[^A-Z0-9]+', '_', regex=True)
    poses = x[['drug_id','receptor','pose_id']].drop_duplicates(); qc['pose_counts_by_drug'] = ';'.join(f'{k}:{v}' for k,v in poses.groupby('drug_id').size().items())
    present = x[['drug_id','receptor','pose_id','feature']].drop_duplicates()
    b = pd.crosstab(present.drug_id, present.feature).clip(upper=1).reset_index(); counts = present.groupby(['drug_id','feature']).size().unstack(fill_value=0)
    f = counts.div(poses.groupby('drug_id').size(), axis=0).reset_index()
    cols = [c for c in b.columns if c != 'drug_id' and b[c].nunique(dropna=False) > 1]
    qc['dropped_zero_variance_features'] = int(len(b.columns) - 1 - len(cols)); b = b[['drug_id'] + cols]; f = f[['drug_id'] + [c for c in cols if c in f.columns]]
    return b, f, qc
