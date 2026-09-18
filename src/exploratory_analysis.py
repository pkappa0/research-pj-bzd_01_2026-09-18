import numpy as np
import pandas as pd
def _numeric(df): return df.set_index('drug_id').select_dtypes(include='number') if not df.empty else pd.DataFrame()
def correlations(binary, pharm, phen):
    cols=['comparison','method','n_pairs','coefficient','p_value','note']
    if len(binary)<3: return pd.DataFrame([{'comparison':'all','method':'none','n_pairs':0,'coefficient':np.nan,'p_value':np.nan,'note':'At least 3 drugs are required; no inferential statistic was computed.'}],columns=cols)
    ids=list(binary.drug_id); x=binary.set_index('drug_id').reindex(ids).fillna(0).to_numpy(); fs=[]
    for i in range(len(ids)):
        for j in range(i+1,len(ids)):
            u=((x[i]>0)|(x[j]>0)).sum(); fs.append(((x[i]>0)&(x[j]>0)).sum()/u if u else np.nan)
    rows=[]
    for label,other in [('pharmacology',pharm),('phenotype',phen)]:
        n=_numeric(other).reindex(ids)
        if n.empty: rows.append({'comparison':f'fingerprint_similarity_vs_{label}','method':'none','n_pairs':0,'coefficient':np.nan,'p_value':np.nan,'note':'No numeric input columns.'}); continue
        vals=[]
        for i in range(len(ids)):
            for j in range(i+1,len(ids)):
                if pd.notna(n.iloc[i]).any() and pd.notna(n.iloc[j]).any():
                    vals.append(np.linalg.norm(n.iloc[i].fillna(n.iloc[i].mean()).to_numpy()-n.iloc[j].fillna(n.iloc[j].mean()).to_numpy()))
        if len(vals)<3: rows.append({'comparison':f'fingerprint_similarity_vs_{label}','method':'none','n_pairs':len(vals),'coefficient':np.nan,'p_value':np.nan,'note':'Fewer than 3 complete drug pairs.'}); continue
        from scipy.stats import spearmanr
        m=min(len(fs),len(vals)); r,p=spearmanr(fs[:m],vals[:m]); rows.append({'comparison':f'fingerprint_similarity_vs_{label}','method':'Spearman','n_pairs':m,'coefficient':r,'p_value':p,'note':'Exploratory pairwise association; no causal interpretation.'})
    return pd.DataFrame(rows,columns=cols)
