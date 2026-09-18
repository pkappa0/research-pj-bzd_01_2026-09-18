import numpy as np
import pandas as pd
from itertools import product

def similarity_matrix(binary, freq):
    drugs = list(binary.get('drug_id', []))
    if not drugs: return pd.DataFrame(columns=['drug_id_a','drug_id_b','jaccard_similarity','frequency_euclidean_distance','cluster_order'])
    b = binary.set_index('drug_id').reindex(drugs).fillna(0); f = freq.set_index('drug_id').reindex(drugs).fillna(0) if len(freq.columns)>1 else pd.DataFrame(0.0, index=drugs, columns=[])
    rows=[]
    for a, bid in product(drugs, drugs):
        av, bv = b.loc[a].to_numpy(), b.loc[bid].to_numpy(); union = np.sum((av > 0) | (bv > 0)); inter = np.sum((av > 0) & (bv > 0))
        rows.append({'drug_id_a':a,'drug_id_b':bid,'jaccard_similarity':float(inter / union) if union else np.nan,'frequency_euclidean_distance':float(np.linalg.norm(f.loc[a].to_numpy(dtype=float)-f.loc[bid].to_numpy(dtype=float)))})
    out = pd.DataFrame(rows)
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
        mat=np.zeros((len(drugs),len(drugs)))
        for i in range(len(drugs)):
            for j in range(len(drugs)):
                u=np.sum((b.iloc[i]>0)|(b.iloc[j]>0)); inter=np.sum((b.iloc[i]>0)&(b.iloc[j]>0)); mat[i,j]=1-inter/u if u else 0
        order=leaves_list(linkage(squareform(mat,checks=False),method='average')) if len(drugs)>1 else np.array([0]); order_map={drugs[int(i)]:int(k) for k,i in enumerate(order)}
    except Exception: order_map={d:i for i,d in enumerate(drugs)}
    out['cluster_order']=out.drug_id_a.map(order_map); return out
