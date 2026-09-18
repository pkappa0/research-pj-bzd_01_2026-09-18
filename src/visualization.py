from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
def _heat(df,path,title,order=None):
    if df.empty or len(df.columns)<=1:
        fig,ax=plt.subplots(figsize=(6,2)); ax.text(.5,.5,'No input data available',ha='center'); ax.axis('off')
    else:
        x=df.set_index('drug_id'); x=x.reindex([d for d in order if d in x.index]) if order else x; x=x.select_dtypes(include='number')
        if x.empty: return _heat(pd.DataFrame(),path,title,order)
        fig,ax=plt.subplots(figsize=(max(6,x.shape[1]*.25),max(2,x.shape[0]*.3))); im=ax.imshow(x.to_numpy(dtype=float),aspect='auto',cmap='viridis'); ax.set_yticks(range(len(x.index)),x.index); ax.set_xticks(range(len(x.columns)),x.columns,rotation=90,fontsize=7); ax.set_title(title); fig.colorbar(im,ax=ax)
    fig.tight_layout(); fig.savefig(path,dpi=160); plt.close(fig)
def make_figures(binary,freq,sim,pharm,phen,out,order=None):
    out=Path(out); out.mkdir(parents=True,exist_ok=True); _heat(binary,out/'fingerprint_binary_heatmap.png','Binary interaction fingerprint',order); _heat(freq,out/'fingerprint_frequency_heatmap.png','Interaction frequency',order)
    if not sim.empty:
        similarity_plot = sim.pivot(index='drug_id_a',columns='drug_id_b',values='jaccard_similarity').reset_index().rename(columns={'drug_id_a':'drug_id'})
        _heat(similarity_plot,out/'fingerprint_similarity_heatmap.png','Fingerprint similarity',order)
    else: _heat(sim,out/'fingerprint_similarity_heatmap.png','Fingerprint similarity',order)
    _heat(pharm,out/'pharmacology_heatmap.png','Pharmacology',order); _heat(phen,out/'phenotype_heatmap.png','Phenotype',order)
    pieces=[d.set_index('drug_id').select_dtypes(include='number') for d in [binary,pharm,phen] if not d.empty]; integrated=pd.concat(pieces,axis=1).reset_index() if pieces else pd.DataFrame(); _heat(integrated,out/'integrated_heatmap.png','Integrated view (fixed fingerprint order)',order)
