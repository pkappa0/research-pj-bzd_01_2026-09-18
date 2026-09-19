"""Layout-only derivatives of frozen tables: no refitting or regrouping."""
from .standardized_redocking import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

def main():
    run=runpath();f=json.loads((run/'STRUCTURAL_RESULTS_FROZEN.json').read_text());assert all(sha(run/p)==h for p,h in f['artifacts'].items())
    out=run/'figures/reviewed';out.mkdir(exist_ok=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
    def savefig(fig,name):
        for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=150,bbox_inches='tight')
        plt.close(fig)
    C={b:read(run/'tables'/f'fingerprint_type_frequency_{b}.csv') for b in BLOCKS}
    for rep,name in [('B_alpha1','alpha1_10drug_heatmap'),('B_alpha2','alpha2_10drug_heatmap'),('E','multireceptor_10drug_heatmap')]:
        rows=C[rep[2:]] if rep!='E' else [{**a,**b} for a,b in zip(C['alpha1'],C['alpha2'])]
        cols=[k for k in rows[0] if k!='drug_id'];X=np.array([[float(r[c]) for c in cols] for r in rows]);Z=np.array([[float(r[k]) for k in ['left','right','distance','n_leaves']] for r in read(run/'tables'/f'linkage_{rep}.csv')])
        fig=plt.figure(figsize=(14,max(9,len(cols)/5)),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,9],width_ratios=[40,1]);top=fig.add_subplot(gs[0,0]);order=dendrogram(Z,ax=top,no_labels=True)['leaves'];top.set_title(rep+' · Euclidean / average linkage · frozen drug order');top.set_xticks([])
        ax=fig.add_subplot(gs[1,0]);im=ax.imshow(X[order].T,vmin=0,vmax=1,aspect='auto',cmap='YlGnBu');ax.set_xticks(range(10),np.array(DRUGS)[order],rotation=45,ha='right');ax.set_yticks(range(len(cols)),[c.replace('alpha1|','α1 ').replace('alpha2|','α2 ').replace('BZD_GAMMA2_','γ2 ').replace('BZD_SITE_','α ').replace('hydrophobic_interaction','hyd').replace('hydrogen_bond','H-bond').replace('|frequency','') for c in cols],fontsize=7)
        fig.colorbar(im,cax=fig.add_subplot(gs[1,1]),label='Unique contact poses / actual evaluable poses')
        fig.supxlabel('Type-resolved rows; dendrogram from residue-union B/E. α2: provisional local construct. Layout revision only.',fontsize=9);savefig(fig,name)
    scores=read(run/'tables/pca_scores.csv');clusters=read(run/'tables/cluster_assignments.csv');fig,axes=plt.subplots(1,3,figsize=(20,7),layout='constrained')
    for ax,rep in zip(axes,['B_alpha1','B_alpha2','E']):
        rows=[r for r in scores if r['representation']==rep];labels={r['drug_id']:int(r['cluster']) for r in clusters if r['representation']==rep}
        x=np.array([float(r['PC1']) for r in rows]);y=np.array([float(r['PC2']) for r in rows]);ax.scatter(x,y,c=[labels[r['drug_id']] for r in rows],cmap='tab10',s=65)
        for i,r in enumerate(rows):
            # Label offsets affect typography only; points and loadings untouched.
            dy={'triazolam':-15,'temazepam':12,'clonazepam':2}.get(r['drug_id'],6)
            dx=5
            if r['drug_id']=='triazolam':dx,dy=(-65,8) if rep=='B_alpha2' else (-65,-3) if rep=='B_alpha1' else (-20,-18)
            ax.annotate(r['drug_id'],(x[i],y[i]),xytext=(dx,dy),textcoords='offset points',fontsize=9,arrowprops=dict(arrowstyle='-',color='#888',lw=.5))
        ax.margins(.30);ax.set(title=rep,xlabel=f"PC1 ({float(rows[0]['PC1_explained_variance']):.1%})",ylabel=f"PC2 ({float(rows[0]['PC2_explained_variance']):.1%})");ax.grid(alpha=.15)
    fig.suptitle('Frozen structural-only PCA · k=3 structural colors · layout-only derivative');savefig(fig,'fingerprint_pca')
    save(run/'reports/PRESENTATION_DERIVATIVES.json',{'created_at':datetime.now().astimezone().isoformat(),'reason':'Align dendrogram and heatmap panel widths; separate nearby PCA labels. No numeric data, cluster, metric, threshold or PCA change. Original frozen figures retained.','frozen_source_sha256':sha(run/'STRUCTURAL_RESULTS_FROZEN.json'),'artifacts':{str(p.relative_to(run)):sha(p) for p in out.iterdir()}})
    assert all(sha(run/p)==h for p,h in f['artifacts'].items())
if __name__=='__main__':main()
