"""Structural-only matrices, clustering, PCA and computational seed sensitivity."""
from .standardized_redocking import *
import numpy as np
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage,fcluster,dendrogram
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def mapping_index(run):
    source=ROOT/'runs/20260917_1706_structure_mouse_bridge/results/tables/corrected_residue_mapping.csv'
    rows=read(source);idx={};qc=[]
    for b in BLOCKS:
        pdb=(run/'raw/receptors'/f'{b}_prepared_H.pdb').read_text();res={(l[21],int(l[22:26])):l[17:20].strip() for l in pdb.splitlines() if l.startswith('ATOM')}
        for r in rows:
            if r['qc_required'].lower()=='true' or not r['mapping_status'].startswith('mapped'):continue
            try:key=(r[b+'_chain'],int(float(r[b+'_residue_number'])))
            except (ValueError,KeyError):continue
            expected=r[b+'_residue'];present=res.get(key)==expected
            qc.append(dict(receptor=b,common_position=r['common_position'],chain=key[0],residue_number=key[1],expected_residue=expected,observed_residue=res.get(key),mapping_valid=present))
            if present:idx[(b,*key)]=(r['common_position'],expected)
    table(run/'tables/residue_mapping_qc.csv',qc)
    for b in BLOCKS:
        for cp in ['BZD_GAMMA2_058','BZD_GAMMA2_060','BZD_GAMMA2_077','BZD_SITE_102','BZD_SITE_156','BZD_SITE_205']:
            assert any(r['receptor']==b and r['common_position']==cp and r['mapping_valid'] for r in qc),cp
    return idx

def matrices():
    run=runpath();idx=mapping_index(run);contacts=json.loads((WORK/'all_contacts.json').read_text());poses=json.loads((WORK/'all_poses.json').read_text())
    evaluable={(d,b,s):{r['pose_id'] for r in poses if r['drug_id']==d and r['receptor']==b and r['seed']==s and r['status']=='success'} for d in DRUGS for b in BLOCKS for s in SEEDS}
    groups={};unmapped=[]
    for c in contacts:
        key=(c['receptor'],c['residue_chain'],int(c['residue_number']))
        if key not in idx:unmapped.append(c);continue
        cp,name=idx[key];feature=(c['receptor'],cp,name,c['interaction_type']);groups.setdefault(feature,[]).append(c)
    table(run/'tables/unmapped_contacts.csv',unmapped or [dict(status='none')])
    schema=sorted(groups);seedrows=[];raw=[];featurecols=[]
    for b,cp,name,typ in schema:
        base=f'{b}|{cp}|{name}|{typ}';cs=groups[(b,cp,name,typ)]
        metrics=['frequency']+(['P_frequency','T_frequency','centdist_mean','angle_mean','offset_mean'] if typ=='pi_stack' else [])
        for metric in metrics:featurecols.append(dict(feature=base+'|'+metric,receptor=b,common_position=cp,residue_name=name,interaction_type=typ,metric=metric,unit='degrees' if metric=='angle_mean' else 'angstrom' if metric.endswith('_mean') else 'fraction'))
        for d in DRUGS:
            drugcontacts=[c for c in cs if c['drug_id']==d];totalvalid=sum(len(evaluable[(d,b,s)]) for s in SEEDS)
            for metric in metrics:
                key=base+'|'+metric
                if metric.endswith('frequency'):
                    subset=drugcontacts if metric=='frequency' else [c for c in drugcontacts if c['raw_attributes'].get('type')==metric[0]]
                    freqs=[]
                    for s in SEEDS:
                        valid=evaluable[(d,b,s)];count=len({c['pose_id'] for c in subset if c['seed']==s} & valid);value=count/len(valid) if valid else None
                        seedrows.append(dict(drug_id=d,receptor=b,feature=key,seed=s,count=count,n_poses=len(valid),frequency=value));freqs.append(value)
                    count=len({c['pose_id'] for c in subset});value=count/totalvalid if totalvalid else None;observed=[f for f in freqs if f is not None]
                    raw.append(dict(drug_id=d,receptor=b,feature=key,feature_type='frequency',total_count=count,total_evaluable_poses=totalvalid,value=value,geometry_n_observed=None,seed_range=max(observed)-min(observed) if observed else None,seed_sd=float(np.std(observed,ddof=1)) if len(observed)>1 else None,unstable_feature=bool(observed and max(observed)-min(observed)>.4)))
                else:
                    values=[float(c['raw_attributes'][metric[:-5]]) for c in drugcontacts if c['raw_attributes'].get(metric[:-5]) is not None]
                    raw.append(dict(drug_id=d,receptor=b,feature=key,feature_type='geometry',total_count=None,total_evaluable_poses=totalvalid,value=float(np.mean(values)) if values else None,geometry_n_observed=len(values),seed_range=None,seed_sd=None,unstable_feature=None))
    # B collapses types through pose union, never summing marginal frequencies.
    B={};Bseed={};Bcols={};unionrows=[]
    for b in BLOCKS:
        sites=sorted({(cp,name) for bb,cp,name,t in schema if bb==b});cols=[f'{b}|{cp}|{name}|any_contact_frequency' for cp,name in sites];Bcols[b]=cols
        B[b]=np.zeros((10,len(sites)));Bseed[b]={s:np.zeros_like(B[b]) for s in SEEDS}
        for j,(cp,name) in enumerate(sites):
            cs=[c for (bb,cc,nn,t),rr in groups.items() if (bb,cc,nn)==(b,cp,name) for c in rr]
            for i,d in enumerate(DRUGS):
                count=0;den=0;freqs=[]
                for s in SEEDS:
                    valid=evaluable[(d,b,s)];poseids={c['pose_id'] for c in cs if c['drug_id']==d and c['seed']==s} & valid;n=len(poseids);v=n/len(valid) if valid else np.nan
                    count+=n;den+=len(valid);Bseed[b][s][i,j]=v;freqs.append(v)
                    seedrows.append(dict(drug_id=d,receptor=b,feature=cols[j],seed=s,count=n,n_poses=len(valid),frequency=None if np.isnan(v) else v))
                B[b][i,j]=count/den if den else np.nan
                unionrows.append(dict(drug_id=d,receptor=b,feature=cols[j],total_count=count,total_evaluable_poses=den,aggregate_frequency=None if not den else count/den,seed_range=float(np.nanmax(freqs)-np.nanmin(freqs)),seed_sd=float(np.nanstd(freqs,ddof=1)),unstable_feature=bool(np.nanmax(freqs)-np.nanmin(freqs)>.4)))
    def wide(arr,cols):return [dict(drug_id=d,**{c:None if np.isnan(v) else float(v) for c,v in zip(cols,arr[i])}) for i,d in enumerate(DRUGS)]
    E=np.concatenate([B[b] for b in BLOCKS],axis=1);ecols=sum([Bcols[b] for b in BLOCKS],[])
    for b in BLOCKS:table(run/'tables'/f'fingerprint_10drug_frequency_{b}.csv',wide(B[b],Bcols[b]))
    table(run/'tables/fingerprint_10drug_multireceptor.csv',wide(E,ecols));table(run/'tables/fingerprint_10drug_binary.csv',wide(np.where(np.isnan(E),np.nan,(E>0).astype(float)),ecols))
    table(run/'tables/fingerprint_seed_level.csv',seedrows);table(run/'tables/fingerprint_aggregate_long.csv',raw)
    table(run/'tables/feature_seed_robustness.csv',[dict(r,representation='C') for r in raw if r['feature_type']=='frequency']+[dict(r,representation='B') for r in unionrows])
    table(run/'tables/feature_schema.csv',featurecols)
    rawindex={(r['drug_id'],r['feature']):r['value'] for r in raw};allcols=[c['feature'] for c in featurecols]
    allraw=np.array([[np.nan if rawindex[(d,c)] is None else rawindex[(d,c)] for c in allcols] for d in DRUGS]);table(run/'tables/fingerprint_10drug_raw.csv',wide(allraw,allcols));table(run/'tables/fingerprint_missingness.csv',wide(np.isnan(allraw).astype(float),allcols))
    C={b:[c['feature'] for c in featurecols if c['receptor']==b and c['metric'].endswith('frequency')] for b in BLOCKS}
    Cm={b:np.array([[rawindex[(d,c)] for c in C[b]] for d in DRUGS]) for b in BLOCKS}
    for b in BLOCKS:table(run/'tables'/f'fingerprint_type_frequency_{b}.csv',wide(Cm[b],C[b]))
    return B,Bseed,Bcols,E,ecols,C,Cm,allraw,featurecols,raw

def structural_analysis():
    run=runpath();B,Bseed,Bcols,E,ecols,C,Cm,allraw,schema,raw=matrices()
    reps={**{f'B_{b}':(B[b],Bcols[b]) for b in BLOCKS},'E':(E,ecols)}
    outputs={};scores=[];loadings=[];clusters=[];distrows=[];cosrows=[];jacrows=[];nrows=[];stability=[];supp=[];exclusions=[];candidate=[]
    for rep,(X,cols) in reps.items():
        complete=np.isfinite(X).all(axis=0);x=X[:,complete];usecols=np.array(cols)[complete];assert len(usecols)>1
        for c,ok in zip(cols,complete):exclusions.append(dict(representation=rep,feature=c,included_in_clustering_pca=bool(ok),reason='all_drugs_observed' if ok else 'missing_in_one_or_more_drugs'))
        dist=squareform(pdist(x));cos=1-squareform(pdist(x,'cosine'));binary=x>0;jac=1-squareform(pdist(binary,'jaccard'))
        for i,d in enumerate(DRUGS):
            distrows.append(dict(representation=rep,drug_id=d,**dict(zip(DRUGS,dist[i]))));cosrows.append(dict(representation=rep,drug_id=d,**dict(zip(DRUGS,cos[i]))));jacrows.append(dict(representation=rep,drug_id=d,**dict(zip(DRUGS,jac[i]))))
            for j,e in enumerate(DRUGS):nrows.append(dict(representation=rep,drug_i=d,drug_j=e,comparable_features=int(np.sum(np.isfinite(X[i])&np.isfinite(X[j]))),analysis_common_features=len(usecols)))
        Z=linkage(pdist(x),method='average');labels=fcluster(Z,3,criterion='maxclust');variable=np.ptp(x,axis=0)>0;xc=x[:,variable]-x[:,variable].mean(axis=0);U,S,V=np.linalg.svd(xc,full_matrices=False)
        for k in range(len(S)):
            j=np.argmax(abs(V[k]));sign=1 if V[k,j]>=0 else -1;V[k]*=sign;U[:,k]*=sign
        pc=U[:,:2]*S[:2];variance=S*S/np.sum(S*S)
        for i,d in enumerate(DRUGS):
            scores.append(dict(representation=rep,drug_id=d,PC1=pc[i,0],PC2=pc[i,1],PC1_explained_variance=variance[0],PC2_explained_variance=variance[1]));clusters.append(dict(representation=rep,drug_id=d,cluster=int(labels[i]),k_requested=3,k_actual=len(set(labels)),method='euclidean_average'))
        for j,c in enumerate(usecols[variable]):
            loadings.append(dict(representation=rep,feature=c,PC1_loading=V[0,j],PC2_loading=V[1,j],drug_range=float(np.ptp(x[:,variable][:,j])),small_range_le_2_over_9=bool(np.ptp(x[:,variable][:,j])<=2/9)))
        table(run/'tables'/f'linkage_{rep}.csv',[dict(left=int(a),right=int(b),distance=c,n_leaves=int(d)) for a,b,c,d in Z])
        outputs[rep]=(x,usecols,Z,labels,pc,variance,dist,cos)
        for metric in ['euclidean','cosine']:
            z=linkage(pdist(x,metric),method='average')
            for k in [2,3,4]:
                lab=fcluster(z,k,criterion='maxclust');supp.append(dict(representation=rep,metric=metric,k=k,aggregate_primary_ARI=adjusted_rand_score(labels,lab),assignments=dict(zip(DRUGS,map(int,lab)))))
        for s in SEEDS:
            sx=Bseed[rep[2:]][s] if rep.startswith('B_') else np.concatenate([Bseed[b][s] for b in BLOCKS],axis=1)
            sx=sx[:,complete];sd=pdist(sx);sl=fcluster(linkage(sd,method='average'),3,criterion='maxclust')
            stability.append(dict(representation=rep,seed=s,aggregate_distance_spearman=float(spearmanr(pdist(x),sd).statistic),aggregate_cluster_ARI=float(adjusted_rand_score(labels,sl)),seed_assignments=dict(zip(DRUGS,map(int,sl))),interpretation='computational_search_sensitivity_not_biological_replicates'))
        for j,c in enumerate(usecols):
            means={int(k):float(x[labels==k,j].mean()) for k in set(labels)};difference=max(means.values())-min(means.values())
            candidate.append(dict(representation=rep,feature=c,cluster_mean_range=difference,drug_range=float(np.ptp(x[:,j])),small_drug_range=bool(np.ptp(x[:,j])<=2/9),cluster_means=means))
    # All prespecified supplementary representations, never chosen by annotation.
    for b in BLOCKS:
        for rep,x in [('A',B[b]>0),('C',Cm[b]),('D',allraw[:,[s['receptor']==b for s in schema]])]:
            valid=np.isfinite(x).all(axis=0)
            if rep=='D':valid &= np.ptp(x,axis=0)>0
            xx=x[:,valid]
            if rep=='D':xx=(xx-xx.mean(axis=0))/xx.std(axis=0)
            metric='jaccard' if rep=='A' else 'euclidean';labs=fcluster(linkage(pdist(xx,metric),method='average'),3,criterion='maxclust')
            supp.append(dict(representation=rep+'_'+b,metric=metric,k=3,aggregate_primary_ARI=adjusted_rand_score(outputs['B_'+b][3],labs),assignments=dict(zip(DRUGS,map(int,labs))),complete_nonconstant_features=int(valid.sum()),excluded_missing_or_constant=int((~valid).sum())))
    for filename,rows in [('similarity_euclidean',distrows),('similarity_cosine',cosrows),('similarity_jaccard',jacrows),('similarity_comparable_feature_counts',nrows),('pca_scores',scores),('pca_loadings',loadings),('cluster_assignments',clusters),('seed_cluster_stability',stability),('supplementary_clustering',supp),('analysis_feature_inclusion',exclusions),('cluster_feature_contrasts',candidate)]:table(run/'tables'/(filename+'.csv'),rows)
    table(run/'tables/pca_loading_candidates.csv',sorted(loadings,key=lambda r:(r['representation'],-max(abs(r['PC1_loading']),abs(r['PC2_loading'])),r['feature'])))
    plot_structures(run,outputs,C,Cm,stability)
    save(run/'reports/STRUCTURAL_SUMMARY.json',dict(representations={k:dict(features=len(v[1]),PC1=float(v[5][0]),PC2=float(v[5][1]),clusters={d:int(c) for d,c in zip(DRUGS,v[3])}) for k,v in outputs.items()},seed_stability=stability,unmapped_contacts=len(read(run/'tables/unmapped_contacts.csv')),frequency_features=sum(r['metric'].endswith('frequency') for r in schema),geometry_features=sum(r['metric'].endswith('mean') for r in schema)))
    save(run/'STRUCTURAL_RESULTS_FROZEN.json',dict(frozen_at=datetime.now().astimezone().isoformat(),external_annotation_loaded=False,config_sha256=sha(run/'config/unsupervised_analysis_config.json'),artifacts={str(p.relative_to(run)):sha(p) for directory in ['tables','figures'] for p in sorted((run/directory).iterdir()) if p.is_file()},code_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'src/standardized_redocking.py',ROOT/'src/standardized_plip.py']}))
    print('STRUCTURE FROZEN',run,flush=True)

def plot_structures(run,outputs,C,Cm,stability):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
    def savefig(fig,name):
        for ext in ['png','pdf']:fig.savefig(run/'figures'/(name+'.'+ext),dpi=150,bbox_inches='tight')
        plt.close(fig)
    for rep,(x,cols,z,lab,pc,var,dist,cos) in outputs.items():
        fig=plt.figure(figsize=(14, max(7,len(C['alpha1'])/5) if rep=='B_alpha1' else max(7,len(C['alpha2'])/5) if rep=='B_alpha2' else max(10,(len(C['alpha1'])+len(C['alpha2']))/6)))
        gs=fig.add_gridspec(2,1,height_ratios=[1,7]);ax=fig.add_subplot(gs[0]);order=dendrogram(z,labels=DRUGS,ax=ax)['leaves'];ax.set_title(rep+' · fixed structural block(s) · Euclidean / average linkage');ax.tick_params(axis='x',bottom=False,labelbottom=False)
        mat=Cm[rep[2:]] if rep.startswith('B_') else np.concatenate(list(Cm.values()),axis=1);names=C[rep[2:]] if rep.startswith('B_') else sum(list(C.values()),[])
        ax=fig.add_subplot(gs[1]);im=ax.imshow(mat[order].T,aspect='auto',vmin=0,vmax=1,cmap='YlGnBu');ax.set_xticks(range(10),np.array(DRUGS)[order],rotation=45,ha='right');ax.set_yticks(range(len(names)),[n.replace('BZD_GAMMA2_','γ2 ').replace('BZD_SITE_','α ').replace('|frequency','').replace('hydrophobic_interaction','hyd').replace('hydrogen_bond','H-bond') for n in names],fontsize=7)
        fig.colorbar(im,ax=ax,label='Unique contact poses / evaluable poses');fig.text(.1,.01,'Rows: type-resolved C; drug order: primary residue-union B/E. α2 provisional local construct.',fontsize=9)
        savefig(fig,{'B_alpha1':'alpha1_10drug_heatmap','B_alpha2':'alpha2_10drug_heatmap','E':'multireceptor_10drug_heatmap'}[rep])
    fig,axes=plt.subplots(1,3,figsize=(19,6),layout='constrained')
    for ax,(rep,v) in zip(axes,outputs.items()):
        im=ax.imshow(v[7],vmin=0,vmax=1,cmap='viridis');ax.set(title=rep+' cosine similarity',xticks=range(10),yticks=range(10),xticklabels=DRUGS,yticklabels=DRUGS);ax.tick_params(axis='x',rotation=90);fig.colorbar(im,ax=ax,shrink=.6)
    savefig(fig,'structural_similarity_matrix')
    fig,axes=plt.subplots(1,3,figsize=(19,6),layout='constrained')
    for ax,(rep,v) in zip(axes,outputs.items()):
        pc=v[4];var=v[5];ax.scatter(pc[:,0],pc[:,1],c=v[3],cmap='tab10',s=65)
        for i,d in enumerate(DRUGS):ax.annotate(d,pc[i],xytext=(4,5),textcoords='offset points',fontsize=8)
        ax.margins(.25);ax.set(title=rep,xlabel=f'PC1 ({var[0]:.1%})',ylabel=f'PC2 ({var[1]:.1%})');ax.axhline(0,color='#ddd',lw=.5);ax.axvline(0,color='#ddd',lw=.5)
    fig.suptitle('Structural-only PCA · centered frequencies · color = fixed k=3 structural cluster');savefig(fig,'fingerprint_pca')
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for rep in outputs:
        rows=[r for r in stability if r['representation']==rep]
        for ax,key in zip(axes,['aggregate_distance_spearman','aggregate_cluster_ARI']):ax.plot(range(1,6),[r[key] for r in rows],'-o',label=rep);ax.set(xlabel='Docking seed index',ylabel=key,ylim=(-1.05,1.05),xticks=range(1,6));ax.legend()
    fig.suptitle('Seed sensitivity, not biological confidence intervals');savefig(fig,'seed_robustness_summary')

if __name__=='__main__':structural_analysis()
