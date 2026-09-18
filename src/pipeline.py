from pathlib import Path
import pandas as pd
from .parse_interactions import load_interactions, resolve_input
from .residue_mapping import load_mapping
from .build_fingerprint import build_fingerprints
from .similarity import similarity_matrix
from .process_pharmacology import process_pharmacology
from .process_phenotype import process_phenotype
from .visualization import make_figures
from .exploratory_analysis import correlations
def run(cfg):
    out={k:Path(v) for k,v in cfg['outputs'].items()}; [p.mkdir(parents=True,exist_ok=True) for p in out.values()]
    interactions,iqc=load_interactions(cfg['inputs']['interactions']); mapping,mqc=load_mapping(cfg['inputs']['residue_mapping']); binary,freq,fqc=build_fingerprints(interactions,mapping); pharm=process_pharmacology(cfg['inputs']['pharmacology']); phen=process_phenotype(cfg['inputs']['phenotype'])
    for name,df in [('fingerprint_binary',binary),('fingerprint_frequency',freq),('residue_mapping',mapping),('pharmacology_processed',pharm),('phenotype_processed',phen)]: df.to_csv(out['processed']/f'{name}.csv',index=False)
    sim=similarity_matrix(binary,freq); sim.to_csv(out['processed']/'fingerprint_similarity.csv',index=False); drugs=set(binary.get('drug_id',[]))|set(pharm.get('drug_id',[]))|set(phen.get('drug_id',[])); integrated=pd.DataFrame({'drug_id':sorted(drugs)})
    for d in [pharm,phen]:
        if not d.empty: integrated=integrated.merge(d,on='drug_id',how='left')
    order=list(sim.drop_duplicates('drug_id_a').sort_values('cluster_order').drug_id_a) if not sim.empty else sorted(drugs); integrated['cluster_order']=integrated.drug_id.map({d:i for i,d in enumerate(order)}); integrated=integrated.sort_values(['cluster_order','drug_id']).drop(columns='cluster_order'); integrated.to_csv(out['processed']/'integrated_dataset.csv',index=False)
    make_figures(binary,freq,sim,pharm,phen,out['figures'],order); correlations(binary,pharm,phen).to_csv(out['tables']/'exploratory_correlations.csv',index=False)
    rows=[{'input':'interactions',**iqc},{'input':'residue_mapping',**mqc},{'input':'fingerprint','status':'loaded' if not interactions.empty else 'no_records','path':'derived','missing_columns':'','records':len(interactions),'unmapped_rows':fqc['unmapped_interaction_rows'],'n_drugs':fqc['n_drugs'],'n_poses':fqc['n_poses'],'pose_counts_by_drug':fqc['pose_counts_by_drug']}]
    for label,d,path,name in [('pharmacology',pharm,cfg['inputs']['pharmacology'],'pharmacology.csv'),('phenotype',phen,cfg['inputs']['phenotype'],'phenotype.csv')]:
        resolved = resolve_input(path, name)
        status = 'loaded' if not d.empty else ('empty' if resolved is not None else 'missing')
        rows.append({'input':label,'status':status,'path':str(path),'missing_columns':'','records':len(d)})
    pd.DataFrame(rows).to_csv(out['tables']/'qc_summary.csv',index=False); return {'order':order,'qc':rows}
