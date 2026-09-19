"""Ten-drug standardized computation. No pharmacology imports in this module."""
from pathlib import Path
import csv,json,hashlib,subprocess,sys,os,re,gzip,shutil,platform,importlib.metadata as md
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT.parent/'standardized_redocking'
ENV=ROOT.parent/'standardized_env/venv'
BASE=Path('/opt/homebrew/Caskroom/miniforge/base/envs/plip_env')
VINA=Path('/opt/homebrew/bin/vina').resolve()
DRUGS=['diazepam','alprazolam','triazolam','zolpidem','lorazepam','clonazepam','midazolam','temazepam','zopiclone','zaleplon']
SEEDS=[2026091901,2026091902,2026091903,2026091904,2026091905]
BLOCKS={'alpha1':('6HUP','ABCDE'),'alpha2':('9CTJ','CDE')}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,default=str)+'\n')
def read(p):return list(csv.DictReader(Path(p).open()))
def table(p,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,keys);w.writeheader()
        for r in rows:w.writerow({k:'NA' if r.get(k) is None else json.dumps(r[k],sort_keys=True) if isinstance(r[k],(list,dict)) else r[k] for k in keys})
def runpath():return Path((WORK/'active_run.txt').read_text().strip())
def configure_plip():
    from plip.basic import config
    config.NOHYDRO=True;config.NOFIX=True;config.NOFIXFILE=True;config.BREAKCOMPOSITE=True;config.XML=True;config.TXT=True
    config.PYMOL=False;config.PICS=False;config.MAXTHREADS=1
    return config

def initialize():
    run=ROOT/'runs'/(datetime.now().astimezone().strftime('%Y%m%d_%H%M')+'_10drug_standardized_redocking');run.mkdir(exist_ok=False)
    for d in ['tables','figures','reports','config','raw/receptors','raw/ligands','raw/archives','logs','code']:(run/d).mkdir(parents=True)
    (WORK/'active_run.txt').write_text(str(run));(WORK/'jobs').mkdir(exist_ok=True)
    packages={d.metadata['Name']:d.version for d in md.distributions()}
    from openbabel import openbabel as ob
    versions={'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'architecture':platform.machine(),'packages':packages,'vina_version':subprocess.check_output([str(VINA),'--version'],text=True).strip(),'vina_binary_sha256':sha(VINA),'vina_binary':str(VINA),'openbabel_library_version':ob.OBReleaseVersion(),'openbabel_conda_package_version':'3.1.1','plip_version':'3.0.0','recorded_before_docking':datetime.now().astimezone().isoformat()}
    save(run/'config/software_versions.json',versions)
    freeze=subprocess.check_output([sys.executable,'-m','pip','freeze','--all'],text=True);(run/'config/complete_environment_lock.txt').write_text(freeze)
    exports=[]
    for p in sorted((BASE/'conda-meta').glob('*.json')):
        r=json.loads(p.read_text());exports.append(dict(name=r['name'],version=r['version'],build=r['build'],url=r.get('url'),sha256=r.get('sha256'),md5=r.get('md5')))
    save(run/'config/conda_base_lock.json',exports)
    (run/'config/conda_explicit_lock.txt').write_text('@EXPLICIT\n'+'\n'.join(r['url']+'#'+r['md5'] for r in exports if r['url'] and r['md5'])+'\n')
    (run/'environment_manifest.txt').write_text('Fixed before docking. Layered environment: exact conda plip_env base + isolated venv --system-site-packages pip lock. Base is read-only; no global environment changed.\nPython '+sys.version+'\nVina '+versions['vina_version']+'\nOpenBabel package 3.1.1, library '+versions['openbabel_library_version']+'\nSee config/software_versions.json, conda_explicit_lock.txt, complete_environment_lock.txt. Vina binary SHA256 recorded; executable is an existing arm64 Mach-O Vina 1.2.7 at /opt/homebrew/bin/vina; installation provenance unknown; exact binary SHA256 identifies it.\n')
    cfg=dict(drugs=DRUGS,blocks=BLOCKS,seeds=SEEDS,exhaustiveness=8,num_modes=9,energy_range=4,cpu=1,concurrent_jobs=4,min_rmsd=1.0,scoring='vina',boxes=read(ROOT/'outputs/docking_boxes.csv'),ligand=dict(source='runs/20260918_1117_chembl_acquisition/processed/compounds.csv',protonation='as-recorded ChEMBL formal charge/tautomer; no pH microstate or tautomer enumeration',stereochemistry='preserve specified; unspecified remains source-unspecified; one deterministic sampled 3D configuration, not a known stereoisomer or racemate ensemble',AddHs=True,embedding='ETKDGv3, randomSeed20260919, numThreads1, enforceChiralityTrue, useRandomCoordsFalse',minimization='MMFF94s maxIters2000; no fallback; convergence audited',pdbqt='Meeko 0.7.1 gasteiger, rigid_macrocycles=True, flexible_amides=False; standard acyclic torsions',seed=20260919),receptor=dict(source='archived PDB ATOM selected chains; model1, altloc blank/A, no waters or heteroatoms',preparation='OpenBabel add hydrogens once, no pH titration, Gasteiger; rigid PDBQT with preserved names, default nonpolar-H merging; all-H PDB for PLIP; coordinates fixed per block'),pose_policy='all returned <=9 per seed, no padding/rerun based on score; denominator actual successfully extracted poses',plip='NOHYDRO=True; NOFIX=True; prepared H retained; raw XML/TXT per pose',unstable_feature_rule='seed frequency range > 0.40 (descriptive preregistered threshold)',policy_revision='All ten drugs newly redocked in common environment; legacy runs historical only, no value matching',external_annotation_allowed=False)
    save(run/'config/docking_config.json',cfg)
    unsup=json.loads((ROOT/'config/ten_drug_unsupervised_preregistered.json').read_text());unsup.update(stage='standardized ten-drug multiseed preregistered',docking=cfg,primary={**unsup['primary'],'cluster_cut_k':3},stability='for each of 5 seeds, same feature schema and metrics; compare aggregate distance ranks (Spearman) and k=3 adjusted Rand; all seeds retained',frequency_denominator='actual evaluable poses, up to45; invalid extraction not zero',unstable_feature_range_threshold=.4)
    save(run/'config/unsupervised_analysis_config.json',unsup)
    pc=configure_plip();shutil.copy2(pc.__file__,run/'config/plip_config.py')
    save(run/'config/plip_effective_config.json',{k:v for k,v in vars(pc).items() if k.isupper() and isinstance(v,(int,float,str,bool,list,dict,type(None)))})
    shutil.copy2(ROOT/'ANALYSIS_CONTRACT.md',run/'config/ANALYSIS_CONTRACT.md')
    (run/'config/ANALYSIS_CONTRACT_ADDENDUM.md').write_text('Primary: same receptor → different drugs. All ten drugs newly redocked under common frozen environment, five computational search seeds; not biological replicates. Legacy runs retained and not combined. This authorized revision supersedes the preflight requirement to recover the legacy environment. Structural results freeze precedes external annotation; no supervised training.\n')
    save(run/'config/PREREGISTRATION_FROZEN.json',{'at':datetime.now().astimezone().isoformat(),'artifacts':{str(p.relative_to(run)):sha(p) for p in (run/'config').iterdir() if p.is_file()},'parent_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()})
    print(run,flush=True)

def prepare():
    from rdkit import Chem
    from rdkit.Chem import AllChem,rdMolDescriptors
    from meeko import MoleculePreparation,PDBQTWriterLegacy
    from openbabel import openbabel as ob
    run=runpath();cfg=json.loads((run/'config/docking_config.json').read_text())
    compounds={r['drug_id']:r for r in read(ROOT/cfg['ligand']['source'])};qc=[]
    for drug in DRUGS:
        r=compounds[drug];m=Chem.MolFromSmiles(r['canonical_smiles']);assert Chem.MolToInchiKey(m)==r['standard_inchikey']
        stereo=Chem.FindMolChiralCenters(m,includeUnassigned=True);unassigned=[i for i,x in stereo if x=='?'];rot=rdMolDescriptors.CalcNumRotatableBonds(m)
        m=Chem.AddHs(m);p=AllChem.ETKDGv3();p.randomSeed=20260919;p.numThreads=1;p.enforceChirality=True;p.useRandomCoords=False
        assert AllChem.EmbedMolecule(m,p)==0;assert AllChem.MMFFHasAllMoleculeParams(m)
        converged=AllChem.MMFFOptimizeMolecule(m,mmffVariant='MMFF94s',maxIters=2000);assert converged==0
        m.SetProp('_Name',drug);sdf=run/'raw/ligands'/f'{drug}.sdf';w=Chem.SDWriter(str(sdf));w.write(m);w.close()
        setup=MoleculePreparation(charge_model='gasteiger',rigid_macrocycles=True,flexible_amides=False).prepare(m)[0]
        txt,ok,err=PDBQTWriterLegacy.write_string(setup);assert ok,err
        (run/'raw/ligands'/f'{drug}.pdbqt').write_text(txt)
        branch=sum(l.startswith('BRANCH ') for l in txt.splitlines());tors=int(next(l.split()[1] for l in txt.splitlines() if l.startswith('TORSDOF')))
        assert branch==tors;assert branch>0 if rot>0 else True
        qc.append(dict(drug_id=drug,molecule_chembl_id=r['molecule_chembl_id'],canonical_smiles=r['canonical_smiles'],standard_inchikey=r['standard_inchikey'],source_file=r['source_file'],source_sha256=r['source_sha256'],source_run=cfg['ligand']['source'],rdkit_rotatable_bonds=rot,pdbqt_branch_count=branch,torsdof=tors,torsion_qc='consistent' if rot==tors else 'definition_difference_requires_review',unassigned_stereocentres=unassigned,stereochemistry_flag='source_unspecified_single_sample_not_known_isomer' if unassigned else 'no_unassigned_atom_stereocentre',formal_charge=Chem.GetFormalCharge(m),mmff_converged=True,pdbqt_sha256=sha(run/'raw/ligands'/f'{drug}.pdbqt')))
        print('PREP ligand',drug,rot,branch,tors,flush=True)
    table(run/'tables/ligand_preparation_qc.csv',qc)
    recs=[]
    for block,(pdb,chains) in BLOCKS.items():
        source=ROOT/f'data/raw/structures/strict3/{pdb}.pdb';lines=[]
        for l in source.read_text().splitlines():
            if l.startswith('ENDMDL'):break
            if l.startswith('ATOM') and l[21] in chains and l[16] in ' A':lines.append(l[:16]+' '+l[17:])
        selected=run/'raw/receptors'/f'{block}_selected.pdb';selected.write_text('\n'.join(lines)+'\nEND\n')
        conv=ob.OBConversion();conv.SetInFormat('pdb');mol=ob.OBMol();assert conv.ReadFile(mol,str(selected));mol.AddHydrogens()
        charge=ob.OBChargeModel.FindType('gasteiger');assert charge.ComputeCharges(mol)
        conv.SetOutFormat('pdb');hfile=run/'raw/receptors'/f'{block}_prepared_H.pdb';assert conv.WriteFile(mol,str(hfile));conv.CloseOutFile()
        conv.SetOutFormat('pdbqt')
        for opt in ['r','c','n']:conv.AddOption(opt,ob.OBConversion.OUTOPTIONS)
        path=run/'raw/receptors'/f'{block}.pdbqt';assert conv.WriteFile(mol,str(path));conv.CloseOutFile()
        assert not any(l.startswith(('BRANCH','ROOT','TORSDOF')) for l in path.read_text().splitlines())
        recs.append(dict(block=block,pdb=pdb,chains=chains,pdbqt_sha256=sha(path),prepared_H_sha256=sha(hfile),source_sha256=sha(source),atoms=mol.NumAtoms(),partial_charge_model='gasteiger',alpha2_provisional=block=='alpha2'))
        print('PREP receptor',block,mol.NumAtoms(),flush=True)
    table(run/'tables/receptor_preparation_qc.csv',recs)
    save(run/'config/PREPARED_INPUTS_FROZEN.json',{'at':datetime.now().astimezone().isoformat(),'files':{str(p.relative_to(run)):sha(p) for folder in ['raw/receptors','raw/ligands'] for p in (run/folder).iterdir() if p.is_file()}})

def dock_job(args):
    drug,block,seed=args;run=runpath();cfg=json.loads((run/'config/docking_config.json').read_text());box=cfg['boxes'][0 if block=='alpha1' else 1]
    job=WORK/'jobs'/drug/block/str(seed);job.mkdir(parents=True,exist_ok=True);out=job/'poses.pdbqt';meta=job/'docking.json'
    if meta.exists():return json.loads(meta.read_text())
    cmd=[str(VINA),'--receptor',str(run/'raw/receptors'/f'{block}.pdbqt'),'--ligand',str(run/'raw/ligands'/f'{drug}.pdbqt'),'--out',str(out)]
    for a in ['center_x','center_y','center_z','size_x','size_y','size_z']:cmd+=['--'+a,str(box[a])]
    for a in ['exhaustiveness','num_modes','energy_range','cpu','min_rmsd','scoring']:cmd+=['--'+a,str(cfg[a])]
    cmd+=['--seed',str(seed)];start=datetime.now().astimezone().isoformat()
    proc=subprocess.run(cmd,text=True,capture_output=True);(job/'vina.log').write_text(proc.stdout+'\n'+proc.stderr)
    poses=len(re.findall(r'^MODEL\s',out.read_text(),re.M)) if out.exists() else 0
    row=dict(drug_id=drug,receptor=block,seed=seed,started_at=start,finished_at=datetime.now().astimezone().isoformat(),returncode=proc.returncode,n_poses=poses,status='success' if proc.returncode==0 and 0<poses<=9 else 'failed',receptor_sha256=sha(run/'raw/receptors'/f'{block}.pdbqt'),ligand_sha256=sha(run/'raw/ligands'/f'{drug}.pdbqt'),command=cmd,output_sha256=sha(out) if out.exists() else None)
    save(meta,row);return row

def docking():
    tasks=[(d,b,s) for d in DRUGS for b in BLOCKS for s in SEEDS];rows=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(dock_job,t) for t in tasks]):
            r=f.result();rows.append(r);print('DOCK',len(rows),'/100',r['drug_id'],r['receptor'],r['seed'],r['status'],r['n_poses'],flush=True)
    table(runpath()/'tables/docking_run_manifest.csv',sorted(rows,key=lambda r:(DRUGS.index(r['drug_id']),r['receptor'],r['seed'])))
    assert all(r['status']=='success' for r in rows)

if __name__=='__main__':
    {'init':initialize,'prepare':prepare,'dock':docking}[sys.argv[1]]()
