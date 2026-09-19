"""All-pose PLIP extraction using fixed prepared protein/ligand hydrogens."""
from .standardized_redocking import *
import xml.etree.ElementTree as ET

def complex_pdb(protein, molecule, conf):
    from rdkit import Chem
    lines=[]
    for line in protein.splitlines():
        if line.startswith(('ATOM','HETATM')):lines.append(f'{line[:6]}{len(lines)+1:5d}{line[11:]}')
    offset=len(lines)
    for i,a in enumerate(molecule.GetAtoms()):
        info=Chem.AtomPDBResidueInfo();info.SetName(f'{a.GetSymbol()}{i+1}'[:4].ljust(4));info.SetResidueName('LIG');info.SetResidueNumber(1);info.SetChainId('Z');info.SetIsHeteroAtom(True);a.SetMonomerInfo(info)
    ligand=Chem.MolToPDBBlock(molecule,confId=conf)
    for line in ligand.splitlines():
        if line.startswith(('ATOM','HETATM')):lines.append(f'{line[:6]}{int(line[6:11])+offset:5d}{line[11:]}')
        elif line.startswith('CONECT'):
            ids=[int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
            lines.append('CONECT'+''.join(f'{i+offset:5d}' for i in ids))
    return '\n'.join(lines)+'\nEND\n'

def worker(drug,block,seed):
    from rdkit import Chem
    from meeko import PDBQTMolecule,RDKitMolCreate
    from plip.structure.preparation import PDBComplex
    from plip.exchange.report import StructureReport
    configure_plip();run=runpath();job=WORK/'jobs'/drug/block/str(seed)
    if (job/'plip_done.json').exists():return
    dock=json.loads((job/'docking.json').read_text());assert dock['status']=='success'
    text=(job/'poses.pdbqt').read_text();mols=RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text,skip_typing=True));assert len(mols)==1 and mols[0] is not None
    mol=mols[0];assert mol.GetNumConformers()==dock['n_poses']
    protein=(run/'raw/receptors'/f'{block}_prepared_H.pdb').read_text();results=[];contacts=[]
    for i in range(mol.GetNumConformers()):
        poseid=f'{seed}:{i+1}';out=job/f'pose_{i+1:02}';out.mkdir(exist_ok=True)
        writer=Chem.SDWriter(str(out/'ligand.sdf'));writer.write(mol,confId=i);writer.close()
        pdb=out/'complex.pdb';pdb.write_text(complex_pdb(protein,mol,i));inputhash=sha(pdb)
        start=datetime.now().astimezone().isoformat()
        try:
            c=PDBComplex();c.output_path=str(out)+'/';c.load_pdb(str(pdb));c.analyze();rep=StructureReport(c);rep.write_xml();rep.write_txt()
            xml=ET.parse(out/'report.xml');sites=[s for s in xml.findall('bindingsite') if s.findtext('identifiers/hetid')=='LIG' and s.findtext('identifiers/chain')=='Z'];assert len(sites)==1,'Missing or ambiguous LIG:Z binding site'
            site=sites[0];n=0
            for group in site.findall('interactions/*'):
                for node in group:
                    attrs={child.tag:child.text.strip() if child.text else None for child in node};attrs.update(node.attrib)
                    contacts.append(dict(drug_id=drug,receptor=block,seed=seed,pose_rank=i+1,pose_id=poseid,interaction_type=node.tag,residue_chain=attrs.get('reschain'),residue_number=attrs.get('resnr'),residue_name=attrs.get('restype'),raw_attributes=attrs,xml_relative=f'{drug}/{block}/{seed}/pose_{i+1:02}/report.xml',xml_sha256=sha(out/'report.xml')));n+=1
            row=dict(drug_id=drug,receptor=block,seed=seed,pose_rank=i+1,pose_id=poseid,status='success',n_interactions=n,xml_sha256=sha(out/'report.xml'),complex_input_sha256=inputhash,ligand_sdf_sha256=sha(out/'ligand.sdf'),started_at=start,finished_at=datetime.now().astimezone().isoformat(),ligand_heavy_atoms=mol.GetNumHeavyAtoms())
        except Exception as exc:
            row=dict(drug_id=drug,receptor=block,seed=seed,pose_rank=i+1,pose_id=poseid,status='failed',error=repr(exc),complex_input_sha256=inputhash)
        save(out/'pose_status.json',row);results.append(row)
        # Reconstructable large protein+ligand concatenation is scratch only.
        # Retain prepared protein, per-pose ligand SDF, XML/TXT, hash and code.
        pdb.unlink()
        print('PLIP',drug,block,seed,i+1,row['status'],flush=True)
    save(job/'pose_manifest.json',results);save(job/'interactions.json',contacts);save(job/'plip_done.json',{'n':len(results),'success':sum(r['status']=='success' for r in results)})

def all_jobs():
    tasks=[(d,b,s) for d in DRUGS for b in BLOCKS for s in SEEDS]
    def call(t):
        job=WORK/'jobs'/t[0]/t[1]/str(t[2]);log=job/'plip.log'
        import time
        while not (job/'docking.json').exists():time.sleep(1)
        with log.open('a') as f:r=subprocess.run([sys.executable,'-m','src.standardized_plip','worker',*map(str,t)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
        if r.returncode:raise RuntimeError(f'PLIP job failed {t}; see {log}')
        return t
    with ThreadPoolExecutor(max_workers=4) as pool:
        for i,f in enumerate(as_completed([pool.submit(call,t) for t in tasks]),1):print('PLIP_JOB',i,'/100',f.result(),flush=True)
    poses=[];contacts=[]
    for t in tasks:
        job=WORK/'jobs'/t[0]/t[1]/str(t[2]);poses+=json.loads((job/'pose_manifest.json').read_text());contacts+=json.loads((job/'interactions.json').read_text())
    for i,r in enumerate(contacts):r['source_contact_id']=i
    table(runpath()/'tables/pose_manifest.csv',poses);table(runpath()/'tables/plip_interactions_raw.csv',contacts)
    save(WORK/'all_contacts.json',contacts);save(WORK/'all_poses.json',poses)
    assert all(r['status']=='success' for r in poses), 'Extraction failure: do not infer absence'

if __name__=='__main__':
    if sys.argv[1]=='worker':worker(sys.argv[2],sys.argv[3],int(sys.argv[4]))
    else:all_jobs()
