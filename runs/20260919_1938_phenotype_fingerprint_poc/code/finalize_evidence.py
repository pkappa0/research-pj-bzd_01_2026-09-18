from pathlib import Path
import csv,json,hashlib,subprocess,shutil,xml.etree.ElementTree as ET
from datetime import datetime
root=Path.cwd();w=root.parent/'phenotype_poc';r=Path((w/'active_run.txt').read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,o):p.write_text(json.dumps(o,indent=2,ensure_ascii=False)+'\n')
def table(p,rows):
    with p.open('w') as h:
        z=csv.DictWriter(h,fieldnames=list(rows[0]));z.writeheader();z.writerows(rows)
paths=subprocess.check_output(['git','ls-files'],text=True).splitlines();audit=[]
for p in paths:
    if not any(s in p.lower() for s in ['phenotype','structure_mouse_bridge','process_phenotype']):continue
    f=root/p
    status='derived or lineage duplicate; not re-counted as independent evidence'
    if p.endswith('data/processed/mouse_phenotype_raw.csv'):status='canonical 44-row archived source; all rows retained'
    elif 'template' in p or 'fixtures' in p:status='template/test fixture; not scientific observations'
    elif p.endswith('phenotype_processed.csv'):status='header-only; no observations'
    audit.append(dict(path=p,sha256=sha(f),audit_disposition=status))
table(r/'tables/repository_phenotype_audit.csv',audit)
queries=['zolpidem zopiclone zaleplon diazepam mice rotarod traction locomotor comparative','lorazepam clonazepam midazolam temazepam mice rotarod muscle relaxation comparative','temazepam midazolam alprazolam mouse rotarod locomotor righting reflex comparative ED50','temazepam mice locomotor rotarod pharmacological profile','temazepam mouse righting reflex locomotor activity pharmacology study','temazepam mice locomotor -site:benchchem.com -site:wikipedia.org','clonazepam mice rotarod locomotor impairment study PubMed']
save(r/'raw/search_log.json',dict(date='2026-09-19',search_engine='web search; primary literature inclusion only',queries=queries,scope='targeted feasibility search, not systematic review',identity_exclusions=['eszopiclone is not zopiclone','tetrazepam is not temazepam'],fulltext_attempts={'Sanger1996':'publisher HTTP403; abstract only','Stanley2005':'publisher restricted access; abstract only','Nishino2008':'PDF downloaded and Table1 rendered','Tanaka2008':'PDF downloaded and Table1 rendered','Shi2024':'PMC indexed fulltext/Table1 and methods verified'}))
raw=list(csv.DictReader((r/'tables/phenotype_evidence_raw.csv').open()));shutil.copy2(r/'tables/phenotype_evidence_raw.csv',r/'raw/source_evidence_transcription.csv')
sources=[]
for study in dict.fromkeys(z['study_id'] for z in raw):
    ss=[z for z in raw if z['study_id']==study];p=w/f'{study}.pdf'
    sources.append(dict(study_id=study,source_url=ss[0]['source_url'],doi=ss[0]['doi'],pmid=ss[0]['pmid'],accessed_date='2026-09-19',source_file_sha256=sha(p) if p.exists() else 'NA',extraction_scope=';'.join(sorted(set(z['source_table_figure'] for z in ss))),provenance='archived extraction + original PDF audit' if study in ['Nishino2008','Tanaka2008'] else 'archived record' if study in ['Bourin1992','Bayley1996','Henauer1984','Lopez1988'] else 'primary public article / abstract',numeric_rows=sum(z['raw_value']!='NA' for z in ss)))
table(r/'raw/source_manifest.csv',sources)
table(r/'tables/out_of_scope_evidence.csv',[
dict(source='PMID10065909; DOI10.1177/026988119801200405',finding='anxiolytic comparisons of triazolam/zopiclone/zolpidem/zaleplon',disposition='anxiety outside four primary axes; no anxiety score imported'),
dict(source='PMID1425360',finding='rat rotarod PDD50: midazolam/lorazepam/diazepam',disposition='rat dose-based endpoint, excluded from mouse primary; no cross-species or dose-only comparison'),
dict(source='PMC3871832',finding='rat rotarod and sleep including eszopiclone',disposition='rat study; eszopiclone not substituted for zopiclone'),
dict(source='PMID2875452',finding='mouse benzodiazepine comparative study including lorazepam/clonazepam',disposition='abstract numeric values refer to chlordesmethyldiazepam outside panel; not assigned to target drugs'),
dict(source='DOI10.1038/sj.bjp.0701513',finding='temazepam mouse light-dark disinhibition',disposition='anxiety-related endpoint outside four primary axes'),
dict(source='PMC11559241',finding='Y-maze/cognition after repeated exposure; remimazolam comparator',disposition='cognition outside axes; remimazolam outside fixed panel')])
norm=list(csv.DictReader((r/'tables/phenotype_normalized_within_study.csv').open()));assert len(raw)==69
assert len({z['drug'] for z in norm if z['tier']=='1'})==2
assert all(z['value']=='NA' for z in norm if z['raw_value']=='NA')
expected={'diazepam':{2:[1,3,4,3],5:[2,4,5,4],10:[4,7,9,8],20:[9,10,10,9]},'triazolam':{.5:[2,3,3,9],1:[3,4,4,2],2:[5,8,7,5],5:[8,9,9,4]}}
for z in raw:
    if z['study_id']=='Nishino2008':assert float(z['raw_value'])==expected[z['drug']][float(z['dose'])][[15,30,60,90].index(int(float(z['observation_time'])))]
joined=list(csv.DictReader((r/'tables/fingerprint_phenotype_joined.csv').open()));assert len(joined)==69
xf=[k for k in joined[0] if 'any_contact_frequency' in k];assert len(xf)==48
base=json.loads((w/'baseline.json').read_text());assert all(sha(root/p)==h for p,h in base.items())
tests=ET.parse(w/'tests.xml');ts=list(tests.iter('testsuite'));counts={k:sum(int(s.get(k,0)) for s in ts) for k in ['tests','failures','errors']};assert counts['errors']==counts['failures']==0
q=json.loads((r/'reports/QC_REPORT.json').read_text());q.update(tests=counts,all_Nishino_counts_checked_against_rendered_table=True,source_PDF_tables_visually_reviewed=True,raw_processed_separated=True,raw69_joined69_no_join_expansion=True,source_audit_records=len(audit));save(r/'reports/QC_REPORT.json',q)
for name in ['tests.xml','tests.log']:shutil.copy2(w/name,r/'reports'/name)
shutil.copy2(root/'tests/test_phenotype_fingerprint_poc.py',r/'code/test_phenotype_fingerprint_poc.py')
shutil.copy2(__file__,r/'code/finalize_evidence.py')
report=r/'reports/FINAL_REPORT.md';report.write_text(report.read_text()+f'\n## Validation\n{counts["tests"]} tests passed; archived {len(base)} files unchanged. All 32 Nishino count cells checked against rendered original Table1. Full 48-feature X retained. Source tables, search/access limits and duplicate/template exclusions are recorded.\n')
save(r/'run_manifest.json',dict(completed_at=datetime.now().astimezone().isoformat(),artifacts={str(p.relative_to(r)):sha(p) for p in sorted(r.rglob('*')) if p.is_file() and p.name!='run_manifest.json'}))
print(json.dumps(q,indent=2))
