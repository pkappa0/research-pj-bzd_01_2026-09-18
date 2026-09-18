"""Evidence-only preflight; deliberately cannot execute docking or clustering."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
from datetime import datetime
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return list(csv.DictReader(path.open()))
def write(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def save(path,obj):path.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')

def check_legacy(root=ROOT):
    base=root/'runs/20260917_1706_structure_mouse_bridge'
    contacts=read(base/'data/processed/reassigned_interactions.csv')
    poses=read(base/'data/raw/lineage/docking_results_4drug.csv')
    old=read(base/'results/tables/candidate_residue_features.csv')
    counts=Counter((r['drug_id'],r['receptor_id']) for r in poses)
    checks=[]
    for row in old:
        subset=[r for r in contacts if all(r[k]==row[k] for k in ['drug_id','receptor','common_position','interaction_type']) and r['mapped_residue']==row['residue_name']]
        n=len({r['pose_id'] for r in subset});den=counts[(row['drug_id'],row['receptor_id'])]
        checks.append(dict(drug_id=row['drug_id'],receptor=row['receptor'],common_position=row['common_position'],interaction_type=row['interaction_type'],archived_frequency=row['frequency'],recounted_frequency=n/den,match=abs(n/den-float(row['frequency']))<1e-12,n_poses=den))
    if not all(r['match'] for r in checks):raise ValueError('Legacy frequency mismatch')
    if len(counts)!=8 or set(counts.values())!={9}:raise ValueError('Legacy pose count mismatch')
    return checks

def gate(environment_evidence, extraction_evidence):
    required={'vina_version','rdkit_version','plip_version','plip_effective_config','openbabel_version'}
    return required <= set(environment_evidence) and all(environment_evidence[k] for k in required) and extraction_evidence is True

def run():
    cfg=json.loads((ROOT/'config/ten_drug_unsupervised_preregistered.json').read_text())
    run=ROOT/'runs'/(datetime.now().astimezone().strftime('%Y%m%d_%H%M')+'_10drug_unsupervised_fingerprint');run.mkdir(exist_ok=False)
    for d in ['tables','reports','config','logs/provenance','logs/tests','logs/docking','logs/PLIP','code']:(run/d).mkdir(parents=True)
    checks=check_legacy();write(run/'tables/legacy_four_drug_frequency_qc.csv',checks)
    fields=['drug_id','molecule_chembl_id','canonical_smiles','standard_inchikey','identity_status','chirality','source_file','source_json_pointer','source_sha256','retrieved_at_utc','chembl_version']
    compounds=read(ROOT/cfg['compound_source']);bydrug={r['drug_id']:r for r in compounds}
    registry=[{k:bydrug[d][k] for k in fields} | {'role':'legacy' if d in cfg['legacy_drugs'] else 'new','docking_status':'existing_only' if d in cfg['legacy_drugs'] else 'not_started_preflight_blocked','stereochemistry_policy':'use recorded SMILES without inventing stereochemistry; unspecified stereocentres require explicit handling before 3D generation'} for d in cfg['drugs']]
    assert len({r['molecule_chembl_id'] for r in registry})==10
    write(run/'tables/compound_structure_registry.csv',registry)
    items=[
      ('receptor structure','documented','6HUP alpha1 full pentamer; 9CTJ provisional alpha2 local','src/strict_poc.py:prepare_receptors'),
      ('receptor preparation','documented','PDB model 0 ATOM protein; remove HETATM; altloc A/blank; simplified atom types; all receptor charges zero; no pH assignment','src/strict_poc.py:prepare_receptors'),
      ('chain selection','documented','alpha1 ABCDE (interface D/C); alpha2 CDE (interface D/E)','src/strict_poc.py:prepare_receptors'),
      ('binding box','documented','archived exact centers and 20 A x 20 A x 20 A sizes; native DZP / gamma2 transferred','outputs/docking_boxes.csv'),
      ('ligand preparation','documented_with_unknown_defaults','ChEMBL SMILES; AddHs; EmbedMolecule randomSeed SHA256(drug) mod (2^31-1); MMFF maxIters500, UFF on exception; RemoveHs; rigid TORSDOF0','src/strict_poc.py:generate_ligands; src/four_drug_poc.py:generate_alprazolam_ligand'),
      ('protonation and charge','documented_limitations','No pH/microstate enumeration; Gasteiger on no-H ligand, charge failure zero fallback; old failure incidence not logged','src/strict_poc.py:generate_ligands'),
      ('software versions','BLOCKER','Historical Vina/RDKit/PLIP/Open Babel versions not found; current Vina 1.2.7 is NOT evidence of historical version','requirements.txt (unpinned, excludes docking dependencies); git initial snapshot'),
      ('docking parameters','documented','exhaustiveness8; num_modes9; energy_range4; seed20260917; CPU/default min_rmsd unspecified','src/strict_poc.py:run_vina; src/four_drug_poc.py:run_alprazolam_docking'),
      ('pose selection','documented','all returned scored MODEL poses; no score-based subset; 72 archived rows, 9 per drug/block','data/raw/lineage/docking_results_4drug.csv'),
      ('PLIP extraction','BLOCKER','CLI -x -t -q --breakcomposite --name report; vendored package/config absent; distance/angle cutoffs and Open Babel protonation defaults cannot be verified','src/strict_poc.py:run_plip; .gitignore'),
      ('raw job evidence','BLOCKER','poses.pdbqt, per-pose PLIP report.xml and Vina logs not present; docking paths in tables refer to former environment; per-pose successful extraction cannot be independently audited','data/raw/strict3/docking and strict4/docking are gitignored'),
      ('residue mapping','documented','receptor-specific chain+residue lookup; alpha D/D, gamma C/E; do not use chain D alone; unknown chains must remain unresolved','runs/20260917_1706_structure_mouse_bridge/results/tables/corrected_residue_mapping.csv'),
      ('common residue key','documented','BZD_SITE_* alpha alignment and BZD_GAMMA2_* gamma alignment, block-qualified','corrected_residue_mapping.csv'),
      ('frequency','documented','unique interacting poses / valid pose count; existing archived candidates recounted without overwriting','src/structure_mouse_bridge.py:_features'),
      ('interaction types','documented_partial_thresholds','hydrophobic_interaction, hydrogen_bond, water_bridge, salt_bridge, pi_stack, pi_cation, halogen_bond, metal_complex parser; preserve all actually observed','src/strict_poc.py:run_plip'),
      ('geometry','documented_partial_thresholds','pi centdist A, angle degrees, offset A, P/T; raw-contact mean/median/std; P/T unique pose frequency; no-contact geometry NA; thresholds unknown','src/structure_mouse_bridge.py:_features')]
    audit=[dict(item=a,status=b,finding=c,evidence=d) for a,b,c,d in items];write(run/'tables/protocol_audit.csv',audit)
    save(run/'config/unsupervised_analysis_config.json',cfg)
    paths=['ANALYSIS_CONTRACT.md','requirements.txt','.gitignore','src/strict_poc.py','src/four_drug_poc.py','src/structure_mouse_bridge.py','outputs/docking_boxes.csv',cfg['compound_source'],cfg['structural_source']+'/data/processed/reassigned_interactions.csv',cfg['structural_source']+'/data/raw/lineage/docking_results_4drug.csv',cfg['structural_source']+'/results/tables/corrected_residue_mapping.csv',cfg['structural_source']+'/results/tables/candidate_residue_features.csv','data/raw/strict3/receptors/alpha1_beta3_gamma2.pdbqt','data/raw/strict3/receptors/alpha2_beta3_gamma2_local_9CTJ.pdbqt']
    save(run/'logs/provenance/source_inventory.json',[dict(path=f,sha256=sha(ROOT/f)) for f in paths])
    (run/'logs/provenance/ANALYSIS_CONTRACT.md').write_bytes((ROOT/'ANALYSIS_CONTRACT.md').read_bytes())
    for name in ['src/ten_drug_preflight.py','tests/test_ten_drug_preflight.py']:(run/'code'/Path(name).name).write_bytes((ROOT/name).read_bytes())
    for d in ['docking','PLIP']:save(run/'logs'/d/'execution_status.json',dict(executed=False,reason='historical_environment_and_job_evidence_unavailable'))
    qc=dict(status='BLOCKED_BEFORE_NEW_DOCKING',gate_passed=gate({},False),legacy_frequency_rows_checked=len(checks),legacy_frequency_match=all(r['match'] for r in checks),legacy_pose_rows=72,new_docking_jobs_executed=0,structural_clustering_executed=False,pca_executed=False,external_annotations_loaded=False,config_frozen_before_analysis=True,compound_structures_available=10,ten_drug_matrix_available=False,missing_values_imputed=False)
    save(run/'reports/QC_REPORT.json',qc)
    report='''# 10-drug unsupervised fingerprint — preflight blocked

This is an audit/preregistration run, NOT a completed ten-drug analysis.

## Audit outcome / 不足項目
追加6剤Dockingは未実行。既存4剤と同一環境をrepositoryから再現できないため、依頼Step1のゲートで停止した。Vina/RDKit/PLIP/Open Babelの旧実行version、vendored PLIPの実効config、元Vina log/pose出力および全poseのPLIP reportが必要。requirements.txtは版固定なしでこれらを含まず、work/とdocking/はgitignore対象。現在Vina v1.2.7が存在することは旧versionの証明ではない。ローカルvenvにはRDKit/PLIP/Open Babelもない。別versionをインストールして同一pipelineと呼ぶことはしていない。

## 確認できた条件
6HUP ABCDE / 9CTJ CDE、保存済みPDBQT、box、seed20260917、exhaustiveness8、num_modes9、energy_range4、全返却pose採択を確認。boxは20 Å角で共通、中心はoutputs/docking_boxes.csvをそのまま使用する設計。ligandはhash seed→AddHs/EmbedMolecule→MMFF（例外時UFF）→RemoveHs→Gasteiger簡易writer→TORSDOF0。受容体電荷は0、pH指定なし。EmbedMoleculeの未指定defaultがversion依存であるため、コメントのETKDGだけから実効設定を推測しない。新薬のために物理的に改善した前処理へ置き換えると旧4剤との条件差になる。

既存4剤100 candidate feature行は元contact表から頻度を再計数して完全一致、72 pose行（各drug/block 9）を確認。ただしこれは保存CSV間の整合性であり、存在しないPLIP元reportによる全pose成功監査の代わりではない。未観測を無条件0とした10剤表は作成していない。

## 化合物構造
10剤のChEMBL ID、SMILES、InChIKey、取得元hashを既存取得runから抽出。追加6剤も別sourceへ切替不要。lorazepam / temazepam / zopicloneの記録SMILESには立体指定がなく、単一生成conformerをracemateの代表と断定できない。未指定立体・microstateの扱いは旧環境回収後も明示する必要がある。

## Preregistration
config/unsupervised_analysis_config.jsonは新規Docking/PCA/clustering前に保存。Bは各残基のtypeをcollapseしたpose union frequency、Cはtype別頻度＋P/T、EはBのα1/α2連結。A binary、D geometry付きも事前定義。primaryはEuclidean＋average linkage、k=3、center-only PCA。secondaryのcosineとk=2/4も事前指定して全結果を報告し、外部annotationに合う条件を選ばない。PCA/clusteringは各表の全10剤共通観測列のみ、pairwise similarityは共通feature数も保存。0補完は禁止。heatmapはtype別C行とprimary B/E dendrogramを表示する計画。geometryには別単位と欠損maskを保持する。

ここでは薬理activityやmouse phenotypeをロードしていない。構造結果hashをfreezeしてからannotation処理を行うゲートを設定。既存会話で知られている薬理知識から解析条件を最適化しない。単一seedのmetric感度を独立Docking再現性の証拠とはしない。

## A. Structural pattern
**indeterminate / 未評価**。追加6剤のIFPがなく、10剤cluster/PCAは実行できない。4剤の一致監査を10剤の結果へ一般化しない。

## B. External correspondence
**indeterminate / 未評価**。構造結果未固定のため、外部annotationの統合は未実行。

## C. ML expansion rationale
**indeterminate**。現段階の不足は再現環境/provenanceの問題であり、科学的仮説を支持も反証もしない。10剤matrix、similarity、cluster assignments、PCA図は未生成。欠損6剤を0で埋めた図やダミー図は作らない。

## 再開に必要な証拠
旧実行環境のpackage lock/conda export/container digest、Vina --version記録、work/plip_vendorと実効config、元strict3/strict4 docking job directories（Vina log、poses.pdbqt、per-pose PLIP XML）。これらが回収不能なら、完全同条件の追加解析という条件を満たせない。条件変更を伴う別設計は、ユーザーによる研究方針の変更を受けてから別runで扱う。
'''
    (run/'reports/FINAL_REPORT.md').write_text(report)
    save(run/'run_manifest.json',dict(status='BLOCKED_BEFORE_NEW_DOCKING',parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),created_at=datetime.now().astimezone().isoformat(),config_sha256=sha(run/'config/unsupervised_analysis_config.json'),artifact_sha256={str(p.relative_to(run)):sha(p) for p in sorted(run.rglob('*')) if p.is_file()}))
    print(run)

if __name__=='__main__':run()
