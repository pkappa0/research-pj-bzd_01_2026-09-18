"""Package immutable raw evidence; validate data and frozen outputs; write report."""
from .standardized_redocking import *
from collections import Counter,defaultdict
import tarfile,xml.etree.ElementTree as ET
import numpy as np

def main():
    run=runpath();f=json.loads((run/'STRUCTURAL_RESULTS_FROZEN.json').read_text());assert all(sha(run/p)==h for p,h in f['artifacts'].items())
    provenance=json.loads((run/'reports/EXTERNAL_ANNOTATION_PROVENANCE.json').read_text());assert provenance['started_at']>f['frozen_at']
    dock=read(run/'tables/docking_run_manifest.csv');poses=read(run/'tables/pose_manifest.csv');raw=read(run/'tables/plip_interactions_raw.csv');freq=read(run/'tables/fingerprint_seed_level.csv')
    assert len(dock)==100 and all(r['status']=='success' and 0<int(r['n_poses'])<=9 for r in dock)
    expected_pose_count=sum(int(r['n_poses']) for r in dock)
    assert len(poses)==expected_pose_count and all(r['status']=='success' for r in poses)
    assert len({(r['drug_id'],r['receptor'],r['pose_id']) for r in poses})==expected_pose_count
    for b in BLOCKS:assert len({r['receptor_sha256'] for r in dock if r['receptor']==b})==1
    for r in dock:
        cmd=json.loads(r['command']);assert cmd[cmd.index('--exhaustiveness')+1]=='8';assert cmd[cmd.index('--num_modes')+1]=='9';assert cmd[cmd.index('--energy_range')+1]=='4';assert int(r['seed']) in SEEDS
    totals=Counter()
    for r in dock:totals[(r['drug_id'],r['receptor'])]+=int(r['n_poses'])
    for r in freq:assert 0<int(r['n_poses'])<=9 and abs(float(r['frequency'])-int(r['count'])/int(r['n_poses']))<1e-12
    for r in read(run/'tables/fingerprint_aggregate_long.csv'):
        assert int(r['total_evaluable_poses'])==totals[(r['drug_id'],r['receptor'])]
        if r['feature_type']=='frequency':assert abs(float(r['value'])-int(r['total_count'])/int(r['total_evaluable_poses']))<1e-12
        elif int(r['geometry_n_observed'])==0:assert r['value']=='NA'
    a=read(run/'tables/fingerprint_10drug_frequency_alpha1.csv');b=read(run/'tables/fingerprint_10drug_frequency_alpha2.csv');e=read(run/'tables/fingerprint_10drug_multireceptor.csv')
    assert all({**ra,**rb}==re for ra,rb,re in zip(a,b,e))
    # Source XML identity and chemistry extraction completeness, one ligand/pose.
    for r in poses:
        p=WORK/'jobs'/r['drug_id']/r['receptor']/r['seed']/f"pose_{int(r['pose_rank']):02}"/'report.xml';assert sha(p)==r['xml_sha256'];xml=ET.parse(p);assert xml.findtext('plipversion')=='3.0.0'
    # Lossless raw archives. No pose selection, all logs and extraction reports.
    archive_rows=[]
    for d in DRUGS:
        path=run/'raw/archives'/f'{d}.tar.gz';base=WORK/'jobs'/d
        with tarfile.open(path,'w:gz') as tar:
            for p in sorted(base.rglob('*')):
                if p.is_file():tar.add(p,arcname=str(p.relative_to(WORK/'jobs')),recursive=False)
        with tarfile.open(path) as tar:
            for m in tar.getmembers():
                if m.isfile():
                    data=tar.extractfile(m).read();original=WORK/'jobs'/m.name;assert hashlib.sha256(data).hexdigest()==sha(original)
                    archive_rows.append(dict(archive=str(path.relative_to(run)),member=m.name,bytes=m.size,sha256=hashlib.sha256(data).hexdigest()))
    table(run/'tables/raw_archive_manifest.csv',archive_rows)
    for name in ['prepare.log','plip_initial.log','docking_progress.log','plip_progress.log','analysis.log','annotation.log','tests.log','tests.xml']:
        p=WORK/name
        if p.exists():shutil.copy2(p,run/'logs'/name)
    for name in ['standardized_redocking.py','standardized_plip.py','standardized_structure_analysis.py','standardized_external.py','standardized_presentation.py','standardized_finalize.py']:
        shutil.copy2(ROOT/'src'/name,run/'code'/name)
    shutil.copy2(ROOT/'tests/test_standardized_redocking.py',run/'code/test_standardized_redocking.py')
    # Mechanical descriptive ranking; all original loadings/contrasts retained.
    loadings=read(run/'tables/pca_loadings.csv');tops=[]
    for rep in ['B_alpha1','B_alpha2','E']:
        for pc in ['PC1','PC2']:
            selected=sorted([r for r in loadings if r['representation']==rep],key=lambda r:(-abs(float(r[pc+'_loading'])),r['feature']))[:10]
            tops += [dict(r,component=pc,rank=i+1) for i,r in enumerate(selected)]
    table(run/'tables/pca_loading_top10_by_pc.csv',tops)
    contrasts=read(run/'tables/cluster_feature_contrasts.csv');table(run/'tables/cluster_contrast_top10.csv',[dict(r,rank=i+1) for rep in ['B_alpha1','B_alpha2','E'] for i,r in enumerate(sorted([r for r in contrasts if r['representation']==rep],key=lambda r:(-float(r['cluster_mean_range']),r['feature']))[:10])])
    baseline=json.loads((WORK/'baseline.json').read_text());assert all(sha(ROOT/p)==h for p,h in baseline.items())
    tests=ET.parse(WORK/'tests.xml').getroot();suites=list(tests.iter('testsuite'));counts={k:sum(int(s.get(k,0)) for s in suites) for k in ['tests','failures','errors','skipped']};assert counts['failures']==counts['errors']==0
    unmapped=[r for r in read(run/'tables/unmapped_contacts.csv') if r.get('drug_id')]
    summary=json.loads((run/'reports/STRUCTURAL_SUMMARY.json').read_text())
    if summary['unmapped_contacts']!=len(unmapped):
        save(run/'logs/STRUCTURAL_SUMMARY_INITIAL.json',summary)
        summary['unmapped_contacts']=len(unmapped);summary['metadata_count_correction']='Excluded status=none CSV placeholder from row count. No protected table/figure/cluster/PCA changes.';save(run/'reports/STRUCTURAL_SUMMARY.json',summary)
    robustness=read(run/'tables/feature_seed_robustness.csv');counts_unstable=Counter(r['representation'] for r in robustness if r['unstable_feature']=='True')
    qc=dict(status='PASS_WITH_METHOD_LIMITATIONS',tests=counts,docking_success=100,docking_jobs=100,poses=expected_pose_count,plip_success=len(poses),raw_contacts=len(raw),unmapped_contacts=len(unmapped),xml_reports_archived=sum(r['member'].endswith('/report.xml') for r in archive_rows),archived_members=len(archive_rows),all_archived_member_hashes_verified=True,five_seeds_fixed=True,same_environment=True,same_prepared_receptor_per_block=True,same_box_and_parameters=True,torsion_branch_equality=True,frequency_denominator_range=[min(totals.values()),max(totals.values())],seed_denominator_distribution=dict(Counter(int(r['n_poses']) for r in dock)),NA_not_zero=True,blocks_concatenated_not_averaged=True,structural_freeze_verified_after_annotation=True,external_loaded_after_freeze=True,existing_files_unchanged=len(baseline),unstable_type_frequency_drug_features=counts_unstable['C'],unstable_union_frequency_drug_features=counts_unstable['B'],visual_review='frozen original figures retained; reviewed layout-only derivatives supplied',method_limitations=['alpha2 provisional local construct','beta3 structure vs beta2 pharmacology','single initial conformer per ligand, fixed rings','source-unspecified stereochemistry for three compounds','no protonation microstate ensemble','OpenBabel full-protein kekulization warnings; binding-site mapping/types separately audited','single receptor coordinates per block','n=10; sparse comparable external data'],ML_rationale='indeterminate')
    save(run/'reports/QC_REPORT.json',qc)
    report(run,summary,provenance,qc,counts_unstable)
    manifest=dict(status='completed_with_limitations',parent_commit=json.loads((run/'config/PREREGISTRATION_FROZEN.json').read_text()).get('parent_commit') or subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),completed_at=datetime.now().astimezone().isoformat(),structural_freeze_sha256=sha(run/'STRUCTURAL_RESULTS_FROZEN.json'),artifact_sha256={str(p.relative_to(run)):sha(p) for p in sorted(run.rglob('*')) if p.is_file() and p.name!='run_manifest.json'})
    save(run/'run_manifest.json',manifest);print(json.dumps(qc,indent=2),flush=True)

def report(run,summary,external,qc,unstable):
    stability=summary['seed_stability'];lines=['# 10-drug standardized redocking / unsupervised fingerprint PoC','',
      '10剤すべてを新しい共通環境で再Dockingした正式dataset。旧4剤の数値は継ぎ足さず、旧runはhistorical PoCとして保持した。主軸は **same receptor → different drugs**。薬理予測モデルは学習していない。',
      '', '## Execution and environment',
      'Python 3.12.13、RDKit 2025.09.6、AutoDock Vina 1.2.7、Meeko 0.7.1、PLIP 3.0.0。Open Babel conda packageは3.1.1、実ライブラリのOBReleaseVersionは3.1.0であり両方を記録。その他依存packageはsoftware_versions.jsonと完全lockに保存。Vinaは既存arm64 binaryの実SHA256を保存し、インストール元は未確認。',
      f"100/100 Docking成功、{qc['poses']} pose、{qc['plip_success']}/{qc['poses']} PLIP成功。全drug/blockで5 seed、各seed最大9 pose。aggregate分母は実数{qc['frequency_denominator_range']}（drug/blockごと）で、9未満のposeを補完しない。raw contacts {qc['raw_contacts']}行。全XMLと元pose PDBQT・logをdrug別tar.gzへ可逆保存し、全member hashを照合。",
      'Seeds: 2026091901–2026091905。exhaustiveness8、num_modes9、energy_range4、CPU1/job、min_rmsd1.0、Vina scoring。boxは旧記録のcenter/sizeを固定、各辺20 Å。各blockで同一prepared receptor hashを100job manifestから確認。',
      'ChEMBL ID/SMILES/InChIKeyを既存取得runから利用。ETKDGv3 seed20260919、AddHs、MMFF94s最大2000反復（全剤収束）、Meeko Gasteiger PDBQT。TORSDOF/BRANCHは全剤一致（1–4）。RDKitとの差はOHやcarbamate、非等価置換tertiary amide等の定義差で、ligand_rotatable_bond_audit.csvに実結合を保存。旧TORSDOF0と異なり全剤を剛体にはしていない。',
      '受容体は6HUP ABCDE、9CTJ CDE。Open BabelでHを一度付加しGasteiger/rigid PDBQTを固定。PLIPはNOHYDRO/NOFIXを固定し同じprepared Hを使用。PLIP実効config・versionを全poseで統一。原子型に関する全蛋白kekulization警告は隠さずprepare.logに保存。',
      '', '## Dataset / analytical contract',
      'A: residue any-contact binary。B: residue単位のtypeをcollapseしたunique-pose union frequency。C: residue×type＋P/T頻度。D: C＋contact条件付きcentroid distance/angle/offset。E: [α1 B | α2 B]。α1 21列、α2 27列、E48列。rawはtype頻度118列＋geometry30列、geometry欠損はNA。観測なし頻度0と混同しない。geometryはcontact-row平均でありpose平均と同じではない。',
      '全観測interaction typeを保持（hydrophobic、π、H-bond、halogen、salt bridge、pi-cation等）。6指定残基を含む修正済みmappingを新prepared PDBへ照合し、未解決contactは0行。STRUCTURAL_SUMMARYの最初の集計で空表placeholderを1行と数えたmetadataを0に訂正したが、freeze対象の表/cluster/PCA/図は変更していない。',
      'PrimaryはEuclidean＋average linkage、k=3、center-only PCA。PCAのzero-variance列除外とSVD符号規則を記録。cosineとbinary Jaccard、事前指定k2/4、C/D表現はsecondaryとして全結果保存。欠損列の採否は表ごと全10剤共通観測列で固定、pairwise比較可能feature数も保存。今回はprimary B/Eに欠損なし。distance matrixのEuclideanは距離で、小さいほど類似。',
      '', '## A. Structural pattern',
      'α1固定blockではzolpidem、zopicloneがそれぞれ別群、他8剤が同群。α2固定blockではdiazepam＋zolpidem、zaleplon単独、他7剤の群。これは各block内の薬剤比較であり、受容体間優劣ではない。Eでは7剤の群／zolpidem＋zaleplon／zopiclone単独に分かれた。k=3は事前指定で、自然な真のcluster数を証明したものではない。',
      '', '| Representation | Features | PC1 | PC2 |','|---|---:|---:|---:|']
    for rep,r in summary['representations'].items():lines.append(f"| {rep} | {r['features']} | {r['PC1']:.2%} | {r['PC2']:.2%} |")
    lines+=['', 'PCA loadingの上位候補は機械的に各PC10件保存し、全loadingも保持。α1ではγ2 PHE77/TYR58/ASN60、α2ではcommon α TYR210、γ2 PHE77、common α PHE100/TYR160などの頻度が大きいloadingを持つ。これらは残基因果効果ではない。小さいdrug range ≤2/9のfeatureをmechanistic findingとして強調しない。',
      '', '## B. Computational robustness','| Representation | Seed vs aggregate distance Spearman | k=3 ARI |','|---|---:|---:|']
    for rep in summary['representations']:
        rr=[r for r in stability if r['representation']==rep];sp=[r['aggregate_distance_spearman'] for r in rr];ari=[r['aggregate_cluster_ARI'] for r in rr];lines.append(f'| {rep} | {min(sp):.3f}–{max(sp):.3f} | {min(ari):.3f}–{max(ari):.3f} |')
    lines += ['',f"Eのk=3所属は5 seed中4でaggregateと一致、1 seedではdiazepamが変化。α2の所属はより変動する。range>0.40の事前閾値では、Cのdrug-feature {unstable['C']}/1180、Bのdrug-feature {unstable['B']}/480をunstableとした。小さい割合だけで全体の頑健性を断定しない。",
      'Representation依存性もある。α1 B対CのARIは0.646、α2 B対Cは1.000。一方binaryやgeometry付きDはBと異なる分割になる（例：α1 AのARI −0.229、α1 D 0.338）。E cosine k3もEuclideanとARI0.643で一致しない。したがって全表現・metricに普遍的なclusterとは言えない。seedの一致は探索計算の感度であり、生物学的独立replicateや統計的信頼区間ではない。',
      '', '## C. External correspondence',
      f"構造freeze日時: {external['structural_freeze_at']}。外部annotation処理開始: {external['started_at']}。freezeの全hashをannotation後にも検証。",
      'Matched pharmacologyはdiazepam/alprazolam/zolpidemの3剤、5 assay context・10 activityのみ。Kiをcontextごとに保持し、functional potency/efficacyのmatched値はNA。unpairedな既存活動は別logへ保存し、群代表値として転用しない。β2薬理とβ3構造はcontext_mismatch。追加6剤とtriazolamのmatched薬理は今回の既存filterではNAであり、文献全体にデータが存在しないという意味ではない。',
      '対応を支持しない例もある。α1でdiazepam/alprazolamは同clusterでもKi14/0.8 nM。diazepam/zolpidemは別clusterでも同一assay Ki16/19 nM。α2でdiazepam/zolpidemは同clusterでもKi20/156 nM。近い構造なら薬理も近いという単純な関係は支持されず、3剤の少数pairから全体関連を検証できない。',
      '既存P1 mouse evidenceはdiazepam/triazolam/zolpidemの3剤39行を保持。dose/time/route/endpoint/studyを固定して比較可能な8条件はdiazepam–triazolamの1薬剤対のみ。Nishino2008のrotarod同用量2/5 mg/kg・各15/30/60/90分では、同E clusterでもtriazolamのpositive countが7条件で高く1条件で同値。8条件を8独立薬剤対と数えない。他studyや未報告条件・qualitative記録は混合しない。mouse phenotype群との一般的対応は判断不能。',
      'freeze後のChEMBL SMILES ring-topology annotationでは、Eの7剤群が7員環2NのBZD-like群に一致した。zolpidem/zaleplonは同群、zopicloneは別群。ただしα1/α2単独ではこの分類が完全再現されず、化学骨格とfingerprintの対応がactivity prediction能力を示すわけでもない。「強い作用」の共通endpointが確保されていないため、その群集も評価不能。',
      '', '## D. ML rationale',
      '**indeterminate**。構造的パターンと一定のseed安定性はあるが、representation/metric依存と外部データ不足・context不一致が残る。追加薬剤のデータ整備を進める実行可能性は示せたが、supervised MLで薬理を予測できる根拠はまだ判断不能。今回のnegative/inconclusiveな対応結果を隠してcluster条件を変更していない。',
      '次は同一receptor composition/species/assay/endpointで外部yのdrug coverageを増やすことが先決。独立drug評価でdocking score/simple contactに対する追加表現の増分を検証する。十分なnの外部検証で増分が再現しなければ、fingerprintの予測情報価値仮説を支持しない。',
      '', '## Major limitations',
      '最大の制約は10剤に対する比較可能な外部薬理が3剤しかなく、しかもβ2/β3 contextが異なること。構造側にもα2 provisional local construct、単一受容体構造、単一初期conformer・固定ring、未指定立体3剤（lorazepam/temazepam/zopiclone）の単一sample、pH/microstate非列挙、Open Babel前処理の警告がある。複数poseは独立drugではない。受容体構造や立体/protonation変更に対するrobustnessは今回未評価。',
      '', '## QC / files',
      f"{qc['tests']['tests']} tests PASS、既存{qc['existing_files_unchanged']}ファイルのSHA256不変、{qc['plip_success']} XML保存/hash一致、archive全member一致、準備file固定、0/NA区別、非平均、外部annotation後も構造freeze一致。",
      'FINAL_REPORT/QC_REPORT、全CSV、完全environment lock、実効PLIP config、raw/archives、再現手順REPRODUCTION.mdを参照。PDF/PNGの読みやすい版はfigures/reviewed/。original frozen figuresはそのまま保持し、dendrogram/heatmap幅・PCA labelだけを修正したderivativeを追加（再計算・refitなし）。',
      '', '![Multi-receptor heatmap](../figures/reviewed/multireceptor_10drug_heatmap.png)',
      '', '![PCA](../figures/reviewed/fingerprint_pca.png)',
      '', '![External overview](../figures/structure_pharmacology_phenotype_overview.png)']
    (run/'reports/FINAL_REPORT.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
