from pathlib import Path
import unittest,json,hashlib,tarfile,re,statistics,xml.etree.ElementTree as ET
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1];ROOT=R.parents[1];D=['diazepam','triazolam','brotizolam','lormetazepam']
class ModelCaseChecks(unittest.TestCase):
 def test_all_archived_pose_values_and_summaries(self):
  raw=pd.read_csv(R/'data/nishino_primary4_vina_all_poses.csv');best=pd.read_csv(R/'data/nishino_primary4_vina_seed_best.csv');summary=pd.read_csv(R/'data/nishino_primary4_vina_summary.csv')
  self.assertEqual(len(raw),354);self.assertEqual(len(best),40);self.assertEqual(len(summary),8)
  for f in pd.read_csv(R/'qc/docking_file_manifest.csv').itertuples():
   arc,member=f.source_pdbqt.split('::')
   with tarfile.open(R/arc) as tf:b=tf.extractfile(member).read()
   self.assertEqual(hashlib.sha256(b).hexdigest(),f.sha256)
   vals=[list(map(float,l.split(':')[1].split())) for l in b.decode().splitlines() if l.startswith('REMARK VINA RESULT:')]
   q=raw[(raw.drug==f.drug)&(raw.receptor_block==f.receptor)&(raw.seed==f.seed)].sort_values('pose_rank')
   np.testing.assert_array_equal(q[['vina_score_kcal_mol','rmsd_lb','rmsd_ub']],vals)
   v=best[(best.drug==f.drug)&(best.receptor_block==f.receptor)&(best.seed==f.seed)].iloc[0];self.assertEqual(v.vina_score_kcal_mol,min(a[0] for a in vals));self.assertEqual(v.pose_rank,1)
  for q in summary.itertuples():
   vals=best[(best.drug==q.drug)&(best.receptor_block==q.receptor_block)].vina_score_kcal_mol.tolist()
   self.assertEqual(q.n_seeds,5);self.assertAlmostEqual(q.median_best_seed_score,statistics.median(vals));self.assertAlmostEqual(q.sd_best_seed_score,statistics.stdev(vals));self.assertAlmostEqual(q.mean_best_seed_score,statistics.mean(vals));self.assertAlmostEqual(q.IQR_best_seed_score,sorted(vals)[3]-sorted(vals)[1])
 def test_plif_from_original_xml(self):
  x=pd.read_csv(R/'data/nishino_primary4_plif.csv').set_index('drug');mapping=pd.read_csv(R/'tables/residue_mapping_qc.csv');ix={(r.receptor,r.chain,int(r.residue_number)):(r.common_position,r.expected_residue) for r in mapping.itertuples() if r.mapping_valid}
  for d in D[2:]:
   hits={c:set() for c in x};den={b:set() for b in ['alpha1','alpha2']}
   with tarfile.open(R/'raw/archives'/f'{d}.tar.gz') as tf:
    for f in tf.getmembers():
     if not f.name.endswith('/report.xml'):continue
     drug,block,seed,pose,_=f.name.split('/');pid=(seed,pose);doc=ET.fromstring(tf.extractfile(f).read());sites=[s for s in doc.findall('bindingsite') if s.findtext('identifiers/hetid')=='LIG' and s.findtext('identifiers/chain')=='Z'];self.assertEqual(len(sites),1);den[block].add(pid)
     for node in sites[0].findall('interactions/*/*'):
      key=(block,node.findtext('reschain'),int(node.findtext('resnr')))
      if key in ix:
       cp,res=ix[key];col=f'{block}|{cp}|{res}|any_contact_frequency'
       if col in hits:hits[col].add(pid)
   for c in x:self.assertAlmostEqual(x.loc[d,c],len(hits[c])/len(den[c.split('|')[0]]),places=12)
 def test_frozen_two_drugs_and_hashes(self):
  x=pd.read_csv(R/'data/nishino_primary4_plif.csv').set_index('drug');old=pd.read_csv(R/'data/frozen_10drug_plif_source.csv').set_index('drug_id');np.testing.assert_array_equal(x.loc[D[:2]],old.loc[D[:2]])
  for p,h in json.loads((R/'config/STRUCTURAL_RESULTS_FROZEN.json').read_text())['files'].items():self.assertEqual(hashlib.sha256((R/p).read_bytes()).hexdigest(),h)
 def test_clinical_counts_mask_and_frozen_reuse(self):
  q=pd.read_csv(R/'qc/nishino_primary4_clinical_coverage.csv');self.assertEqual(len(q),16);self.assertEqual(int(q.estimable.sum()),14)
  for r in q.itertuples():
   a,b,c,d=r.a_target_event,r.b_target_nonevent,r.c_comparator_event,r.d_comparator_nonevent;self.assertEqual(a+b+c+d,r.total_reports)
   cells=np.array([a,b,c,d],float)
   if (cells==0).any():cells+=.5
   v=np.log(cells[0]*cells[3]/(cells[1]*cells[2]));self.assertAlmostEqual(v,r.raw_logROR_before_low_count_mask,places=12)
   self.assertEqual(pd.isna(r.logROR),a<5)
 def test_correspondence_independent_rank_pearson(self):
  x=pd.read_csv(R/'data/nishino_primary4_plif.csv').set_index('drug').loc[D];m=pd.read_csv(R/'data/nishino2008_primary4_in_vivo_M.csv').set_index('drug').loc[D];s=pd.read_csv(R/'data/nishino_primary4_structural_multiview.csv').set_index('drug').loc[D]
  for r in pd.read_csv(R/'data/nishino_primary4_X_M_correspondence.csv').itertuples():
   a=x.iloc[:,int(r.structural_feature[1:])-1] if r.structural_feature.startswith('F') else s[r.structural_feature]
   if a.nunique()<2:self.assertTrue(pd.isna(r.spearman_rho))
   else:self.assertAlmostEqual(r.spearman_rho,a.rank().corr(m.rotarod_potency.rank()),places=12)
   self.assertEqual(r.n_effective,4)
  my=pd.read_csv(R/'data/nishino_primary4_M_Y_correspondence.csv');self.assertEqual(my.n_effective.tolist(),[4,4,2,4]);self.assertTrue(pd.isna(my.iloc[2].spearman_rho))
 def test_raw_in_vivo_and_identity_alignment(self):
  raw=pd.read_csv(R/'data/nishino2008_rotarod_raw.csv');self.assertEqual(len(raw),80);self.assertEqual(set(raw.drug),set(D+['rilmazafone']));self.assertTrue((raw.n_total==10).all());np.testing.assert_allclose(raw.failure_fraction,raw.n_positive/raw.n_total)
  q=raw[(raw.drug=='triazolam')&(raw.dose_mg_kg==.5)&(raw.time_min==90)];self.assertEqual(q.n_positive.iloc[0],9)
  master=pd.read_csv(R/'data/nishino_primary4_XMY_master.csv');self.assertEqual(master.drug.tolist(),D);np.testing.assert_allclose(master.rotarod_potency,-np.log10(master.ed50_mg_kg));self.assertEqual(master.ed50_mg_kg.tolist(),[3.11,1.25,5.76,3.39])
if __name__=='__main__':unittest.main()
