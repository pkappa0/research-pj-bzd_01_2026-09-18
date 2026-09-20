import unittest, statistics, hashlib
from pathlib import Path
import pandas as pd
import numpy as np
from recover_vina_score_layer import parse_pdbqt
R=Path(__file__).resolve().parents[1]
class ParserChecks(unittest.TestCase):
 def parse(self,s):return parse_pdbqt(s,'drug','alpha1',1,'fixture')
 def test_preserve_nonbest(self):
  r,e,n=self.parse('MODEL 1\nREMARK VINA RESULT: -10.1 0 0\nENDMDL\nMODEL 2\nREMARK VINA RESULT: -8.0 1.2 2.3\nENDMDL\n')
  self.assertEqual(len(r),2);self.assertEqual(n,2);self.assertEqual(e,[]);self.assertEqual(r[1]['vina_score_kcal_mol'],-8)
 def test_malformed_retained(self):
  r,e,n=self.parse('MODEL 1\nREMARK VINA RESULT: broken 0 0\nENDMDL\n');self.assertEqual(len(r),1);self.assertTrue(np.isnan(r[0]['vina_score_kcal_mol']));self.assertIn('malformed_vina_result',[i['issue'] for i in e])
 def test_missing_model_flagged(self):
  r,e,n=self.parse('REMARK VINA RESULT: -7.0 0 0\n');self.assertIsNone(r[0]['pose_rank']);self.assertIn('result_without_model_rank',[i['issue'] for i in e])
 def test_invalid_rmsd_flagged(self):
  r,e,n=self.parse('MODEL 1\nREMARK VINA RESULT: -7.0 2 1\nENDMDL\n');self.assertEqual(len(r),1);self.assertIn('malformed_or_invalid_numeric_record',[i['issue'] for i in e])
 def test_missing_result_flagged(self):
  r,e,n=self.parse('MODEL 1\nATOM dummy\nENDMDL\n');self.assertEqual(r,[]);self.assertIn('model_result_count',[i['issue'] for i in e])
 def test_duplicate_result_preserved_flagged(self):
  r,e,n=self.parse('MODEL 1\nREMARK VINA RESULT: -7.0 0 0\nREMARK VINA RESULT: -6.9 0 0\nENDMDL\n');self.assertEqual(len(r),2);self.assertIn('model_result_count',[i['issue'] for i in e])
 def test_saved_results_independently(self):
  raw=pd.read_csv(R/'data/vina_scores_all_poses.csv');best=pd.read_csv(R/'data/vina_scores_seed_best.csv');summ=pd.read_csv(R/'data/vina_scores_drug_receptor_summary.csv');fm=pd.read_csv(R/'data/pdbqt_file_manifest.csv')
  self.assertEqual(len(raw),883);self.assertEqual(len(best),100);self.assertEqual(len(summ),20)
  self.assertFalse(raw.duplicated(['drug','receptor_block','seed','pose_rank']).any())
  for f in fm.itertuples():
   b=(R/f.snapshot_file).read_bytes();self.assertEqual(hashlib.sha256(b).hexdigest(),f.sha256)
   lines=[line for line in b.decode().splitlines() if line.startswith('REMARK VINA RESULT:')]
   cells=raw[raw.source_file==f.source_file].sort_values('record_index');self.assertEqual(len(cells),len(lines));self.assertEqual(len(lines),f.expected_n_poses)
   for line,row in zip(lines,cells.itertuples()):
    fields=line.split(':',1)[1].split();self.assertEqual(float(fields[0]),row.vina_score_kcal_mol);self.assertEqual(float(fields[1]),row.rmsd_lb);self.assertEqual(float(fields[2]),row.rmsd_ub)
  for row in best.itertuples():
   v=raw[(raw.drug==row.drug)&(raw.receptor_block==row.receptor_block)&(raw.seed==row.seed)]
   self.assertEqual(row.vina_score_kcal_mol,min(v.vina_score_kcal_mol));self.assertEqual(row.pose_rank,1)
  for row in summ.itertuples():
   vals=sorted(best[(best.drug==row.drug)&(best.receptor_block==row.receptor_block)].vina_score_kcal_mol)
   self.assertEqual(len(vals),5);self.assertEqual(row.n_seeds_available,5)
   self.assertAlmostEqual(row.median_best_seed_score,statistics.median(vals),places=12);self.assertAlmostEqual(row.mean_best_seed_score,statistics.mean(vals),places=12)
   self.assertEqual(row.min_best_seed_score,vals[0]);self.assertEqual(row.max_best_seed_score,vals[-1]);self.assertAlmostEqual(row.IQR_best_seed_score,vals[3]-vals[1],places=12)
class LayerChecks(unittest.TestCase):
 def test_sd_and_favorability(self):
  best=pd.read_csv(R/'data/vina_scores_seed_best.csv')
  for row in pd.read_csv(R/'data/vina_scores_drug_receptor_summary.csv').itertuples():
   v=best[(best.drug==row.drug)&(best.receptor_block==row.receptor_block)].best_vina_score_kcal_mol
   self.assertAlmostEqual(row.sd_best_seed_score,statistics.stdev(v),places=12)
   self.assertEqual(row.vina_favorability,-row.median_best_seed_score)
 def test_correlations_from_independent_ranks(self):
  y=pd.read_csv(R/'data/frozen_clinical_fingerprint.csv').set_index('drug')
  scores=pd.read_csv(R/'data/vina_scores_drug_receptor_summary.csv')
  for row in pd.read_csv(R/'data/vina_clinical_c1_correspondence.csv').itertuples():
   x=scores[scores.receptor_block==row.receptor_block].set_index('drug').median_best_seed_score
   pair=pd.concat([x.rename('x'),y[row.clinical_term].rename('y')],axis=1).dropna()
   rho=pair.x.rank().corr(pair.y.rank())
   self.assertEqual(row.n_effective,len(pair))
   self.assertAlmostEqual(row.rho_raw_vina_score,rho,places=12)
   self.assertAlmostEqual(row.rho_vina_favorability,-rho,places=12)
 def test_frozen_hashes_and_vector_references(self):
  import json
  frozen=json.loads((R/'VINA_LAYER_FROZEN.json').read_text())
  for f,h in frozen['files'].items():self.assertEqual(hashlib.sha256((R/'data'/f).read_bytes()).hexdigest(),h)
  refs=json.loads((R/'data/plif_vector_reference_manifest.json').read_text())
  self.assertEqual(len(refs),20)
  x=pd.read_csv(R/'data/frozen_plif_fingerprint.csv').set_index('drug_id')
  idx=pd.read_csv(R/'data/structural_multiview_index.csv');self.assertEqual(len(idx),10)
  for row in idx.to_dict('records'):
   for rec,n in [('alpha1',21),('alpha2',27)]:
    ref=refs[row[rec+'_plif_vector_identifier']]
    self.assertEqual(ref['row_key'],row['drug'])
    self.assertEqual(len(ref['feature_names_in_frozen_order']),n)
    self.assertEqual(ref['feature_names_in_frozen_order'],[c for c in x if c.startswith(rec+'|')])
    self.assertEqual(hashlib.sha256((R/ref['source_file']).read_bytes()).hexdigest(),ref['source_sha256'])
if __name__=='__main__':unittest.main()
