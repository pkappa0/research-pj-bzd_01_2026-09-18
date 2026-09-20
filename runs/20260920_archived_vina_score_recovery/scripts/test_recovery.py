import unittest, statistics, hashlib
from pathlib import Path
import pandas as pd
import numpy as np
from recover_vina_scores import parse_pdbqt
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
if __name__=='__main__':unittest.main()
