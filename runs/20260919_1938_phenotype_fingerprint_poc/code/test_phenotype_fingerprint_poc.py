import numpy as np
import pandas as pd
from src.phenotype_fingerprint_poc import normalize, rho_result, KEYS

def test_missing_is_not_zero_and_duration_not_scaled():
    base=dict(raw_value=np.nan,sample_size=10,control_value=0,raw_unit='positive_mice_count',outcome_name='failure')
    assert np.isnan(normalize(base)[0])
    assert normalize(dict(base,raw_value=0))[0]==0
    assert normalize(dict(base,raw_value=5))[0]==.5
    assert normalize(dict(base,raw_value=205.3,raw_unit='seconds_mean'))[0]==205.3

def test_reduction_and_no_clipping():
    base=dict(raw_value=20,sample_size=6,control_value=100,raw_unit='counts',outcome_name='locomotor_activity')
    assert normalize(base)[0]==.8
    assert normalize(dict(base,raw_value=120))[0]==-.2
    assert np.isnan(normalize(dict(base,control_value=0))[0])

def test_effective_drugs_not_rows_and_two_point_rho():
    rho,n,status=rho_result(np.array([1,2]),np.array([2,1]))
    assert np.isnan(rho) and n==2 and status=='n_effective_lt_3'
    rho,n,status=rho_result(np.array([1.,2.,3.]),np.array([3.,2.,1.]))
    assert rho==-1 and n==3
    assert np.isnan(rho_result(np.ones(3),np.arange(3.))[0])

def test_comparability_includes_confounding_conditions():
    for key in ['dose','vehicle','co_treatment','acute_chronic','species','strain','sex','observation_time','administration_route']:
        assert key in KEYS
