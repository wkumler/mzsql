import numpy as np
import pandas as pd
import pytest
from mzsql import *




# Tests for DuckDB file type
def test_get_chrom_duckdb():
    ref_data = get_chrom_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    test_data = get_chrom_duckdb("../demo_data/180205_Poo_TruePoo_Full1.duckdb", 118.0865, 10)
    
    # Align indices before comparison
    ref_data, test_data = ref_data.align(test_data, join='inner')
    
    assert (ref_data["rt"] == test_data["rt"]).all()
    assert (ref_data["mz"] == test_data["mz"]).all()
    assert (ref_data["int"] == test_data["int"]).all()

def test_get_spec_duckdb():
    spec_data_mzml = get_spec_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1_idx.mzML", 1)
    spec_data_duckdb = get_spec_duckdb("../demo_data/180205_Poo_TruePoo_Full1.duckdb", 1)
    
    # Align indices before comparison
    spec_data_mzml, spec_data_duckdb = spec_data_mzml.align(spec_data_duckdb, join='inner')
    
    assert (spec_data_mzml["mz"] == spec_data_duckdb["mz"]).all()
    assert (spec_data_mzml["int"] == spec_data_duckdb["int"]).all()

def test_get_rtrange_duckdb():
    rtrange_data_mzml = get_rtrange_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    rtrange_data_duckdb = get_rtrange_duckdb("../demo_data/180205_Poo_TruePoo_Full1.duckdb", 6.5, 8)
    
    # Align indices before comparison
    rtrange_data_mzml, rtrange_data_duckdb = rtrange_data_mzml.align(rtrange_data_duckdb, join='inner')

    assert (rtrange_data_mzml["rt"] == rtrange_data_duckdb["rt"]).all()
    assert (rtrange_data_mzml["mz"] == rtrange_data_duckdb["mz"]).all()
    assert (rtrange_data_mzml["int"] == rtrange_data_duckdb["int"]).all()
