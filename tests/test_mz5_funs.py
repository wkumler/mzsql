import numpy as np
import pandas as pd
import pytest
from mzsql import *

# Tests for mz5 file type
def test_get_chrom_mz5():
    ref_data = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    test_data = get_chrom_mz5("demo_data/180205_Poo_TruePoo_Full1.mz5", 118.0865, 10)
    assert (ref_data["rt"] == test_data["rt"]).all()
    assert (ref_data["mz"] == test_data["mz"]).all()
    assert (ref_data["int"] == test_data["int"]).all()

def test_get_spec_mz5():
    spec_data_mzml = get_spec_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    spec_data_mz5 = get_spec_mz5("demo_data/180205_Poo_TruePoo_Full1.mz5", 1)
    assert (spec_data_mzml["mz"] == spec_data_mz5["mz"]).all()
    assert (spec_data_mzml["int"] == spec_data_mz5["int"]).all()

def test_get_rtrange_mz5():
    rtrange_data_mzml = get_rtrange_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    rtrange_data_mz5 = get_rtrange_mz5("demo_data/180205_Poo_TruePoo_Full1.mz5", 6.5, 8)
    assert (rtrange_data_mzml["rt"] == rtrange_data_mz5["rt"]).all()
    assert (rtrange_data_mzml["mz"] == rtrange_data_mz5["mz"]).all()
    assert (rtrange_data_mzml["int"] == rtrange_data_mz5["int"]).all()