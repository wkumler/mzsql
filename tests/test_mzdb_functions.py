import numpy as np
import pandas as pd
import pytest
from mzsql import *

def test_get_chrom_mzdb():
    ref_data = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    test_data = get_chrom_mzdb("demo_data/180205_Poo_TruePoo_Full1.raw.mzDB", 118.0865, 10)
    assert (ref_data["rt"] == test_data["rt"]).all()
    assert (ref_data["mz"] == test_data["mz"]).all()
    assert (ref_data["int"] == test_data["int"]).all()

def test_get_spec_mzdb():
    spec_data_mzml = get_spec_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    spec_data_mzdb = get_spec_mzdb("demo_data/180205_Poo_TruePoo_Full1.raw.mzDB", 1)
    assert (spec_data_mzml["mz"] == spec_data_mzdb["mz"]).all()
    assert (spec_data_mzml["int"] == spec_data_mzdb["int"]).all()

def test_get_rtrange_mzdb():
    rtrange_data_mzml = get_rtrange_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    rtrange_data_mzdb = get_rtrange_mzdb("demo_data/180205_Poo_TruePoo_Full1.raw.mzDB", 6.5, 8)
    assert (rtrange_data_mzml["rt"] == rtrange_data_mzdb["rt"]).all()
    assert (rtrange_data_mzml["mz"] == rtrange_data_mzdb["mz"]).all()
    assert (rtrange_data_mzml["int"] == rtrange_data_mzdb["int"]).all()
