import numpy as np
import pandas as pd
import pytest
from mzsql import *

def test_get_chrom_mzdb():
    ref_data = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    test_data = get_chrom_mzdb("demo_data/180205_Poo_TruePoo_Full1.raw.mzDB", 118.0865, 10)
    assert np.allclose(ref_data["rt"], test_data["rt"])
    assert np.allclose(ref_data["mz"], test_data["mz"])
    assert np.allclose(ref_data["int"], test_data["int"])

def test_get_spec_mzdb():
    spec_data_mzml = get_spec_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    spec_data_mzdb = get_spec_mzdb("demo_data/180205_Poo_TruePoo_Full1.raw.mzDB", 1)
    assert np.allclose(spec_data_mzml["mz"], spec_data_mzdb["mz"])
    assert np.allclose(spec_data_mzml["int"], spec_data_mzdb["int"])

def test_get_rtrange_mzdb():
    rtrange_data_mzml = get_rtrange_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    rtrange_data_mzdb = get_rtrange_mzdb("demo_data/180205_Poo_TruePoo_Full1.raw.mzDB", 6.5, 8)
    assert np.allclose(rtrange_data_mzml["rt"], rtrange_data_mzdb["rt"])
    assert np.allclose(rtrange_data_mzml["mz"], rtrange_data_mzdb["mz"])
    assert np.allclose(rtrange_data_mzml["int"], rtrange_data_mzdb["int"])
