import numpy as np
import pandas as pd
import pytest
from mzsql import *

# Tests for mzMLB file type
def test_get_chrom_mzmlb():
    ref_data = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    test_data = get_chrom_mzmlb("demo_data/180205_Poo_TruePoo_Full1.mzMLb", 118.0865, 10)
    assert ref_data["rt"].equals(test_data["rt"])
    assert ref_data["mz"].equals(test_data["mz"])
    assert ref_data["int"].equals(test_data["int"])

def test_get_spec_mzmlb():
    spec_data_mzml = get_spec_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    spec_data_mzmlb = get_spec_mzmlb("demo_data/180205_Poo_TruePoo_Full1.mzMLb", 1)
    assert spec_data_mzml["mz"] == spec_data_mzmlb["mz"]
    assert spec_data_mzml["int"] == spec_data_mzmlb["int"]

def test_get_rtrange_mzmlb():
    rtrange_data_mzml = get_rtrange_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    rtrange_data_mzmlb = get_rtrange_mzmlb("demo_data/180205_Poo_TruePoo_Full1.mzMLb", 6.5, 8)
    assert rtrange_data_mzml["mz"] == rtrange_data_mzmlb["mz"]
    assert rtrange_data_mzml["int"] == rtrange_data_mzmlb["int"]
