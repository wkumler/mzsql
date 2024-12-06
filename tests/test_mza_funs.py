import numpy as np
import pandas as pd
import pytest
from mzsql import *

# Tests for MZA file type
def test_get_chrom_mza():
    ref_data = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    test_data = get_chrom_mza("demo_data/180205_Poo_TruePoo_Full1.mza", 118.0865, 10)
    assert ref_data["rt"].equals(test_data["rt"])
    assert ref_data["mz"].equals(test_data["mz"])
    assert ref_data["int"].equals(test_data["int"])

def test_get_spec_mza():
    spec_data_mzml = get_spec_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    spec_data_mza = get_spec_mza("demo_data/180205_Poo_TruePoo_Full1.mza", 1)
    assert spec_data_mzml["mz"] == spec_data_mza["mz"]
    assert spec_data_mzml["int"] == spec_data_mza["int"]

def test_get_rtrange_mza():
    rtrange_data_mzml = get_rtrange_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    rtrange_data_mza = get_rtrange_mza("demo_data/180205_Poo_TruePoo_Full1.mza", 6.5, 8)
    assert rtrange_data_mzml["mz"] == rtrange_data_mza["mz"]
    assert rtrange_data_mzml["int"] == rtrange_data_mza["int"]
