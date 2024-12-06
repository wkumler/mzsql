import numpy as np
import pandas as pd
import pytest
from mzsql import *

# Tests for SQLite file type
def test_get_chrom_sqlite():
    ref_data = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    test_data = get_chrom_sqlite("demo_data/180205_Poo_TruePoo_Full1_ordered_rt.sqlite", 118.0865, 10)
    assert (ref_data["rt"] == test_data["rt"]).all()
    assert (ref_data["mz"] == test_data["mz"]).all()
    assert (ref_data["int"] == test_data["int"]).all()

def test_get_spec_sqlite():
    spec_data_mzml = get_spec_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    spec_data_sqlite = get_spec_sqlite("demo_data/180205_Poo_TruePoo_Full1_ordered_rt.sqlite", 1)
    assert (spec_data_mzml["mz"] == spec_data_sqlite["mz"]).all()
    assert (spec_data_mzml["int"] == spec_data_sqlite["int"]).all()

def test_get_rtrange_sqlite():
    rtrange_data_mzml = get_rtrange_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    rtrange_data_sqlite = get_rtrange_sqlite("demo_data/180205_Poo_TruePoo_Full1_ordered_rt.sqlite", 6.5, 8)
    assert (rtrange_data_mzml["rt"] == rtrange_data_sqlite["rt"]).all()
    assert (rtrange_data_mzml["mz"] == rtrange_data_sqlite["mz"]).all()
    assert (rtrange_data_mzml["int"] == rtrange_data_sqlite["int"]).all()