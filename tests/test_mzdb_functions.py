import numpy as np
import pandas as pd
from pyteomics import mzml, mzmlb
import pyopenms
import pymzml
from mzsql import *

def test_get_chrom_mzdb():
    ref_data = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    test_data = get_chrom_mzdb("demo_data/180205_Poo_TruePoo_Full1.raw.mzDB", 118.0865, 10)
    test_data["rt"] /= 60  # Adjust rt for mzdb
    assert ref_data["rt"].equals(test_data["rt"])
    assert ref_data["mz"].equals(test_data["mz"])
    assert ref_data["int"].equals(test_data["int"])

def test_get_spec_mzdb():
    spec_data_mzml = get_spec_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    spec_data_mzdb = get_spec_mzdb("demo_data/180205_Poo_TruePoo_Full1.raw.mzDB", 1)
    assert spec_data_mzml["mz"] == spec_data_mzdb["mz"]
    assert spec_data_mzml["int"] == spec_data_mzdb["int"]

def test_get_rtrange_mzdb():
    rtrange_data_mzml = get_rtrange_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    rtrange_data_mzdb = get_rtrange_mzdb("demo_data/180205_Poo_TruePoo_Full1.raw.mzDB", 1)
    assert rtrange_data_mzml["mz"] == rtrange_data_mzdb["mz"]
    assert rtrange_data_mzml["int"] == rtrange_data_mzdb["int"]