import numpy as np
import pandas as pd
from pyteomics import mzml, mzmlb
import pyopenms
import pytest
import pymzml
from mzsql import *

# Tests for DuckDB file type
def test_get_chrom_duckdb():
    ref_data = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    test_data = get_chrom_duckdb("demo_data/180205_Poo_TruePoo_Full1.duckdb", 118.0865, 10)
    assert ref_data["rt"].equals(test_data["rt"])
    assert ref_data["mz"].equals(test_data["mz"])
    assert ref_data["int"].equals(test_data["int"])

def test_get_spec_duckdb():
    spec_data_mzml = get_spec_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    spec_data_duckdb = get_spec_duckdb("demo_data/180205_Poo_TruePoo_Full1.duckdb", 1)
    assert spec_data_mzml["mz"] == spec_data_duckdb["mz"]
    assert spec_data_mzml["int"] == spec_data_duckdb["int"]

def test_get_rtrange_duckdb():
    rtrange_data_mzml = get_rtrange_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    rtrange_data_duckdb = get_rtrange_duckdb("demo_data/180205_Poo_TruePoo_Full1.duckdb", 1)
    assert rtrange_data_mzml["mz"] == rtrange_data_duckdb["mz"]
    assert rtrange_data_mzml["int"] == rtrange_data_duckdb["int"]