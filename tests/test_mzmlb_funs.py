import numpy as np
import pandas as pd
import pytest
from mzsql import *

# Tests for mzMLB file type
def test_get_chrom_mzmlb():
    ref_data = get_chrom_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    test_data = get_chrom_mzmlb("../demo_data/180205_Poo_TruePoo_Full1.mzMLb", 118.0865, 10)

    # Merge data on 'rt' column
    merged_data = pd.merge(ref_data, test_data, on="rt", suffixes=("_ref", "_test"))

    assert (merged_data["mz_ref"] == merged_data["mz_test"]).all()
    assert (merged_data["int_ref"] == merged_data["int_test"]).all()

#expected failure
def test_get_spec_mzmlb():
    spec_data_mzml = get_spec_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    spec_data_mzmlb = get_spec_mzmlb("../demo_data/180205_Poo_TruePoo_Full1.mzMLb", 1)

    # Merge data on 'rt' column
    spec_data_mzml, spec_data_mzmlb = spec_data_mzml.align(spec_data_mzmlb, join='inner')

    assert (spec_data_mzml["mz"] == spec_data_mzmlb["mz"]).all()
    assert (spec_data_mzml["int"] == spec_data_mzmlb["int"]).all()

def test_get_rtrange_mzmlb():
    rtrange_data_mzml = get_rtrange_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    rtrange_data_mzmlb = get_rtrange_mzmlb("../demo_data/180205_Poo_TruePoo_Full1.mzMLb", 6.5, 8)

    rtrange_data_mzml, rtrange_data_mzmlb = rtrange_data_mzml.align(rtrange_data_mzmlb, join='inner')

    assert (rtrange_data_mzml["rt"] == rtrange_data_mzmlb["rt"]).all()
    assert (rtrange_data_mzml["mz"] == rtrange_data_mzmlb["mz"]).all()
    assert (rtrange_data_mzml["int"] == rtrange_data_mzmlb["int"]).all()
