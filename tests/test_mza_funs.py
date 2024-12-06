import numpy as np
import pandas as pd
import pytest
from mzsql import *

# Tests for MZA file type
def test_get_chrom_mza():
    ref_data = get_chrom_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1_idx.mzML", 118.0865, 10)
    test_data = get_chrom_mza("../demo_data/180205_Poo_TruePoo_Full1.mza", 118.0865, 10)

    # Merge data on 'rt' column
    merged_data = pd.merge(ref_data, test_data, on="rt", suffixes=("_ref", "_test"))

    assert (merged_data["mz_ref"] == merged_data["mz_test"]).all()
    assert (merged_data["int_ref"] == merged_data["int_test"]).all()


def test_get_spec_mza():
    spec_data_mzml = get_spec_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1_idx.mzML", 1)
    spec_data_mza = get_spec_mza("../demo_data/180205_Poo_TruePoo_Full1.mza", 1)

    # Align the data for comparison
    spec_data_mza, spec_data_mzml = spec_data_mza.align(spec_data_mzml, join='inner')

    assert (spec_data_mzml["mz"] == spec_data_mza["mz"]).all()
    assert (spec_data_mzml["int"] == spec_data_mza["int"]).all()

def test_get_rtrange_mza():
    rtrange_data_mzml = get_rtrange_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1_idx.mzML", 6.5, 8)
    rtrange_data_mza = get_rtrange_mza("../demo_data/180205_Poo_TruePoo_Full1.mza", 6.5, 8)

    # Align the data for comparison
    rtrange_data_mza, rtrange_data_mzml = rtrange_data_mza.align(rtrange_data_mzml, join='inner')

    assert (rtrange_data_mzml["rt"] == rtrange_data_mza["rt"]).all()
    assert (rtrange_data_mzml["mz"] == rtrange_data_mza["mz"]).all()
    assert (rtrange_data_mzml["int"] == rtrange_data_mza["int"]).all()

