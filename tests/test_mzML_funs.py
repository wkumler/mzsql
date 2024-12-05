import numpy as np
import pandas as pd
from pyteomics import mzml, mzmlb
import pyopenms
import pymzml
from mzsql import *

def test_gcm_pyteomics():
    chrom_data = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    assert chrom_data.shape == (1359, 3)
    assert min(chrom_data["mz"]) >= pmppm(118.0865, 10)[0]
    assert max(chrom_data["mz"]) <= pmppm(118.0865, 10)[1]
    assert max(chrom_data["int"]) < 4e8
    assert min(chrom_data["rt"] > 0)
    assert max(chrom_data["rt"] < 25)
