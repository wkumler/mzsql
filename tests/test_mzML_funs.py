import numpy as np
import pandas as pd
from pyteomics import mzml, mzmlb
import pyopenms
import pytest
import pymzml
from mzsql import *

def test_gcm_pyteomics():
    chrom_data = get_chrom_mzml_pyteomics("../demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    assert chrom_data.shape == (1359, 3)
    assert min(chrom_data["mz"]) >= pmppm(118.0865, 10)[0]
    assert max(chrom_data["mz"]) <= pmppm(118.0865, 10)[1]
    assert max(chrom_data["int"]) < 4e8
    assert min(chrom_data["rt"] > 0)
    assert max(chrom_data["rt"] < 25)

def test_gcm_pymzml():
    chrom_data = get_chrom_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    assert chrom_data.shape == (1359, 3)
    assert min(chrom_data["mz"]) >= pmppm(118.0865, 10)[0]
    assert max(chrom_data["mz"]) <= pmppm(118.0865, 10)[1]
    assert max(chrom_data["int"]) < 4e8
    assert min(chrom_data["rt"] > 0)
    assert max(chrom_data["rt"] < 25)

def test_gcm_pyopenms():
    chrom_data = get_chrom_mzml_pyopenms("../demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    assert chrom_data.shape == (1359, 3)
    assert min(chrom_data["mz"]) >= pmppm(118.0865, 10)[0]
    assert max(chrom_data["mz"]) <= pmppm(118.0865, 10)[1]
    assert max(chrom_data["int"]) < 4e8
    assert min(chrom_data["rt"] > 0)
    assert max(chrom_data["rt"] < 25)

def test_gcm_pyopenms_2d():
    chrom_data = get_chrom_mzml_pyopenms_2DPeak("../demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    assert chrom_data.shape == (1359, 3)
    assert min(chrom_data["mz"]) >= pmppm(118.0865, 10)[0]
    assert max(chrom_data["mz"]) <= pmppm(118.0865, 10)[1]
    assert max(chrom_data["int"]) < 4e8
    assert min(chrom_data["rt"] > 0)
    assert max(chrom_data["rt"] < 25)

def test_intercompare_gcm_data():
    pyopenms_data = get_chrom_mzml_pyopenms("../demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    pyopenms_data_2d = get_chrom_mzml_pyopenms_2DPeak("../demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    pymzml_data = get_chrom_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    pyteomics_data = get_chrom_mzml_pyteomics("../demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    
    assert (pyopenms_data["rt"] == pyopenms_data_2d["rt"]).all()
    assert (pyopenms_data["rt"] == pymzml_data["rt"]).all()
    assert (pyopenms_data["rt"] == pyteomics["rt"]).all()
    assert (pymzml_data["rt"] == pyteomics["rt"]).all()

    assert (pyopenms_data["mz"] == pyopenms_data_2d["mz"]).all()
    assert (pyopenms_data["mz"] == pymzml_data["mz"]).all()
    assert (pyopenms_data["mz"] == pyteomics["mz"]).all()
    assert (pymzml_data["mz"] == pyteomics["mz"]).all()

    assert (pyopenms_data["int"] == pyopenms_data_2d["int"]).all()
    assert (pyopenms_data["int"] == pymzml_data["int"]).all()
    assert (pyopenms_data["int"] == pyteomics["int"]).all()
    assert (pymzml_data["int"] == pyteomics["int"]).all()



def test_gsm_pyteomics():
    spec_data = get_spec_mzml_pyteomics("../demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    assert spec_data.shape == (32, 2)
    assert min(spec_data["mz"]) >= 60
    assert max(spec_data["mz"]) <= 900
    assert max(spec_data["int"]) < 3e3

def test_gsm_pymzml():
    spec_data = get_spec_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1.mzML", 3)
    assert spec_data.shape == (32, 2)
    assert min(spec_data["mz"]) >= 60
    assert max(spec_data["mz"]) <= 900
    assert max(spec_data["int"]) < 3e3

def test_gsm_pyopenms():
    spec_data = get_spec_mzml_pyopenms("../demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    assert spec_data.shape == (32, 2)
    assert min(spec_data["mz"]) >= 60
    assert max(spec_data["mz"]) <= 900
    assert max(spec_data["int"]) < 3e3

def test_intercompare_gsm_data():
    pyopenms_data = get_spec_mzml_pyopenms("../demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    pymzml_data = get_spec_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1.mzML", 3)
    pyteomics_data = get_spec_mzml_pyteomics("../demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    
    assert (pyopenms_data["mz"] == pymzml_data["mz"]).all()
    assert (pyopenms_data["mz"] == pyteomics["mz"]).all()
    assert (pymzml_data["mz"] == pyteomics["mz"]).all()
    assert (pyopenms_data["int"] == pymzml_data["int"]).all()
    assert (pyopenms_data["int"] == pyteomics["int"]).all()
    assert (pymzml_data["int"] == pyteomics["int"]).all()



def test_grtm_pyteomics():
    rtrange_data = get_rtrange_mzml_pyteomics("../demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    assert rtrange_data.shape == (72636, 3)
    assert min(rtrange_data["rt"]) >= 6.5
    assert max(rtrange_data["rt"]) <= 8
    assert min(rtrange_data["mz"]) >= 60
    assert max(rtrange_data["mz"]) <= 900
    assert min(rtrange_data["int"]) >= 0

def test_grtm_pymzml():
    rtrange_data = get_rtrange_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    assert rtrange_data.shape == (72636, 3)
    assert min(rtrange_data["rt"]) >= 6.5
    assert max(rtrange_data["rt"]) <= 8
    assert min(rtrange_data["mz"]) >= 60
    assert max(rtrange_data["mz"]) <= 900
    assert min(rtrange_data["int"]) >= 0

def test_grtm_pyopenms():
    rtrange_data = get_rtrange_mzml_pyopenms("../demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    assert rtrange_data.shape == (72636, 3)
    assert min(rtrange_data["rt"]) >= 6.5
    assert max(rtrange_data["rt"]) <= 8
    assert min(rtrange_data["mz"]) >= 60
    assert max(rtrange_data["mz"]) <= 900
    assert min(rtrange_data["int"]) >= 0

def test_grtm_pyopenms():
    rtrange_data = get_rtrange_mzml_pyopenms_2DPeak("../demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    assert rtrange_data.shape == (72636, 3)
    assert min(rtrange_data["rt"]) >= 6.5
    assert max(rtrange_data["rt"]) <= 8
    assert min(rtrange_data["mz"]) >= 60
    assert max(rtrange_data["mz"]) <= 900
    assert min(rtrange_data["int"]) >= 0

def test_intercompare_grtm_data():
    pyopenms_data = get_rtrange_mzml_pyopenms("../demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    pyopenms_data_2d = get_rtrange_mzml_pyopenms_2DPeak("../demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    pymzml_data = get_rtrange_mzml_pymzml("../demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    pyteomics_data = get_rtrange_mzml_pyteomics("../demo_data/180205_Poo_TruePoo_Full1.mzML", 6.5, 8)
    
    assert (pyopenms_data["rt"] == pyopenms_data_2d["rt"]).all()
    assert (pyopenms_data["rt"] == pymzml_data["rt"]).all()
    assert (pyopenms_data["rt"] == pyteomics["rt"]).all()
    assert (pymzml_data["rt"] == pyteomics["rt"]).all()

    assert (pyopenms_data["mz"] == pyopenms_data_2d["mz"]).all()
    assert (pyopenms_data["mz"] == pymzml_data["mz"]).all()
    assert (pyopenms_data["mz"] == pyteomics["mz"]).all()
    assert (pymzml_data["mz"] == pyteomics["mz"]).all()

    assert (pyopenms_data["int"] == pyopenms_data_2d["int"]).all()
    assert (pyopenms_data["int"] == pymzml_data["int"]).all()
    assert (pyopenms_data["int"] == pyteomics["int"]).all()
    assert (pymzml_data["int"] == pyteomics["int"]).all()

