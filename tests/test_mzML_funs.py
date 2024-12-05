from mzsql import *

def test_gcm_pyteomics():
    chrom_data = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    assert chrom_data.shape == (1359, 3)
    assert min(chrom_data["mz"]) >= pmppm(118.0865, 10)[0]
    assert max(chrom_data["mz"]) <= pmppm(118.0865, 10)[1]
    assert max(chrom_data["int"]) < 4e8
    assert min(chrom_data["rt"] > 0)
    assert max(chrom_data["rt"] < 25)

def test_get_chrom_mzdb():
    chrom_data_mzml = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    chrom_data_mzdb = get_chrom_mzdb("demo_data/180205_Poo_TruePoo_Full1.raw.mzDB", 118.0865, 10)
    assert chrom_data_mzml["rt"] == chrom_data_mzdb["rt"]
    assert chrom_data_mzml["mz"] == chrom_data_mzdb["mz"]
    assert chrom_data_mzml["int"] == chrom_data_mzdb["int"]

def test_get_spec_mzdb():
    spec_data_mzml = get_spec_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 1)
    spec_data_mzdb = get_spec_mzdb("demo_data/180205_Poo_TruePoo_Full1.raw.mzDB", 1)
    assert spec_data_mzml["mz"] == spec_data_mzdb["mz"]
    assert spec_data_mzml["int"] == spec_data_mzdb["int"]
