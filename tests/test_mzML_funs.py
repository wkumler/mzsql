import numpy as np
import pandas as pd
from pyteomics import mzml, mzmlb
import pyopenms
import pymzml

def pmppm(mass, ppm=4):
    return(mass*(1-ppm/1000000), mass*(1+ppm/1000000))

def get_chrom_mzml_pyteomics(file, mz, ppm):
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    for spectrum in mzml.MzML(file):
        rt_val = spectrum['scanList']['scan'][0]['scan start time']
        mz_vals=spectrum['m/z array']
        int_vals = spectrum['intensity array']
        bet_idxs = (mzmin < spectrum["m/z array"]) & (spectrum["m/z array"] < mzmax)
        if(sum(bet_idxs)>0):
            df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
def test_gcm_pyteomics():
    chrom_data = get_chrom_mzml_pyteomics("demo_data/180205_Poo_TruePoo_Full1.mzML", 118.0865, 10)
    assert chrom_data.shape == (1359, 3)
    assert min(chrom_data["mz"]) >= pmppm(118.0865, 10)[0]
    assert max(chrom_data["mz"]) <= pmppm(118.0865, 10)[1]
    assert max(chrom_data["int"]) < 4e8
    assert min(chrom_data["rt"] > 0)
    assert max(chrom_data["rt"] < 25)





def get_chrom_mzml_pyopenms(file, mz, ppm):
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    for spectrum in exp:
        rt_val = spectrum.getRT()
        mz_vals, int_vals = spectrum.get_peaks()
        bet_idxs = mzmin < mz_vals < mzmax
        if(sum(bet_idxs)>0):
            df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
def get_chrom_mzml_pymzml(file, mz, ppm):
    run = pymzml.run.Reader(file)
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    for spectrum in run:
        rt_val = spectrum.scan_time_in_minutes()
        mz_vals = spectrum.mz
        int_vals = spectrum.i
        bet_idxs = mzmin < mz_vals < mzmax

        if(sum(bet_idxs)>0):
            df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
def get_chrom_mzml_pyopenms_2DPeak(file, mz, ppm):
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    exp.updateRanges()
    mzmin, mzmax = pmppm(mz, ppm)
    chrom_data=exp.get2DPeakDataLong(min_mz=mzmin, max_mz=mzmax, min_rt=exp.getMinRT(), max_rt=exp.getMaxRT())
    return(pd.DataFrame({"rt":chrom_data[0], "mz":chrom_data[1], "int":chrom_data[2]}))
