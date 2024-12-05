
import numpy as np
import pandas as pd
import pyteomics.mzml
import pyopenms
import pymzml
from .helpers import pmppm

# pyteomics things
def get_chrom_mzml_pyteomics(file, mz, ppm):
    """
    Retrieves chromatogram data from a mzML file using pyteomics for a specific m/z range.

    Args:
        file (str): Path to the mzML file.
        mz (float): The target m/z value.
        ppm (float): The mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified m/z range.
    """
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    for spectrum in pyteomics.mzml.MzML(file):
        rt_val = spectrum['scanList']['scan'][0]['scan start time']
        mz_vals=spectrum['m/z array']
        int_vals = spectrum['intensity array']
        bet_idxs = (mzmin < spectrum["m/z array"]) & (spectrum["m/z array"] < mzmax)
        if(sum(bet_idxs)>0):
            df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
            scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))
    
def get_spec_mzml_pyteomics(file, scan_num):
    """
    Retrieves spectrum data for a specific scan number from a mzML file using pyteomics.

    Args:
        file (str): Path to the mzML file.
        scan_num (int): The scan number.

    Returns:
        pd.DataFrame: A DataFrame containing m/z and intensity data for the specified scan number.
    """
    file_data = pyteomics.mzml.MzML(file)
    return(pd.DataFrame({"mz":file_data[scan_num]['m/z array'], "int":file_data[scan_num]['intensity array']}))


def get_rtrange_mzml_pyteomics(file, rtstart, rtend):
    """
    Retrieves chromatogram data for a specific retention time range from a mzML file using pyteomics.

    Args:
        file (str): Path to the mzML file.
        rtstart (float): The start retention time (in minutes).
        rtend (float): The end retention time (in minutes).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified retention time range.
    """
    scan_dfs = []
    for spectrum in pyteomics.mzml.MzML(file):
        rt_val = spectrum['scanList']['scan'][0]['scan start time']
        if(rtstart < rt_val < rtend):
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']
            df_scan = pd.DataFrame({'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(mz_vals)})
            scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))



# pyopenms things
def get_chrom_mzml_pyopenms(file, mz, ppm):
    """
    Retrieves chromatogram data from a mzML file using pyopenms for a specific m/z range.

    Args:
        file (str): Path to the mzML file.
        mz (float): The target m/z value.
        ppm (float): The mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified m/z range.
    """
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    for spectrum in exp:
        rt_val = spectrum.getRT()
        mz_vals, int_vals = spectrum.get_peaks()
        bet_idxs = (mzmin < mz_vals) & (mz_vals < mzmax)
        if(sum(bet_idxs)>0):
            df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
            scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))
    
def get_chrom_mzml_pyopenms_2DPeak(file, mz, ppm):
    """
    Retrieves 2D chromatogram data from a mzML file using pyopenms for a specific m/z range.

    Args:
        file (str): Path to the mzML file.
        mz (float): The target m/z value.
        ppm (float): The mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified m/z range in a 2D peak format.
    """
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    exp.updateRanges()
    mzmin, mzmax = pmppm(mz, ppm)
    chrom_data=exp.get2DPeakDataLong(min_mz=mzmin, max_mz=mzmax, min_rt=exp.getMinRT(), max_rt=exp.getMaxRT())
    return(pd.DataFrame({"rt":chrom_data[0], "mz":chrom_data[1], "int":chrom_data[2]}))
    
def get_spec_mzml_pyopenms(file, scan_num):
    """
    Retrieves spectrum data for a specific scan number from a mzML file using pyopenms.

    Args:
        file (str): Path to the mzML file.
        scan_num (int): The scan number.

    Returns:
        pd.DataFrame: A DataFrame containing m/z and intensity data for the specified scan number.
    """
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    spec1_data = exp[scan_num].get_peaks()
    return(pd.DataFrame({"mz":spec1_data[0], "int":spec1_data[1]}))
    
def get_rtrange_mzml_pyopenms(file, rtstart, rtend):
    """
    Retrieves chromatogram data for a specific retention time range from a mzML file using pyopenms.

    Args:
        file (str): Path to the mzML file.
        rtstart (float): The start retention time (in minutes).
        rtend (float): The end retention time (in minutes).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified retention time range.
    """
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    exp.updateRanges()
    scan_dfs = []
    for spectrum in exp:
        rt_val = spectrum.getRT()
        if(rtstart*60 < rt_val < rtend*60):
            mz_vals, int_vals = spectrum.get_peaks()
            df_scan = pd.DataFrame({'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(int_vals)})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
    
def get_rtrange_mzml_pyopenms_2DPeak(file, rtstart, rtend):
    """
    Retrieves 2D chromatogram data from a mzML file using pyopenms for a specific retention time range.

    Args:
        file (str): Path to the mzML file.
        rtstart (float): The start retention time (in minutes).
        rtend (float): The end retention time (in minutes).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified retention time range in a 2D peak format.
    """
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    exp.updateRanges()
    rtrange_data=exp.get2DPeakDataLong(min_mz=exp.getMinMZ(), max_mz=exp.getMaxMZ(), min_rt=rtstart*60, max_rt=rtend*60)
    return(pd.DataFrame({"rt":rtrange_data[0], "mz":rtrange_data[1], "int":rtrange_data[2]}))




# pymzml things
# Example available at https://github.com/pymzml/pymzML/blob/dev/example_scripts/extract_ion_chromatogram.py
# but does not show how to add a ppm tolerance
def get_chrom_mzml_pymzml(file, mz, ppm):
     """
    Retrieves chromatogram data from a mzML file using pymzml for a specific m/z range.

    Args:
        file (str): Path to the mzML file.
        mz (float): The target m/z value.
        ppm (float): The mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified m/z range.
    """
    run = pymzml.run.Reader(file, build_index_from_scratch=True)
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    for spectrum in run:
        rt_val = spectrum.scan_time_in_minutes()
        mz_vals = spectrum.mz
        int_vals = spectrum.i
        bet_idxs = (mzmin < mz_vals) & (mz_vals < mzmax)
    
        if(sum(bet_idxs)>0):
            df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))

def get_spec_mzml_pymzml(file, scan_num):
    # fails for the first scan when given a non-indexed mzML even if build_index is True
    # run = pymzml.run.Reader("demo_data/180205_Poo_TruePoo_Full1.mzML")
    # run = pymzml.run.Reader("demo_data/180205_Poo_TruePoo_Full1.mzML", build_index_from_scratch=True)
    # run = pymzml.run.Reader("demo_data/180205_Poo_TruePoo_Full1_idx.mzML")
    # pymzml seems to reference their spectra by scan number
    # spec1_data = run[1].peaks("raw")
    # print(spec1_data)
    # plt.stem(spec1_data[:,0], spec1_data[:,1])
    run = pymzml.run.Reader(file, build_index_from_scratch=True)
    spec1_data = run[scan_num].peaks("raw")
    return(pd.DataFrame({"mz":spec1_data[:,0], "int":spec1_data[:,1]}))

def get_rtrange_mzml_pymzml(file, rtstart, rtend):
    """
    Retrieves chromatogram data for a specific retention time range from a mzML file using pymzml.

    Args:
        file (str): Path to the mzML file.
        rtstart (float): The start retention time (in minutes).
        rtend (float): The end retention time (in minutes).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified retention time range.
    """
    run = pymzml.run.Reader(file, build_index_from_scratch=True)  
    scan_dfs = []
    for spectrum in run:
        rt_val = spectrum.scan_time_in_minutes()
        if(rtstart < rt_val < rtend):
            df_scan = pd.DataFrame({'mz':spectrum.mz, 'int':spectrum.i, 'rt':[rt_val]*len(spectrum.i)})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))


