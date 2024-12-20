
import numpy as np
import pandas as pd
import pyteomics.mzml
from .helpers import pmppm


def get_chrom_mzml_idx(idx_file, mz, ppm):
    """
    Retrieves chromatogram data from an indexed mzML file for a specific m/z range.

    Args:
        idx_file (str): Path to the indexed mzML file.
        mz (float): The target m/z value.
        ppm (float): The mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified m/z range.
    """
    raise Exception("Indexed mzML files not currently supported")
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    # file_data = file_data=pyteomics.mzml.PreIndexedMzML(idx_file).build_byte_index()
    # with pyteomics.mzml.read(idx_file, use_index=True) as reader
    for spectrum in file_data:
        rt_val = spectrum['scanList']['scan'][0]['scan start time']
        mz_vals=spectrum['m/z array']
        int_vals = spectrum['intensity array']
        bet_idxs = (mzmin < spectrum["m/z array"]) & (spectrum["m/z array"] < mzmax)
        if(sum(bet_idxs)>0):
            df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
            scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))


def get_spec_mzml_idx(idx_file, scan_num):
    """
    Retrieves spectrum data for a specific scan number from an indexed mzML file.

    Args:
        idx_file (str): Path to the indexed mzML file.
        scan_num (int): The scan number.

    Returns:
        pd.DataFrame: A DataFrame containing m/z and intensity data for the specified scan number.
    """
    raise Exception("Indexed mzML files not currently supported")
    # file_data=pyteomics.mzml.PreIndexedMzML(idx_file).build_byte_index()
    # with pyteomics.mzml.read(idx_file, use_index=True) as reader:
    return(pd.DataFrame({"mz":file_data[scan_num]['m/z array'], "int":file_data[scan_num]['intensity array']}))
    
def get_rtrange_mzml_idx(idx_file, rtstart, rtend):
    """
    Retrieves chromatogram data from an indexed mzML file for a specific retention time range.

    Args:
        idx_file (str): Path to the indexed mzML file.
        rtstart (float): The start retention time (in minutes).
        rtend (float): The end retention time (in minutes).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified retention time range.
    """
    raise Exception("Indexed mzML files not currently supported")
    scan_dfs = []
    file_data=pyteomics.mzml.PreIndexedMzML(idx_file).build_byte_index()
    for spectrum in file_data:
        rt_val = spectrum['scanList']['scan'][0]['scan start time']
        if(rtstart < rt_val < rtend):
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']
            df_scan = pd.DataFrame({'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(mz_vals)})
            scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))