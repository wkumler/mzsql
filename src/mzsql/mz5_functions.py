
import numpy as np
import pandas as pd
import h5py
from .helpers import pmppm


def get_chrom_mz5(file, mz, ppm):
    """
    Retrieves chromatogram data from an mz5 file for a specified m/z range.

    Args:
        file (str): Path to the mz5 file to be processed.
        mz (float): The target m/z value.
        ppm (float): The allowed mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing m/z, intensity, and retention time data 
                      for the specified m/z range.
    """
    mz5_file = h5py.File(file, 'r')
    bounds_df = pd.DataFrame({"lower":np.concatenate(([0], mz5_file["SpectrumIndex"][...][:-1])),
                              "upper":mz5_file["SpectrumIndex"][...]})
    
    scan_dfs = []
    for index, row in bounds_df.iterrows():
        scan_df = pd.DataFrame({
            "rt": mz5_file["ChomatogramTime"][...][index],
            "mz": np.cumsum(mz5_file["SpectrumMZ"][row["lower"]:row["upper"]]),
            "int": mz5_file["SpectrumIntensity"][row["lower"]:row["upper"]]
        })
        scan_dfs.append(scan_df)
    file_df = pd.concat(scan_dfs, ignore_index=True)
    
    mzmin, mzmax = pmppm(mz, ppm)
    bet_df = file_df[(file_df["mz"]>mzmin) & (file_df["mz"]<mzmax)]
    
    mz5_file.close()
    return(bet_df)

def get_spec_mz5(file, scan_num):
    """
    Retrieves spectrum data for a specific scan number from an mz5 file.

    Args:
        file (str): Path to the mz5 file.
        scan_num (int): The scan number to retrieve.

    Returns:
        pd.DataFrame: A DataFrame containing m/z and intensity values for the specified scan.
    """
    # Double-check that pyteomics doesn't have a function for this already
    mz5_file = h5py.File(file, 'r')
    scan_idxs = np.concatenate(([0], mz5_file["SpectrumIndex"][...]))
    lower_bound = scan_idxs[scan_num]
    upper_bound = scan_idxs[scan_num+1]
    mz_vals = np.cumsum(mz5_file["SpectrumMZ"][lower_bound:upper_bound])
    int_vals = mz5_file["SpectrumIntensity"][lower_bound:upper_bound]
    return(pd.DataFrame({"mz":mz_vals, "int":int_vals}))

def get_rtrange_mz5(file, rtstart, rtend):
    """
    Retrieves spectrum data within a specified retention time range from an mz5 file.
    
    Note:
        This function may not yet be fully functional and currently returns scans
        without relevant spectra.

    Args:
        file (str): Path to the mz5 file.
        rtstart (float): The start retention time value.
        rtend (float): The end retention time value.

    Returns:
        pd.DataFrame: A DataFrame containing m/z and intensity values for spectra 
                      within the specified retention time range.
    """
    print("Warning: get_rtrange_mz5 does not yet seem to be functional")
    print("Warning: returns scans without relevant spectra")
    mz5_file = h5py.File(file, 'r')
    scan_idxs = np.concatenate(([0], mz5_file["SpectrumIndex"][...]))
    scan_dfs = []
    for index, rt_val in enumerate(mz5_file["ChomatogramTime"][...]):
        if(rtstart < rt_val < rtend):
            scan_df = pd.DataFrame({
                "rt": rt_val,
                "mz": np.cumsum(mz5_file["SpectrumMZ"][scan_idxs[index]:scan_idxs[index+1]]),
                "int": mz5_file["SpectrumIntensity"][scan_idxs[index]:scan_idxs[index+1]]
            })
            scan_dfs.append(scan_df)
    file_df = pd.concat(scan_dfs, ignore_index=True)
    mz5_file.close()
    return(file_df)