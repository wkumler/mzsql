
import numpy as np
import pandas as pd
import pyteomics.mzmlb
from .helpers import pmppm

def get_chrom_mzmlb(file, mz, ppm):
    """
    Retrieves chromatogram data from an mzMLb file for a specific m/z range.

    Args:
        file (str): Path to the mzMLb file.
        mz (float): The target m/z value.
        ppm (float): The mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified m/z range.
    """
    scan_dfs = []
    for spectrum in pyteomics.mzmlb.MzMLb(file):
        if spectrum['ms level'] == 1:
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']

            mzmin, mzmax = pmppm(mz, ppm)
            bet_idxs = (mzmin < spectrum["m/z array"]) & (spectrum["m/z array"] < mzmax)
            if(sum(bet_idxs)>0):
                df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
                scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))

def get_spec_mzmlb(file, scan_num):
    """
    Retrieves spectrum data for a specific scan number from an mzMLb file.

    Args:
        file (str): Path to the mzMLb file.
        scan_num (int): The scan number.

    Returns:
        pd.DataFrame: A DataFrame containing m/z and intensity data for the specified scan number.
    """
    file_data = pyteomics.mzmlb.MzMLb(file)
    return(pd.DataFrame({"mz":file_data[scan_num]['m/z array'], "int":file_data[scan_num]['intensity array']}))


def get_rtrange_mzmlb(file, rtstart, rtend):
    """
    Retrieves chromatogram data from an mzMLb file for a specific retention time range.

    Args:
        file (str): Path to the mzMLb file.
        rtstart (float): The start retention time (in minutes).
        rtend (float): The end retention time (in minutes).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified retention time range.
    """
    scan_dfs = []
    for spectrum in pyteomics.mzmlb.MzMLb(file):
        rt_val = spectrum['scanList']['scan'][0]['scan start time']
        if(rtstart < rt_val < rtend):
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']
            df_scan = pd.DataFrame({'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(mz_vals)})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))

def get_MS2scan_mzmlb(mzmlb_file, scan_num):
    file_data = pyteomics.mzmlb.MzMLb(mzmlb_file)
    rt_val = file_data[scan_num]['scanList']['scan'][0]['scan start time']
    premz_val = file_data[scan_num]["precursorList"]["precursor"][0]["isolationWindow"]['isolation window target m/z']
    return(pd.DataFrame({"rt":rt_val, "premz":premz_val, "fragmz":file_data[scan_num]['m/z array'], "int":file_data[scan_num]['intensity array']}))

def get_MS2premz_mzmlb(mzmlb_file, precursor_mz, ppm):
    mzmin, mzmax = pmppm(precursor_mz, ppm_acc)
    scan_dfs = []
    for spectrum in pyteomics.mzmlb.MzMLb(mzmlb_file):
        if spectrum['ms level'] == 2:
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            premz_val = spectrum["precursorList"]["precursor"][0]["isolationWindow"]['isolation window target m/z']
            if mzmin < premz_val < mzmax:
                mz_vals=spectrum['m/z array']
                int_vals = spectrum['intensity array']
                df_scan = pd.DataFrame({'rt':rt_val, 'premz':premz_val, 'fragmz':mz_vals, 'int':int_vals})
                scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))

def get_MS2fragmz_mzmlb(mzmlb_file, fragment_mz, ppm):
    mzmin, mzmax = pmppm(fragment_mz, ppm_acc)
    scan_dfs = []
    for spectrum in pyteomics.mzmlb.MzMLb(mzmlb_file):
        if spectrum['ms level'] == 2:
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            premz_val = spectrum["precursorList"]["precursor"][0]["isolationWindow"]['isolation window target m/z']
            chosen_frag_idxs = (mzmin < spectrum['m/z array']) & (spectrum['m/z array'] < mzmax)
            mz_vals=spectrum['m/z array'][chosen_frag_idxs]
            int_vals = spectrum['intensity array'][chosen_frag_idxs]
            df_scan = pd.DataFrame({'rt':rt_val, 'premz':premz_val, 'fragmz':mz_vals, 'int':int_vals})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))

def get_MS2nloss_mzmlb(mzmlb_file, neutral_loss, ppm):
    mzmin, mzmax = pmppm(neutral_loss, ppm_acc)
    scan_dfs = []
    for spectrum in pyteomics.mzmlb.MzMLb(mzmlb_file):
        if spectrum['ms level'] == 2:
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            premz_val = spectrum["precursorList"]["precursor"][0]["isolationWindow"]['isolation window target m/z']
            chosen_frag_idxs = (mzmin < premz_val-spectrum['m/z array']) & (premz_val-spectrum['m/z array'] < mzmax)
            mz_vals=spectrum['m/z array'][chosen_frag_idxs]
            int_vals = spectrum['intensity array'][chosen_frag_idxs]
            df_scan = pd.DataFrame({'rt':rt_val, 'premz':premz_val, 'fragmz':mz_vals, 'int':int_vals})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
