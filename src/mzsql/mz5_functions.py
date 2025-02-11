
import numpy as np
import pandas as pd
import h5py
from .helpers import pmppm


def get_chrom_mz5(file, mz, ppm):
    mz5_file = h5py.File(file, 'r')
    bounds_df = pd.DataFrame({"lower":np.concatenate(([0], mz5_file["SpectrumIndex"][...][:-1])),
                              "upper":mz5_file["SpectrumIndex"][...]})
    
    mzmin, mzmax = pmppm(mz, ppm)
    has_precursor = [len(item)>0 for item in mz5_file["SpectrumMetaData"]["precursors"]]
    scan_dfs = []
    for index, row in bounds_df.iterrows():
        if(not(has_precursor[index])):
            scan_df = pd.DataFrame({
                "rt": mz5_file["ChomatogramTime"][index],
                "mz": np.cumsum(mz5_file["SpectrumMZ"][row["lower"]:row["upper"]]),
                "int": mz5_file["SpectrumIntensity"][row["lower"]:row["upper"]]
            })
            bet_df = scan_df[(scan_df["mz"]>mzmin) & (scan_df["mz"]<mzmax)]
            scan_dfs.append(bet_df)
    file_df = pd.concat(scan_dfs, ignore_index=True)
    mz5_file.close()
    return(file_df)

def get_spec_mz5(file, scan_num):
    mz5_file = h5py.File(file, 'r')
    mz5_meta = mz5_file["SpectrumMetaData"]["id", "index"]
    idx = np.where(np.char.find(mz5_meta['id'].astype(str), f'scan={scan_num}') != -1)[0][0]
    scan_idxs = np.concatenate(([0], mz5_file["SpectrumIndex"][...]))
    lower_bound = scan_idxs[idx-1]
    upper_bound = scan_idxs[idx]
    mz_vals = np.cumsum(mz5_file["SpectrumMZ"][lower_bound:upper_bound])
    int_vals = mz5_file["SpectrumIntensity"][lower_bound:upper_bound]
    return(pd.DataFrame({"mz":mz_vals, "int":int_vals}))

def get_rtrange_mz5(file, rtstart, rtend):
    mz5_file = h5py.File(file, 'r')
    scan_idxs = np.concatenate(([0], mz5_file["SpectrumIndex"][...]))
    scan_dfs = []
    has_precursor = [len(item)>0 for item in mz5_file["SpectrumMetaData"]["precursors"]]
    
    for index, rt_val in enumerate(mz5_file["ChomatogramTime"][...]):
        if(rtstart < rt_val < rtend):
            if(not(has_precursor[index])):
                scan_df = pd.DataFrame({
                    "rt": rt_val,
                    "mz": np.cumsum(mz5_file["SpectrumMZ"][scan_idxs[index]:scan_idxs[index+1]]),
                    "int": mz5_file["SpectrumIntensity"][scan_idxs[index]:scan_idxs[index+1]]
                })
                scan_dfs.append(scan_df)
    file_df = pd.concat(scan_dfs, ignore_index=True)
    mz5_file.close()
    return(file_df)