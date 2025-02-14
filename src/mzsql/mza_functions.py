
import numpy as np
import pandas as pd
import h5py
import mzapy
from .helpers import pmppm

def get_chrom_mza(file, mz, ppm):
    """
    Retrieves chromatogram data from an mzapy MZA file for a specified m/z range.

    Args:
        file (str): Path to the mzapy file to be processed.
        mz (float): The target m/z value.
        ppm (float): The allowed mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing retention time and intensity data 
                      for the specified m/z range.
    """
    mzmin, mzmax = pmppm(mz, ppm)
    mza=mzapy.MZA(file)
    xic_rt, xic_int = mza.collect_xic_arrays_by_mz(mzmin, mzmax)
    mza.close()
    chrom_data = pd.DataFrame({"rt":xic_rt, "int":xic_int})
    return(chrom_data)
    
def get_rtrange_mza(file, rtstart, rtend):
    """
    Retrieves spectrum data within a specified retention time range from an mzapy MZA file.

    Args:
        file (str): Path to the mzapy file.
        rtstart (float): The start retention time value.
        rtend (float): The end retention time value.

    Returns:
        pd.DataFrame: A DataFrame containing retention time, m/z, and intensity values for 
                      spectra within the specified retention time range.
    """
    mza=mzapy.MZA(file)
    rtrange_data = mza.collect_ms1_df_by_rt(rtstart, rtend)
    rtrange_data = rtrange_data[["rt", "mz", "intensity"]]
    rtrange_data.columns = ["rt", "mz", "int"]
    mza.close()
    return(rtrange_data)
    
def get_spec_mza(file, spectrum_idx):
    """
    Retrieves spectrum data for a specific scan index from an MZA file.

    Args:
        file (str): Path to the MZA file.
        spectrum_idx (int): The scan index of the spectrum to retrieve.

    Returns:
        pd.DataFrame: A DataFrame containing m/z and intensity values for the specified spectrum.
    """
    mza = h5py.File(file, 'r')
    intensities = mza["Arrays_intensity/"+str(spectrum_idx)][:]
    mz = mza["Arrays_mz/"+str(spectrum_idx)][:]
    mza.close()
    spec_df = pd.DataFrame({"mz":mz, "int":intensities})
    return(spec_df)

def get_MS2premz_mza(mza_file, precursor_mz, ppm_acc):
    mzmin, mzmax = pmppm(precursor_mz, ppm_acc)
    scan_dfs = []
    mza = h5py.File(mza_file, 'r')
    for entry in mza["Metadata"]:
        if entry['MSLevel'] == 2 and mzmin <= entry['IsolationWindowTargetMz'] <= mzmax:
            scan_df=pd.DataFrame({
                    "rt": entry['RetentionTime'],
                    "premz": entry['IsolationWindowTargetMz'],
                    "fragmz": mza["Arrays_mz/"+str(entry["Scan"])][:],
                    "int": mza["Arrays_intensity/"+str(entry["Scan"])][:]
                })
            scan_dfs.append(scan_df)
    mza.close()
    spec_data = pd.concat(scan_dfs, ignore_index=True)
    return(spec_data)

def get_MS2fragmz_mza(mza_file, fragment_mz, ppm_acc):
    mzmin, mzmax = pmppm(fragment_mz, ppm_acc)
    scan_dfs = []
    mza = h5py.File(mza_file, 'r')

    file_keys = sorted(mza["Arrays_intensity"].keys(), key=lambda x: int(x))
    for index, scan_num in enumerate(file_keys):
        if mza["Metadata"][index]["MSLevel"] == 2:
            frag_mz_idxs = (mzmin < mza[f"Arrays_mz/{scan_num}"][...]) & (mza[f"Arrays_mz/{scan_num}"][...] < mzmax)
            mz_vals = mza[f"Arrays_mz/{scan_num}"][...][frag_mz_idxs]
            int_vals = mza[f"Arrays_intensity/{scan_num}"][...][frag_mz_idxs]
            
            scan_df=pd.DataFrame({
                "rt": mza["Metadata"][index]["RetentionTime"],
                "premz": mza["Metadata"][index]["IsolationWindowTargetMz"],
                "fragmz": mz_vals,
                "int": int_vals
            })
            scan_dfs.append(scan_df)
    mza.close()
    spec_data = pd.concat(scan_dfs, ignore_index=True)
    return(spec_data)

def get_MS2nloss_mza(mza_file, neutral_loss, ppm_acc):
    scan_dfs = []
    mza = h5py.File(mza_file, 'r')

    file_keys = sorted(mza["Arrays_intensity"].keys(), key=lambda x: int(x))
    for index, scan_num in enumerate(file_keys):
        if mza["Metadata"][index]["MSLevel"] == 2:
            mzmin, mzmax = pmppm(neutral_loss, ppm_acc)
            init_mz_vals = mza["Metadata"][index]["IsolationWindowTargetMz"]-mza[f"Arrays_mz/{scan_num}"][...]
            nloss_idxs = (mzmin < init_mz_vals) & (init_mz_vals < mzmax)
            mz_vals = mza[f"Arrays_mz/{scan_num}"][...][nloss_idxs]
            int_vals = mza[f"Arrays_intensity/{scan_num}"][...][nloss_idxs]
            scan_df=pd.DataFrame({
                "rt": mza["Metadata"][index]["RetentionTime"],
                "premz": mza["Metadata"][index]["IsolationWindowTargetMz"],
                "fragmz": mz_vals,
                "int": int_vals
            })
            scan_dfs.append(scan_df)
    mza.close()
    spec_data = pd.concat(scan_dfs, ignore_index=True)
    return(spec_data)
