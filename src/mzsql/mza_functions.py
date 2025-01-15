
import numpy as np
import pandas as pd
import h5py
import mzapy
from .helpers import pmppm

def get_chrom_mza(file, mz, ppm):
    """
    Retrieves chromatogram data from an MZA file for a specified m/z range.

    Args:
        file (str): Path to the MZA file to be processed.
        mz (float): The target m/z value.
        ppm (float): The allowed mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing m/z, intensity, and retention time data 
                      for the specified m/z range.
    """
    mza = h5py.File(file, 'r')
    
    scan_dfs = []
    file_keys = sorted(mza["Arrays_intensity"].keys(), key=lambda x: int(x))
    for index, scan_num in enumerate(file_keys):
        if mza["Metadata"][index]["MSLevel"] == 1:
            scan_df=pd.DataFrame({
                "rt": mza["Metadata"][index]["RetentionTime"],
                "mz": mza["Arrays_mz/"+scan_num][...],
                "int": mza["Arrays_intensity/"+scan_num][...]
            })
            scan_dfs.append(scan_df)
    mza.close()
    
    file_df = pd.concat(scan_dfs, ignore_index=True)
    
    mzmin, mzmax = pmppm(mz, ppm)
    chrom_data = file_df[(file_df["mz"]>mzmin) & (file_df["mz"]<mzmax)]
    return(chrom_data)
    
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

def get_rtrange_mza(file, rtstart, rtend):
    """
    Retrieves spectrum data within a specified retention time range from an MZA file.

    Args:
        file (str): Path to the MZA file.
        rtstart (float): The start retention time value.
        rtend (float): The end retention time value.

    Returns:
        pd.DataFrame: A DataFrame containing m/z, intensity, and retention time values for 
                      spectra within the specified retention time range.
    """
    mza = h5py.File(file, 'r')
    scan_dfs = []
    file_keys = sorted(mza["Arrays_intensity"].keys(), key=lambda x: int(x))
    for index, scan_num in enumerate(file_keys):
        if mza["Metadata"][index]["MSLevel"] == 1:
            rt_val = mza["Metadata"][index]["RetentionTime"]
            if(rtstart < rt_val < rtend):
                scan_df=pd.DataFrame({
                    "rt": mza["Metadata"][index][6],
                    "mz": mza["Arrays_mz/"+scan_num][...],
                    "int": mza["Arrays_intensity/"+scan_num][...]
                })
                scan_dfs.append(scan_df)
    mza.close()
    rtrange_data = pd.concat(scan_dfs, ignore_index=True)
    return(rtrange_data)



def get_chrom_mzapy(file, mz, ppm):
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
    
def get_spec_mzapy(file, spectrum_idx):
    """
    Retrieves spectrum data for a specific scan index from an mzapy MZA file.

    Args:
        file (str): Path to the mzapy file.
        spectrum_idx (int): The scan index of the spectrum to retrieve.

    Returns:
        pd.DataFrame: A DataFrame containing m/z and intensity values for the specified spectrum.
    """
    mza=mzapy.MZA(file)
    mean_rt_diff = np.mean(np.diff(mza.rt))
    rt_min, rt_max = mza.rt[spectrum_idx-1]-mean_rt_diff/100, mza.rt[spectrum_idx-1]+mean_rt_diff/100
    mz_array, int_array = mza.collect_ms1_arrays_by_rt(rt_min, rt_max)
    mza.close()
    spec_df = pd.DataFrame({"mz":mz_array, "int":int_array})
    return(spec_df)

def get_rtrange_mzapy(file, rtstart, rtend):
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
