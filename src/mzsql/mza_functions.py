
import numpy as np
import pandas as pd
import h5py
import mzapy
from .helpers import pmppm

def get_chrom_mza(file, mz, ppm):
    mza = h5py.File(file, 'r')
    
    scan_dfs = []
    file_keys = sorted(mza["Arrays_intensity"].keys(), key=lambda x: int(x))
    for index, scan_num in enumerate(file_keys):
        scan_df=pd.DataFrame({
            "rt": mza["Metadata"][index][6],
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
    mza = h5py.File(file, 'r')
    intensities = mza["Arrays_intensity/"+str(spectrum_idx)][:]
    mz = mza["Arrays_mz/"+str(spectrum_idx)][:]
    mza.close()
    spec_df = pd.DataFrame({"mz":mz, "int":intensities})
    return(spec_df)

def get_rtrange_mza(file, rtstart, rtend):
    mza = h5py.File(file, 'r')
    scan_dfs = []
    file_keys = sorted(mza["Arrays_intensity"].keys(), key=lambda x: int(x))
    for index, scan_num in enumerate(file_keys):
        rt_val = mza["Metadata"][index][6]
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
    mzmin, mzmax = pmppm(mz, ppm)
    mza=mzapy.MZA(file)
    xic_rt, xic_int = mza.collect_xic_arrays_by_mz(mzmin, mzmax)
    mza.close()
    chrom_data = pd.DataFrame({"rt":xic_rt, "int":xic_int})
    return(chrom_data)
    
def get_spec_mzapy(file, spectrum_idx):
    mza=mzapy.MZA(file)
    mean_rt_diff = np.mean(np.diff(mza.rt))
    rt_min, rt_max = mza.rt[spectrum_idx-1]-mean_rt_diff/100, mza.rt[spectrum_idx-1]+mean_rt_diff/100
    mz_array, int_array = mza.collect_ms1_arrays_by_rt(rt_min, rt_max)
    mza.close()
    spec_df = pd.DataFrame({"mz":mz_array, "int":int_array})
    return(spec_df)

def get_rtrange_mzapy(file, rtstart, rtend):
    mza=mzapy.MZA(file)
    rtrange_data = mza.collect_ms1_df_by_rt(rtstart, rtend)
    rtrange_data = rtrange_data[["rt", "mz", "intensity"]]
    rtrange_data.columns = ["rt", "mz", "int"]
    mza.close()
    return(rtrange_data)
