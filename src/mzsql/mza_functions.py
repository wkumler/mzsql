
import numpy as np
import pandas as pd
import h5py
import mzapy

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
    mza = h5py.File(mzafile, 'r')
    intensities = mza["Arrays_intensity/"+scan_num][:]
    mz = mza["Arrays_mz/"+scan_num][:]

    mza.close()
    spec_df = pd.DataFrame({"mz":mz, "int":intensities})

def get_rtrange_mza(file, rtstart, rtend):
    raise Exception("mza rtrange extraction yet implemented")



def get_chrom_mzapy(file, mz, ppm):
    raise Exception("mza chrom extraction via mzapy not yet implemented")
    
def get_spec_mza(file, spectrum_idx):
    raise Exception("mza spec extraction via mzapy not yet implemented")

def get_rtrange_mza(file, rtstart, rtend):
    raise Exception("mza rtrange extraction via mzapy yet implemented")

