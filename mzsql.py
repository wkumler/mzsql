
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyteomics import mzml, mzmlb
import h5py
import sqlite3
import duckdb

# Helper functions
def pmppm(mass, ppm=4):
    return(mass*(1-ppm/1000000), mass*(1+ppm/1000000))


# mzML things
def get_chrom_mzml(file, mz, ppm):
    scan_dfs = []
    for spectrum in mzml.MzML(file):
        if spectrum['ms level'] == 1:
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']

            mzmin, mzmax = pmppm(mz, ppm)
            bet_idxs = (spectrum["m/z array"]>mzmin) & (spectrum["m/z array"]<mzmax)
            if(sum(bet_idxs)>0):
                df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
                scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))

def get_spectrum_mzml(file, spectrum_idx):
    mzml.MzML(file)[spectrum_idx]

def get_rtrange_mzml(file, rtstart, rtend):
    raise Exception("mzML rtrange extraction yet implemented")


# Indexed mzML things
def get_chrom_mzml_idx(idx_file, mz, ppm):
    get_chrom_mzml(file, mz, ppm)

def get_spectrum_mzml_idx(idx_file, spectrum_idx):
    # I think this one should be fancier than for default mzML extraction bc index exists but pyteomics should handle it well
    mzml.MzML(file)[100]

def get_rtrange_mzml_idx(idx_file, rtstart, rtend):
    raise Exception("Indexed mzML rtrange extraction yet implemented")



# mzMLb things
def get_chrom_mzmlb(file, mz, ppm):
    scan_dfs = []
    for spectrum in mzmlb.MzMLb(file):
        if spectrum['ms level'] == 1:
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']

            mzmin, mzmax = pmppm(mz, ppm)
            bet_idxs = (spectrum["m/z array"]>mzmin) & (spectrum["m/z array"]<mzmax)
            if(sum(bet_idxs)>0):
                df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
                scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))

def get_spectrum_mzml(file, spectrum_idx):
    raise Exception("mzMLb spectrum extraction yet implemented")

def get_rtrange_mzml(file, rtstart, rtend):
    raise Exception("mzMLb rtrange extraction yet implemented")







# mz5 things
def get_chrom_mz5(file, mz, ppm):
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

def get_spectrum_mz5(file, spectrum_idx):
    # I think this one should be fancier than for default mzML extraction bc index exists but pyteomics should handle it well
    raise Exception("mz5 spectrum extraction yet implemented")

def get_rtrange_mz5(file, rtstart, rtend):
    raise Exception("mz5 rtrange extraction yet implemented")




# MZA things
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
    
def get_spectrum_mza(file, spectrum_idx):
    mza = h5py.File(mzafile, 'r')
    intensities = mza["Arrays_intensity/1001"][:]
    mz = mza["Arrays_mz/1001"][:]

    mza.close()
    spec_df = pd.DataFrame({"mz":mz, "int":intensities})

def get_rtrange_mza(file, rtstart, rtend):
    raise Exception("mz5 rtrange extraction yet implemented")




# SQLite things
# Note - can be parallelized across files so the syntax may differ slightly
def get_chrom_sqlite(file, mz, ppm):
    raise Exception("SQLite spectrum extraction yet implemented")

def get_spectrum_sqlite(file, spectrum_idx):
    raise Exception("SQLite spectrum extraction yet implemented")

def get_rtrange_sqlite(file, rtstart, rtend):
    raise Exception("SQLite rtrange extraction yet implemented")




# DuckDB things
# Likely to be very similar to the SQLite thing
def get_chrom_duckdb(file, mz, ppm):
    raise Exception("DuckDB spectrum extraction yet implemented")

def get_spectrum_duckdb(file, spectrum_idx):
    raise Exception("DuckDB spectrum extraction yet implemented")

def get_rtrange_duckdb(file, rtstart, rtend):
    raise Exception("DuckDB rtrange extraction yet implemented")



