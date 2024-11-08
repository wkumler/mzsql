
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyteomics import mzml, mzmlb
import pyopenms
import pymzml
import h5py
import sqlite3
import duckdb

# Helper functions
def pmppm(mass, ppm=4):
    return(mass*(1-ppm/1000000), mass*(1+ppm/1000000))


# mzML things ----------------------------------------------------------------------------
def get_chrom_mzml_pyteomics(file, mz, ppm):
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    for spectrum in mzml.MzML(file):
        rt_val = spectrum['scanList']['scan'][0]['scan start time']
        mz_vals=spectrum['m/z array']
        int_vals = spectrum['intensity array']
        bet_idxs = (spectrum["m/z array"]>mzmin) & (spectrum["m/z array"]<mzmax)
        if(sum(bet_idxs)>0):
            df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
            scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))
def get_chrom_mzml_pyopenms(file, mz, ppm):
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    for spectrum in exp:
        rt_val = spectrum.getRT()
        mz_vals, int_vals = spectrum.get_peaks()
        bet_idxs = (mz_vals>mzmin) & (mz_vals<mzmax)
        if(sum(bet_idxs)>0):
            df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
            scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))
def get_chrom_mzml_pymzml(file, mz, ppm):
    run = pymzml.run.Reader(file)
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    for spectrum in run:
        rt_val = spectrum.scan_time_in_minutes()
        mz_vals = spectrum.mz
        int_vals = spectrum.i
        bet_idxs = (mz_vals>mzmin) & (mz_vals<mzmax)
    
        if(sum(bet_idxs)>0):
            df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
            scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))
def get_chrom_mzml_pyopenms_2DPeak(file, mz, ppm):
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    exp.updateRanges()
    mzmin, mzmax = pmppm(mz, ppm)
    chrom_data=exp.get2DPeakDataLong(min_mz=mzmin, max_mz=mzmax, min_rt=exp.getMinRT(), max_rt=exp.getMaxRT())
    return(pd.DataFrame({"rt":chrom_data[0], "mz":chrom_data[1], "int":chrom_data[2]}))



def get_spec_mzml_pyteomics(file, scan_num):
    file_data = mzml.MzML(file)
    return(pd.DataFrame({"mz":file_data[scan_num]['m/z array'], "int":file_data[scan_num]['intensity array']}))
def get_spec_mzml_pyopenms(file, scan_num):
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    spec1_data = exp[scan_num].get_peaks()
    return(pd.DataFrame({"mz":spec1_data[0], "int":spec1_data[1]}))
def get_spec_mzml_pymzml(file, scan_num):
    run = pymzml.run.Reader(file)
    spec1_data = run[scan_num].peaks("raw")
    return(pd.DataFrame({"mz":spec1_data[:,0], "int":spec1_data[:,1]}))




def get_rtrange_mzml_pyteomics(file, rtstart, rtend):
    scan_dfs = []
    for spectrum in mzml.MzML(file):
        rt_val = spectrum['scanList']['scan'][0]['scan start time']
        if(rtstart < rt_val < rtend):
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']
            df_scan = pd.DataFrame({'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(mz_vals)})
            scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))
def get_rtrange_mzml_pyopenms(file, rtstart, rtend):
    scan_dfs = []
    for spectrum in exp:
        rt_val = spectrum.getRT()
        if(rtstart*60 < rt_val < rtend*60):
            mz_vals, int_vals = spectrum.get_peaks()
            df_scan = pd.DataFrame({'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(int_vals)})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
def get_rtrange_mzml_pymzml(file, rtstart, rtend):
    run = pymzml.run.Reader(file)  
    scan_dfs = []
    for spectrum in run:
        rt_val = spectrum.scan_time_in_minutes()
        if(rtstart<rt_val<rtend):
            df_scan = pd.DataFrame({'mz':spectrum.mz, 'int':spectrum.i, 'rt':[rt_val]*len(spectrum.i)})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
def get_rtrange_mzml_pyopenms_2DPeak(file, rtstart, rtend):
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    exp.updateRanges()
    rtrange_data=exp.get2DPeakDataLong(min_mz=exp.getMinMZ(), max_mz=exp.getMaxMZ(), min_rt=rtstart, max_rt=rtend)
    return(pd.DataFrame({"rt":rtrange_data[0], "mz":rtrange_data[1], "int":rtrange_data[2]}))






def get_rtrange_mzml(file, rtstart, rtend):
    raise Exception("mzML rtrange extraction yet implemented")


# Indexed mzML things --------------------------------------------------------------------------
def get_chrom_mzml_idx(idx_file, mz, ppm):
    get_chrom_mzml(file, mz, ppm)

def get_spectrum_mzml_idx(idx_file, spectrum_idx):
    # I think this one should be fancier than for default mzML extraction bc index exists but pyteomics should handle it well
    mzml.MzML(file)[100]

def get_rtrange_mzml_idx(idx_file, rtstart, rtend):
    raise Exception("Indexed mzML rtrange extraction yet implemented")



# mzMLb things ---------------------------------------------------------------------------------
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







# mz5 things ---------------------------------------------------------------------------------------
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




# MZA things -----------------------------------------------------------------------------------
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




# SQLite things ------------------------------------------------------------------------------------
# Note - can be parallelized across files so the syntax may differ slightly
import sqlite3
import pandas as pd

import sqlite3
import pandas as pd

def turn_mzml_sqlite(file):
    #converts mzml file to sqlite database. Here for reference, and to show change for row id
    conn = sqlite3.connect("msdata.sqlite")
    for spectrum in mzml.MzML(file):
        if spectrum['ms level'] == 1:
            print(spectrum["index"])
            idx = "index"
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            df_scan = pd.DataFrame({'id':idx,'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(mz_vals)})
            df_scan.to_sql("MS1", conn, if_exists="append", index=False)
conn.close()

def get_chrom_sqlite(file, mz, ppm):
    #returns the chromatogram from the given file. 
    #Uses mz and ppm to find mz range that will be rturned.
    
    mz_min, mz_max = pmppm(mz, ppm)
    
    conn = sqlite3.connect(file)
    query = f"SELECT * FROM MS1 WHERE mz BETWEEN {mz_min} AND {mz_max}"
    query_data = pd.read_sql_query(query, conn)
    
    conn.close()
    return query_data

def get_spectrum_sqlite(file, spectrum_idx):
    #return a single spectrum according to ID
    conn = sqlite3.connect(file)

    query = f"SELECT * FROM MS1 WHERE id = {spectrum_idx}"
    spectrum_data = pd.read_sql_query(query, conn)

    conn.close()
    
    return spectrum_data

def get_rtrange_sqlite(file, rtstart, rtend):
    #return table between rtstart and rtend:
    conn = sqlite3.connect(file)
    query = f"SELECT * FROM MS1 WHERE rt >= {rtstart} AND rt <= {rtend}"
    rt_range_data = pd.read_sql_query(query, conn)
    
    conn.close()
    
    return rt_range_data




# DuckDB things
# Likely to be very similar to the SQLite thing
def get_chrom_duckdb(file, mz, ppm):
    raise Exception("DuckDB spectrum extraction yet implemented")

def get_spectrum_duckdb(file, spectrum_idx):
    raise Exception("DuckDB spectrum extraction yet implemented")

def get_rtrange_duckdb(file, rtstart, rtend):
    raise Exception("DuckDB rtrange extraction yet implemented")



