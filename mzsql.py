
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
        bet_idxs = (mzmin < spectrum["m/z array"]) & (spectrum["m/z array"] < mzmax)
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
        bet_idxs = mzmin < mz_vals < mzmax
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
        bet_idxs = mzmin < mz_vals < mzmax
    
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
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    exp.updateRanges()
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
        if(rtstart < rt_val < rtend):
            df_scan = pd.DataFrame({'mz':spectrum.mz, 'int':spectrum.i, 'rt':[rt_val]*len(spectrum.i)})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
def get_rtrange_mzml_pyopenms_2DPeak(file, rtstart, rtend):
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    exp.updateRanges()
    rtrange_data=exp.get2DPeakDataLong(min_mz=exp.getMinMZ(), max_mz=exp.getMaxMZ(), min_rt=rtstart*60, max_rt=rtend*60)
    return(pd.DataFrame({"rt":rtrange_data[0], "mz":rtrange_data[1], "int":rtrange_data[2]}))





# Indexed mzML things --------------------------------------------------------------------------
def get_chrom_mzml_idx(idx_file, mz, ppm):
    raise Exception("Indexed mzML files not currently supported")
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    # file_data = file_data=mzml.PreIndexedMzML(idx_file).build_byte_index()
    # with mzml.read(idx_file, use_index=True) as reader
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
    raise Exception("Indexed mzML files not currently supported")
    # file_data=mzml.PreIndexedMzML(idx_file).build_byte_index()
    # with mzml.read(idx_file, use_index=True) as reader:
    return(pd.DataFrame({"mz":file_data[scan_num]['m/z array'], "int":file_data[scan_num]['intensity array']}))


def get_rtrange_mzml_idx(idx_file, rtstart, rtend):
    raise Exception("Indexed mzML files not currently supported")
    scan_dfs = []
    file_data=mzml.PreIndexedMzML(idx_file).build_byte_index()
    for spectrum in file_data:
        rt_val = spectrum['scanList']['scan'][0]['scan start time']
        if(rtstart < rt_val < rtend):
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']
            df_scan = pd.DataFrame({'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(mz_vals)})
            scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))



# mzMLb things ---------------------------------------------------------------------------------
def get_chrom_mzmlb(file, mz, ppm):
    scan_dfs = []
    for spectrum in mzmlb.MzMLb(file):
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
    file_data = mzmlb.MzMLb(file)
    return(pd.DataFrame({"mz":file_data[scan_num]['m/z array'], "int":file_data[scan_num]['intensity array']}))


def get_rtrange_mzmlb(file, rtstart, rtend):
    scan_dfs = []
    for spectrum in mzmlb.MzMLb(file):
        rt_val = spectrum['scanList']['scan'][0]['scan start time']
        if(rtstart < rt_val < rtend):
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']
            df_scan = pd.DataFrame({'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(mz_vals)})
            scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))







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

def get_spec_mz5(file, scan_num):
    # Double-check that pyteomics doesn't have a function for this already
    mz5_file = h5py.File(file, 'r')
    scan_idxs = np.concatenate(([0], mz5_file["SpectrumIndex"][...]))
    lower_bound = scan_idxs[scan_num]
    upper_bound = scan_idxs[scan_num+1]
    mz_vals = np.cumsum(mz5_file["SpectrumMZ"][lower_bound:upper_bound])
    int_vals = mz5_file["SpectrumIntensity"][lower_bound:upper_bound]
    return(pd.DataFrame({"mz":mz_vals, "int":int_vals}))

def get_rtrange_mz5(file, rtstart, rtend):
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
def turn_mzml_sqlite(file, outfile):
    #converts mzml file to sqlite database. Here for reference, and to show change for row id
    conn = sqlite3.connect(outfile)
    for spectrum in mzml.MzML(file):
        if spectrum['ms level'] == 1:
            #print(spectrum["index"])
            idx = int(spectrum['id'].split("scan=")[-1].split()[0])
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            df_scan = pd.DataFrame({'id':idx,'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(mz_vals)})
            df_scan.to_sql("MS1", conn, if_exists="append", index=False)
    conn.close()
    return(outfile)

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



