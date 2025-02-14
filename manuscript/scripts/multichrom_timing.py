import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mzsql import *
import sqlite3
import os
import random

# Create databases from mzML files if they don't already exist
if(not(os.path.exists("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.sqlite"))):
    turn_mzml_sqlite("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.mzML", "E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.sqlite")
    turn_mzml_duckdb("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.mzML", "E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.duckdb")
    turn_mzml_parquet("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.mzML", "E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01_pqds")

conn = sqlite3.connect("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.sqlite")
top_int_df = pd.read_sql_query("SELECT * FROM MS1 ORDER BY int DESC LIMIT 30000", conn)
top_masses = []
top_rts = []
for i in range(100):
    top_rts.append(top_int_df["rt"].iloc[0])
    top_mz = top_int_df["mz"].iloc[0]
    top_masses.append(top_mz)
    mzmin, mzmax = pmppm(top_mz, 10)
    top_int_df = top_int_df[((top_int_df["mz"]<mzmin) | (top_int_df["mz"]>mzmax))]

def get_multichrom_mzml_pyteomics(file, mzs, ppm):
    file_data = pyteomics.mzml.MzML(file)
    scan_dfs = []
    for spectrum in file_data:
        if spectrum['ms level'] == 1:
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            mz_vals = np.array(spectrum['m/z array'])
            int_vals = np.array(spectrum['intensity array'])
            for mz_i in mzs:
                mzmin, mzmax = pmppm(mz_i, ppm)
                bet_idxs = (mz_vals > mzmin) & (mz_vals < mzmax)
                if np.any(bet_idxs):
                    scan_dfs.append(pd.DataFrame({
                        'mz_i': mz_i,
                        'mz': mz_vals[bet_idxs],
                        'int': int_vals[bet_idxs],
                        'rt': rt_val
                    }))
    return pd.concat(scan_dfs, ignore_index=True)
def get_multichrom_mzmlb(file, mzs, ppm):
    file_data = pyteomics.mzmlb.MzMLb(file)
    scan_dfs = []
    for spectrum in file_data:
        if spectrum['ms level'] == 1:
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            mz_vals = np.array(spectrum['m/z array'])
            int_vals = np.array(spectrum['intensity array'])
            for mz_i in mzs:
                mzmin, mzmax = pmppm(mz_i, ppm)
                bet_idxs = (mz_vals > mzmin) & (mz_vals < mzmax)
                if np.any(bet_idxs):
                    scan_dfs.append(pd.DataFrame({
                        'mz_i': mz_i,
                        'mz': mz_vals[bet_idxs],
                        'int': int_vals[bet_idxs],
                        'rt': rt_val
                    }))
    return pd.concat(scan_dfs, ignore_index=True)
def get_multichrom_mzml_pyopenms(file, mzs, ppm):
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    scan_dfs = []
    for mz_i in mzs:
        mzmin, mzmax = pmppm(mz_i, ppm)
        for spectrum in exp:
            if(spectrum.getMSLevel()==1):
                rt_val = spectrum.getRT()
                mz_vals, int_vals = spectrum.get_peaks()
                bet_idxs = (mzmin < mz_vals) & (mz_vals < mzmax)
                if(sum(bet_idxs)>0):
                    df_scan = pd.DataFrame({'mz_i':mz_i, 'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
                    scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
def get_multichrom_mzml_pyopenms_2Dpeak(file, mzs, ppm):
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(file, exp)
    exp.updateRanges()
    scan_dfs = []
    for mz_i in mzs:
        mzmin, mzmax = pmppm(mz_i, ppm)
        chrom_data=exp.get2DPeakDataLong(min_mz=mzmin, max_mz=mzmax, min_rt=exp.getMinRT(), max_rt=exp.getMaxRT())
        chrom_data = pd.DataFrame({"mz_i":mz_i, "rt":chrom_data[0], "mz":chrom_data[1], "int":chrom_data[2]})
        scan_dfs.append(chrom_data)
    return(pd.concat(scan_dfs, ignore_index=True))
def get_multichrom_mzml_pymzml(file, mzs, ppm):
    run = pymzml.run.Reader(file, build_index_from_scratch=True)
    scan_dfs = []
    for mz_i in mzs:
        mzmin, mzmax = pmppm(mz_i, ppm)
        for spectrum in run:
            if(spectrum.ms_level==1):
                rt_val = spectrum.scan_time_in_minutes()
                mz_vals = spectrum.mz
                int_vals = spectrum.i
                bet_idxs = (mzmin < mz_vals) & (mz_vals < mzmax)
            
                if(sum(bet_idxs)>0):
                    df_scan = pd.DataFrame({'mz_i':mz_i, 'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
                    scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
def get_multichrom_mz5(file, mzs, ppm):
    mz5_file = h5py.File(file, 'r')
    bounds_df = pd.DataFrame({"lower":np.concatenate(([0], mz5_file["SpectrumIndex"][...][:-1])),
                              "upper":mz5_file["SpectrumIndex"][...]})
    scan_dfs = []
    for mz_i in mzs:
        mzmin, mzmax = pmppm(mz_i, ppm)
        has_precursor = [len(item)>0 for item in mz5_file["SpectrumMetaData"]["precursors"]]
        for index, row in bounds_df.iterrows():
            if(not(has_precursor[index])):
                scan_df = pd.DataFrame({
                    "mz_i": mz_i,
                    "rt": mz5_file["ChomatogramTime"][index],
                    "mz": np.cumsum(mz5_file["SpectrumMZ"][row["lower"]:row["upper"]]),
                    "int": mz5_file["SpectrumIntensity"][row["lower"]:row["upper"]]
                })
                bet_df = scan_df[(scan_df["mz"]>mzmin) & (scan_df["mz"]<mzmax)]
                scan_dfs.append(bet_df)
    file_df = pd.concat(scan_dfs, ignore_index=True)
    mz5_file.close()
    return(file_df)
def get_multichrom_mza(file, mzs, ppm):
    mza=mzapy.MZA(file)
    scan_dfs = []
    for mz_i in mzs:
        mzmin, mzmax = pmppm(mz_i, ppm)
        xic_rt, xic_int = mza.collect_xic_arrays_by_mz(mzmin, mzmax)
        scan_dfs.append(pd.DataFrame({"mz_i":mz_i, "rt":xic_rt, "int":xic_int}))
    file_df = pd.concat(scan_dfs, ignore_index=True)
    mza.close()
    return(file_df)
def get_multichrom_mzdb(file, mzs, ppm):
    connection = sqlite3.connect(file)
    cursor = connection.cursor()
    bb_query = """
        SELECT bounding_box.id, begin_mz, end_mz, first_spectrum_id, time AS first_spectrum_time, data
        FROM run_slice, bounding_box, spectrum
        WHERE bounding_box.run_slice_id = run_slice.id
        AND bounding_box.first_spectrum_id = spectrum.id
        AND spectrum.ms_level = 1
        AND begin_mz = ?
        """
    all_chroms = []
    for mz_i in mzs:
        mzmin, mzmax = pmppm(mz_i, ppm)
        spec_bb_query = "SELECT MAX(begin_mz) FROM run_slice WHERE begin_mz < ?"
        bb_id_for_chrom = cursor.execute(spec_bb_query, (mzmin,)).fetchone()[0]
        scanid_rt_pd = pd.read_sql("SELECT id AS scan_id, time AS rt, ms_level FROM spectrum", connection)
        bb_dataframe = pd.read_sql(bb_query, connection, params=(bb_id_for_chrom,))
        unpacked_bb_list = [unpack_raw_bb(bb_data) for bb_data in bb_dataframe["data"]]
        bb_chrom = pd.concat(unpacked_bb_list).merge(scanid_rt_pd)
        bb_chrom["rt"] /= 60
        chrom_data = bb_chrom[(mzmin < bb_chrom["mz"]) & (bb_chrom["mz"] < mzmax)].copy()
        chrom_data["mz_i"] = mz_i
        all_chroms.append(chrom_data)
    connection.close()
    
    return(pd.concat(all_chroms, ignore_index=True))
def get_multichrom_duckdb_loop(file, mzs, ppm):
    conn = duckdb.connect(file)
    all_chroms = []
    for mz_i in mzs:
        mz_min, mz_max = pmppm(mz_i, ppm)
        query = "SELECT mz, int, rt FROM MS1 WHERE mz BETWEEN ? AND ?"
        query_data = conn.execute(query, (mz_min, mz_max)).fetchdf()
        query_data["mz_i"] = mz_i
        all_chroms.append(query_data)
    conn.close()
    return(pd.concat(all_chroms, ignore_index=True))
def get_multichrom_sqlite_loop(file, mzs, ppm):
    conn = sqlite3.connect(file)
    all_chroms = []
    for mz_i in mzs:
        mz_min, mz_max = pmppm(mz_i, ppm)
        query = "SELECT mz, int, rt FROM MS1 WHERE mz BETWEEN ? AND ?"
        query_data = pd.read_sql_query(query, conn, params=(mz_min, mz_max))
        query_data["mz_i"] = mz_i
        all_chroms.append(query_data)
    conn.close()
    return(pd.concat(all_chroms, ignore_index=True))
def get_multichrom_parquet_loop(pqds_dir, mzs, ppm):
    scan_data =[]
    for mz_i in mzs:
        mz_min, mz_max = pmppm(mz_i, ppm)
        dataset = ds.dataset(f"{pqds_dir}/MS1", format="parquet")
        bet_df = dataset.to_table(filter=((ds.field("mslevel")=="MS1") & (ds.field('mz') >= mz_min) & (ds.field('mz') <= mz_max))).to_pandas()
        bet_df["mz_i"] = mz_i
        scan_data.append(bet_df)
    return pd.concat(scan_data, ignore_index=True)

function_list = [
#    ("pyteomics", "mzml_pyteomics", ".mzML"),
#    ("pyopenms", "mzml_pyopenms", ".mzML"),
#    ("pyopenms_2DPeak", "mzml_pyopenms_2DPeak", ".mzML"),
#    ("pymzml", "mzml_pymzml", ".mzML"),
#    ("mzMLb", "mzmlb", ".mzMLb"),
#    ("mzDB", "mzdb", ".raw.mzDB"),
#    ("MZA", "mza", ".mza"),
#    ("mz5", "mz5", ".mz5"),
    ("SQLite", "sqlite_loop", ".sqlite"),
    ("DuckDB", "duckdb_loop", ".duckdb"),
    ("Parquet", "parquet_loop", "_pqds")
]

def time_multichrom(fun_suffix, file_ending, mz_list, verbose=True):
    if(verbose):
        start_time = time.time()
        print(f"Running {fun_suffix}", end="\r")
    rep_function = f"get_multichrom_{fun_suffix}('E:/mzsql/MTBLS10066/20220923_LEAP-POS_QC04{file_ending}', {mz_list}, 5)"
    time_vals = timeit.repeat(rep_function, globals=globals(), number=1, repeat=3)
    if(verbose):
        elapsed = time.time() - start_time
        print(f"Running {fun_suffix}... ({elapsed:.2f}s)")
    return(time_vals)

random.seed(123)
for n_chroms in [1, 3, 10, 30, 100]:
    print(n_chroms)
    mz_list = random.sample(top_masses, n_chroms)
    for name, suffix, filetype in function_list:
        time_taken = time_multichrom(suffix, filetype, mz_list, verbose=True)
        time_df = pd.DataFrame({"n_chrom":n_chroms, "method":name, "time":time_taken})
        time_df.to_csv('data/multichrom_times.csv', mode='a', header=False, index=False)
