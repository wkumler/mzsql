
import numpy as np
import pandas as pd
import os
import duckdb
import pyteomics.mzml
from .helpers import pmppm

def turn_mzml_duckdb(files, outfile, ordered=None):
    conn = duckdb.connect(outfile)
    conn.execute("DROP TABLE IF EXISTS MS1")
    conn.execute("DROP TABLE IF EXISTS MS2")
    empty_ms1 = pd.DataFrame(columns=["filename", "id", "mz", "int", "rt"]).astype(
        {"filename": "string", "id": "int64", "mz": "float64", "int": "float64", "rt": "float64"}
    )
    conn.execute("CREATE TABLE MS1 AS SELECT * FROM empty_ms1")
    empty_ms2 = pd.DataFrame(columns=["filename", "id", "premz", "fragmz", "int", "rt"]).astype(
        {"filename": "string", "id": "int64", "premz": "float64", "fragmz": "float64", "int": "float64", "rt": "float64"}
    )
    conn.execute("CREATE TABLE MS2 AS SELECT * FROM empty_ms2")

    if isinstance(files, str):
        files = [files]
    for file in files:
        if not os.path.exists(file):
            raise FileNotFoundError(f"Input mzML file '{file}' does not exist.")
        MS1_dfs = []
        MS2_dfs = []
        for spectrum in pyteomics.mzml.MzML(file):
            if spectrum['ms level'] == 1:
                idx = int(spectrum['id'].split("scan=")[-1].split()[0])
                mz_vals = spectrum['m/z array']
                int_vals = spectrum['intensity array']
                rt_val = spectrum['scanList']['scan'][0]['scan start time']
                df_scan = pd.DataFrame({'filename': os.path.basename(file), 'id': idx, 
                                        'mz': mz_vals, 'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
                MS1_dfs.append(df_scan)
            if spectrum["ms level"] == 2:
                idx = int(spectrum['id'].split("scan=")[-1].split()[0])
                premz_val = spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
                mz_vals = spectrum['m/z array']
                int_vals = spectrum['intensity array']
                rt_val = spectrum['scanList']['scan'][0]['scan start time']
                df_scan = pd.DataFrame({'filename': os.path.basename(file), 'id': idx, 'premz': premz_val, 
                                        'fragmz': mz_vals, 'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
                MS2_dfs.append(df_scan)
    
        all_MS1 = pd.concat(MS1_dfs, ignore_index=True)
        all_MS2 = pd.concat(MS2_dfs, ignore_index=True)
        if ordered is not None:
            if ordered == "rt":
                all_MS1.sort_values("rt", inplace=True)
                all_MS2.sort_values("rt", inplace=True)
            if ordered == "mz":
                all_MS1.sort_values("mz", inplace=True)
            if ordered in ["fragmz", "premz"]:
                all_MS1.sort_values(ordered, inplace=True)
        conn.execute("INSERT INTO MS1 SELECT * FROM all_MS1")
        conn.execute("INSERT INTO MS2 SELECT * FROM all_MS2")
    conn.close()

    return outfile

def get_chrom_duckdb(file, mz, ppm):
    conn = duckdb.connect(file)
    mz_min, mz_max = pmppm(mz, ppm)
    query = "SELECT filename, mz, int, rt FROM MS1 WHERE mz BETWEEN ? AND ?"
    query_data = conn.execute(query, (mz_min, mz_max)).fetchdf()    
    conn.close()
    return query_data

def get_spec_duckdb(file, spectrum_idx):
    conn = duckdb.connect(file)
    query = "SELECT id, mz, int FROM MS1 WHERE id = ? UNION ALL SELECT id, fragmz AS mz, int FROM MS2 WHERE id = ?"
    spectrum_data = conn.execute(query, (spectrum_idx, spectrum_idx)).fetchdf()
    conn.close()
    return spectrum_data

def get_rtrange_duckdb(file, rtstart, rtend):
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS1 WHERE rt BETWEEN ? AND ?"
    rt_range_data = conn.execute(query, (rtstart, rtend)).fetchdf()
    conn.close()
    return rt_range_data

def get_MS2fragmz_duckdb(file, fragment_mz, ppm_acc):
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS2 WHERE fragmz BETWEEN ? AND ?"
    spectrum_data = conn.execute(query, pmppm(fragment_mz, ppm_acc)).fetchdf()
    conn.close()
    return(spectrum_data)

def get_MS2premz_duckdb(file, precursor_mz, ppm_acc):
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS2 WHERE premz BETWEEN ? AND ?"
    spectrum_data = conn.execute(query, pmppm(precursor_mz, ppm_acc)).fetchdf()
    conn.close()
    return(spectrum_data)

def get_MS2nloss_duckdb(file, nloss_mz, ppm_acc):
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS2 WHERE (premz-fragmz) BETWEEN ? AND ?"
    spectrum_data = conn.execute(query, pmppm(nloss_mz, ppm_acc)).fetchdf()
    conn.close()
    return(spectrum_data)
