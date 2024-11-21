
import numpy as np
import pandas as pd
import duckdb
import pyteomics.mzml
from .helpers import pmppm

# def turn_mzml_duckdb(file, outfile):
#     conn = duckdb.connect(outfile)
#     conn.execute("DROP TABLE IF EXISTS MS1")
#     for spectrum in pyteomics.mzml.MzML(file):
#         if spectrum['ms level'] == 1:
#             print(spectrum["index"])
#             idx = int(spectrum['id'].split("scan=")[-1].split()[0])
#             mz_vals=spectrum['m/z array']
#             int_vals = spectrum['intensity array']
#             rt_val = spectrum['scanList']['scan'][0]['scan start time']
#             df_scan = pd.DataFrame({'id':idx,'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(mz_vals)})
#             df_scan.to_sql("MS1", conn, if_exists="append", index=False)
#     conn.close()
#     return(outfile)

def turn_mzml_duckdb(file, outfile, ordered=False):
    conn = duckdb.connect(outfile)
    conn.execute("DROP TABLE IF EXISTS MS1")

    scan_dfs = []
    for spectrum in pyteomics.mzml.MzML(file):
        if spectrum['ms level'] == 1:
            idx = int(spectrum['id'].split("scan=")[-1].split()[0])
            mz_vals = spectrum['m/z array']
            int_vals = spectrum['intensity array']
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            df_scan = pd.DataFrame({'id': idx, 'mz': mz_vals, 'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
            scan_dfs.append(df_scan)

    all_pds = pd.concat(scan_dfs, ignore_index=True)
    if ordered:
        all_pds.sort_values(by="mz", inplace=True)

    conn.execute("CREATE TABLE MS1 AS SELECT * FROM all_pds")
    conn.close()

    return outfile

def get_chrom_duckdb(file, mz, ppm):
    conn = duckdb.connect(file)
    mz_min, mz_max = pmppm(mz, ppm)
    query = "SELECT mz, int, rt FROM MS1 WHERE mz BETWEEN ? AND ?"
    query_data = conn.execute(query, (mz_min, mz_max)).fetchdf()    
    conn.close()
    return query_data

def get_spec_duckdb(file, spectrum_idx):
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS1 WHERE id = ?"
    spectrum_data = conn.execute(query, (spectrum_idx,)).fetchdf()
    conn.close()
    return spectrum_data

def get_rtrange_duckdb(file, rtstart, rtend):
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS1 WHERE rt BETWEEN ? AND ?"
    rt_range_data = conn.execute(query, (rtstart, rtend)).fetchdf()
    conn.close()
    return rt_range_data