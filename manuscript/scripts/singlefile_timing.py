
# Setup things
import numpy as np
import pandas as pd
import sqlite3
import os
from mzsql import *

# Create databases from mzML files if they don't already exist
if(not(os.path.exists("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.sqlite"))):
    turn_mzml_sqlite("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.mzML", "E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.sqlite")
    turn_mzml_duckdb("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.mzML", "E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.duckdb")
    turn_mzml_parquet("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.mzML", "E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01_pqds")

# Identify the bits of information to be extracted
conn = sqlite3.connect("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.sqlite")
# Get 10 random spectra
rand_scan_query = """
SELECT id 
FROM (SELECT DISTINCT id FROM {0}) 
ORDER BY RANDOM() 
LIMIT 10
"""
rand_MS1_scans = pd.read_sql_query(rand_scan_query.format("MS1"), conn)["id"].tolist()
rand_MS2_scans = pd.read_sql_query(rand_scan_query.format("MS2"), conn)["id"].tolist()
# Get the top 10ish masses by intensity and their retention times
top_int_df = pd.read_sql_query("SELECT * FROM MS1 ORDER BY int DESC LIMIT 30000", conn)
top_masses = []
top_rts = []
for i in range(10):
    top_rts.append(top_int_df["rt"].iloc[0])
    top_mz = top_int_df["mz"].iloc[0]
    top_masses.append(top_mz)
    mzmin, mzmax = pmppm(top_mz, 10)
    rtmin, rtmax = top_int_df["rt"].iloc[0]+(-1, 1)
    top_int_df = top_int_df[((top_int_df["mz"]<mzmin) | (top_int_df["mz"]>mzmax)) & ((top_int_df["rt"]<rtmin) | (top_int_df["rt"]>rtmax))]
# Get the top 10ish fragments and precursors
top_fragments = []
top_frag_df = pd.read_sql_query("SELECT * FROM MS2 ORDER BY int DESC LIMIT 1000", conn)
for i in range(10):
    top_frag = top_frag_df["fragmz"].iloc[0]
    top_fragments.append(top_frag)
    mzmin, mzmax = pmppm(top_frag, 10)
    top_frag_df = top_frag_df[((top_frag_df["fragmz"]<mzmin) | (top_frag_df["fragmz"]>mzmax))]

top_precursors = []
top_frag_df = pd.read_sql_query("SELECT * FROM MS2 ORDER BY int DESC LIMIT 1000", conn)
for i in range(10):
    top_frag = top_frag_df["premz"].iloc[0]
    top_precursors.append(top_frag)
    mzmin, mzmax = pmppm(top_frag, 10)
    top_frag_df = top_frag_df[((top_frag_df["premz"]<mzmin) | (top_frag_df["premz"]>mzmax))]

# Create a list of functions to loop over
function_list = [
    ("pyteomics", "mzml_pyteomics", ".mzML"),
    ("pyopenms", "mzml_pyopenms", ".mzML"),
    ("pyopenms_2DPeak", "mzml_pyopenms_2DPeak", ".mzML"),
    ("pymzml", "mzml_pymzml", ".mzML"),
    ("mzMLb", "mzmlb", ".mzMLb"),
    ("mzDB", "mzdb", ".raw.mzDB"),
    ("MZA", "mza", ".mza"),
    ("mz5", "mz5", ".mz5"),
    ("SQLite", "sqlite", ".sqlite"),
    ("DuckDB", "duckdb", ".duckdb"),
    ("Parquet", "parquet", "_pqds"),
    ("MZTree", "mztree", "http://127.0.0.1:4568"),
    ("mzMD", "mzMD", "http://127.0.0.1:4567")
]

# Create a CSV to store the data in
init_df = pd.DataFrame(columns=["query", "target", "method","time"])
init_df.to_csv('data/singlefile_times.csv', index=False)

# And perform the timings!
for target_spec in rand_MS1_scans:
    print(target_spec)
    for name, suffix, filetype in [function_list[i] for i in [0, 1, 3, 4,5,6,8,9,10,11]]:
        time_taken = time_spec(suffix, filetype, target_spec, verbose=True)
        time_df = pd.DataFrame({"query":["ms1_spec"], "target":[target_spec], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

for target_mz in top_masses:
    print(target_mz)
    for name, suffix, filetype in function_list:
        time_taken = time_chrom(suffix, filetype, target_mz, 10, verbose=True)
        time_df = pd.DataFrame({"query":["chrom"], "target":[target_mz], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

for target_rt in top_rts:
    print(target_rt)
    for name, suffix, filetype in function_list:
        time_taken = time_rtrange(suffix, filetype, target_rt, verbose=True)
        time_df = pd.DataFrame({"query":["rtrange"], "target":[target_rt], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

for target_spec in rand_MS2_scans:
    print(target_spec)
    for name, suffix, filetype in [function_list[i] for i in [0, 1, 3, 4,5,6,8,9,10,11]]:
        time_taken = time_spec(suffix, filetype, target_spec, verbose=True)
        time_df = pd.DataFrame({"query":["ms2_spec"], "target":[target_spec], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

for target_pre in top_precursors:
    print(target_pre)
    for name, suffix, filetype in [function_list[i] for i in [0, 1, 3, 4,5,6,9,10,11]]:
        time_taken = time_premz(suffix, filetype, target_pre, 10, verbose=True)
        time_df = pd.DataFrame({"query":["premz"], "target":[target_pre], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

for target_frag in top_fragments:
    print(target_frag)
    for name, suffix, filetype in [function_list[i] for i in [0, 1, 3, 4,5,6,9,10,11]]:
        time_taken = time_fragmz(suffix, filetype, target_frag, 10, verbose=True)
        time_df = pd.DataFrame({"query":["fragmz"], "target":[target_frag], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

print("Success!")

file_sizes = {}
for name, suffix, filetype in function_list:
    if(name == "mzMD"):
        file_sizes[suffix] = os.path.getsize("E:/mzsql/MTBLS10066/result.mzMD")
    elif(name == "MZTree"):
        file_sizes[suffix] = os.path.getsize("E:/mzsql/MTBLS10066/01-30-2025_10-35-06.mzTree")
    elif(name == "Parquet"):
        dirsize = sum(
            os.path.getsize(os.path.join(root, f)) 
            for root, _, files in os.walk("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01_pqds") 
            for f in files
        )
        file_sizes[suffix] = dirsize
    else:
        file_sizes[suffix] = os.path.getsize(f"E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01{filetype}")

pd.DataFrame(file_sizes.items(), columns=["method", "file_size"]).to_csv("data/file_sizes.csv", index=False)
