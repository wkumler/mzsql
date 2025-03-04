
# Setup things
import numpy as np
import pandas as pd
import sqlite3
import os
import glob
import random
import timeit
from mzsql import *

random.seed(123)
basename=random.sample(glob.glob("E:/mzsql/MTBLS10066/*.mzML"), 1)[0].replace(".mzML", "").replace("\\", "/")
print(basename)
# 20220923_LEAP-POS_QC04

# Create databases from mzML files if they don't already exist
if(not(os.path.exists(f"{basename}.sqlite"))):
    turn_mzml_sqlite(f"{basename}.mzML", f"{basename}.sqlite")
    turn_mzml_duckdb(f"{basename}.mzML", f"{basename}.duckdb")
    turn_mzml_parquet(f"{basename}.mzML", f"{basename}_pqds")

# # Identify the bits of information to be extracted
conn = sqlite3.connect(f"{basename}.sqlite")
# # Get 10 random spectra
# rand_scan_query = """
# SELECT id 
# FROM (SELECT DISTINCT id FROM {0}) 
# ORDER BY RANDOM() 
# LIMIT 10
# """
# rand_MS1_scans = pd.read_sql_query(rand_scan_query.format("MS1"), conn)["id"].tolist()
# rand_MS2_scans = pd.read_sql_query(rand_scan_query.format("MS2"), conn)["id"].tolist()
# print(rand_MS1_scans)
# print(rand_MS2_scans)
# Manually specify scans here since SQLite doesn't handle seed setting very well:
rand_MS1_scans = (9044, 187, 7532, 2396, 10409, 87, 5259, 69, 5157, 7988)
rand_MS2_scans = (3778, 3815, 5685, 6859, 4223, 3372, 6948, 6354, 3140, 912)

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

conn.close()


# Define timing functions
def time_chrom(basename, fun_suffix, file_ending, target_mz, ppm_acc, verbose=True):
    if(fun_suffix in ["mztree", "mzMD"]):
        rep_function = f"get_chrom_{fun_suffix}('{file_ending}', {target_mz}, {ppm_acc})"
    else:
        rep_function = f"get_chrom_{fun_suffix}('{basename}{file_ending}', {target_mz}, {ppm_acc})"
    return timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
def time_spec(basename, fun_suffix, file_ending, spec_id, verbose=True):
    rep_function = f"get_spec_{fun_suffix}('{basename}{file_ending}', {spec_id})"
    return timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
def time_rtrange(basename, fun_suffix, file_ending, rt, verbose=True):
    if(fun_suffix in ["mztree", "mzMD"]):
        rep_function = f"get_rtrange_{fun_suffix}('{file_ending}', {rt-0.2}, {rt+0.2})"
    else:
        rep_function = f"get_rtrange_{fun_suffix}('{basename}{file_ending}', {rt-0.2}, {rt+0.2})"
    return timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
def time_premz(basename, fun_suffix, file_ending, target_premz, ppm_acc, verbose=True):
    rep_function = f"get_MS2premz_{fun_suffix}('{basename}{file_ending}', {target_premz}, {ppm_acc})"
    return timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
def time_fragmz(basename, fun_suffix, file_ending, target_fragmz, ppm_acc, verbose=True):
    rep_function = f"get_MS2fragmz_{fun_suffix}('{basename}{file_ending}', {target_fragmz}, {ppm_acc})"
    return timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)

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
    ("MZTree", "mztree", "http://127.0.0.1:4568"),
    ("mzMD", "mzMD", "http://127.0.0.1:4567"),
    ("SQLite", "sqlite", ".sqlite"),
    ("DuckDB", "duckdb", ".duckdb"),
    ("Parquet", "parquet", "_pqds")
]

# Create a CSV to store the data in
init_df = pd.DataFrame(columns=["query", "target", "method","time"])
init_df.to_csv('data/singlefile_times.csv', index=False)

# And perform the timings!
for target_spec in rand_MS1_scans:
    print(target_spec)
    for name, suffix, filetype in [function_list[i] for i in [0, 1, 3, 4,5,6,7,10,11,12]]:
        time_taken = time_spec(basename, suffix, filetype, target_spec, verbose=True)
        time_df = pd.DataFrame({"query":["ms1_spec"], "target":[target_spec], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

for target_mz in top_masses:
    print(target_mz)
    for name, suffix, filetype in function_list:
        time_taken = time_chrom(basename, suffix, filetype, target_mz, 10, verbose=True)
        time_df = pd.DataFrame({"query":["chrom"], "target":[target_mz], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

for target_rt in top_rts:
    print(target_rt)
    for name, suffix, filetype in function_list:
        time_taken = time_rtrange(basename, suffix, filetype, target_rt, verbose=True)
        time_df = pd.DataFrame({"query":["rtrange"], "target":[target_rt], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

for target_spec in rand_MS2_scans:
    print(target_spec)
    for name, suffix, filetype in [function_list[i] for i in [0, 1, 3, 4, 5, 6, 7, 10, 11, 12]]:
        time_taken = time_spec(basename, suffix, filetype, target_spec, verbose=True)
        time_df = pd.DataFrame({"query":["ms2_spec"], "target":[target_spec], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

for target_pre in top_precursors:
    print(target_pre)
    for name, suffix, filetype in [function_list[i] for i in [0, 1, 3, 4, 5, 6, 10, 11, 12]]:
        time_taken = time_premz(basename, suffix, filetype, target_pre, 10, verbose=True)
        time_df = pd.DataFrame({"query":["premz"], "target":[target_pre], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

for target_frag in top_fragments:
    print(target_frag)
    for name, suffix, filetype in [function_list[i] for i in [0, 1, 3, 4, 5, 6, 10, 11, 12]]:
        time_taken = time_fragmz(basename, suffix, filetype, target_frag, 10, verbose=True)
        time_df = pd.DataFrame({"query":["fragmz"], "target":[target_frag], "method":[name], "time":time_taken})
        time_df.to_csv('data/singlefile_times.csv', mode='a', header=False, index=False)

file_sizes = {}
for name, suffix, filetype in function_list:
    if(name == "mzMD"):
        file_sizes[name] = os.path.getsize(f"{basename}.mzMD")
    elif(name == "MZTree"):
        total_size = sum([
            os.path.getsize(f"{basename}.mzTree"),
            os.path.getsize(f"{basename}.mzTree-points")
        ])
        file_sizes[name] = total_size
    elif(name == "Parquet"):
        dirsize = sum(
            os.path.getsize(os.path.join(root, f)) 
            for root, _, files in os.walk(f"{basename}_pqds") 
            for f in files
        )
        file_sizes[name] = dirsize
    else:
        file_sizes[name] = os.path.getsize(f"{basename}{filetype}")

pd.DataFrame(file_sizes.items(), columns=["method", "file_size"]).to_csv("data/file_sizes.csv", index=False)

print("Success!")
