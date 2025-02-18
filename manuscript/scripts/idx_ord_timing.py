
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

if(not(os.path.exists("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01_ord.sqlite"))):
    turn_mzml_sqlite("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.mzML", "E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01_ord.sqlite", ordered="mz")
    turn_mzml_duckdb("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.mzML", "E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01_ord.duckdb", ordered="mz")
    turn_mzml_parquet("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.mzML", "E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01_pqds_ord", ordered="mz")

conn = sqlite3.connect("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01.sqlite")
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

function_list = [
    ("SQLite", "sqlite", ".sqlite"),
    ("SQLite (indexed)", "sqlite", "_ord.sqlite"),
    ("DuckDB", "duckdb", ".duckdb"),
    ("DuckDB (ordered)", "duckdb", "_ord.duckdb"),
    ("Parquet", "parquet", "_pqds"),
    ("Parquet (ordered)", "parquet", "_pqds_ord")
]

for target_mz in top_masses:
    print(target_mz)
    for name, suffix, filetype in function_list:
        time_taken = time_chrom(suffix, filetype, target_mz, 10, verbose=True)
        time_df = pd.DataFrame({"query":["chrom"], "target":[target_mz], "method":[name], "time":time_taken})
        time_df.to_csv('data/idx_ord_times.csv', mode='a', header=False, index=False)

file_sizes = {}
for name, suffix, filetype in function_list:
    if(suffix == "parquet"):
        dirsize = sum(
            os.path.getsize(os.path.join(root, f)) 
            for root, _, files in os.walk("E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01_pqds") 
            for f in files
        )
        file_sizes[name] = dirsize
    else:
        file_sizes[name] = os.path.getsize(f"E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01{filetype}")

pd.DataFrame(file_sizes.items(), columns=["method", "file_size"]).to_csv("data/idx_sizes.csv", index=False)
