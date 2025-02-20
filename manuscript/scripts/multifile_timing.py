import numpy as np
import pandas as pd
import seaborn as sns
import sqlite3
import os
import glob
import random
from mzsql import *

random.seed(123)
multifile_subset = random.sample(glob.glob("E:/mzsql/MTBLS10066/202*.mzML"), 100)
multifile_subset = [f.replace("\\", "/") for f in multifile_subset]

# Create individual sqlite, duckdb, and parquet files
# for file_i in multifile_subset:
#     turn_mzml_sqlite(file_i, outfile=file_i.replace(".mzML", ".sqlite"), ordered="mz")
#     turn_mzml_duckdb(file_i, outfile=file_i.replace(".mzML", ".duckdb"), ordered="mz")
#     turn_mzml_parquet(file_i, outdir=file_i.replace(".mzML", "_pqds"), ordered="mz")

# Create multifile dbs
# turn_mzml_sqlite(multifile_subset[0], outfile="E:/mzsql/MTBLS10066/multifile_1.sqlite", ordered="mz")
# turn_mzml_duckdb(multifile_subset[0], outfile="E:/mzsql/MTBLS10066/multifile_1.duckdb", ordered="mz")
# turn_mzml_parquet(multifile_subset[0], outdir="E:/mzsql/MTBLS10066/multifile_1_pqds", ordered="mz")
# turn_mzml_sqlite(multifile_subset[:10], outfile="E:/mzsql/MTBLS10066/multifile_10.sqlite", ordered="mz")
# turn_mzml_duckdb(multifile_subset[:10], outfile="E:/mzsql/MTBLS10066/multifile_10.duckdb", ordered="mz")
# turn_mzml_parquet(multifile_subset[:10], outdir="E:/mzsql/MTBLS10066/multifile_10_pqds", ordered="mz")
# turn_mzml_sqlite(multifile_subset, outfile="E:/mzsql/MTBLS10066/multifile_100.sqlite", ordered="mz")
# turn_mzml_duckdb(multifile_subset, outfile="E:/mzsql/MTBLS10066/multifile_100.duckdb", ordered="mz")
# turn_mzml_parquet(multifile_subset, outdir="E:/mzsql/MTBLS10066/multifile_100_pqds", ordered="mz")

# Get the largest masses to query
conn = sqlite3.connect("E:/mzsql/MTBLS10066/multifile_1.sqlite")
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
conn.close()

singlefile_timings = []
for mz_i in top_masses:
    print(mz_i)
    for n_files in [1, 10, 100]:
        print(n_files)
        this_multifile_subset = multifile_subset[:n_files]
        total_time = {"duckdb":0, "sqlite":0, "parquet":0}
        for mzml_file in this_multifile_subset:
            for db_name, db_ending in [("duckdb", ".duckdb"), ("sqlite", ".sqlite"), ("parquet", "_pqds")]:
                db_path = mzml_file.replace(".mzML", db_ending)
                rep_function = f"get_chrom_{db_name}('{db_path}', {mz_i}, 5)"
                time_val = timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
                total_time[db_name] += time_val[0]
        for db_name, time_val in total_time.items():
            time_data = pd.DataFrame([{
                "method": db_name,
                "mz_target": mz_i,
                "n_files": n_files,
                "type": "singlefile",
                "time": time_val
            }])
            singlefile_timings.append(time_data)

multifile_timings = []
for mz_i in top_masses:
    print(mz_i)
    for n_files in [1, 10, 100]:
        print(n_files)
        for db_name, db_ending in [("duckdb", ".duckdb"), ("sqlite", ".sqlite"), ("parquet", "_pqds")]:
            db_path = f"E:/mzsql/MTBLS10066/multifile_{n_files}{db_ending}"
            print(db_path)
            rep_function = f"get_chrom_{db_name}('{db_path}', {mz_i}, 5)"
            time_vals = timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
            time_data = pd.DataFrame({"method":db_name, "mz_target":mz_i, "n_files":n_files, "type":"multifile",
                                      "mzml_file":os.path.basename(mzml_file), 
                                      "time":time_vals})
            multifile_timings.append(time_data)

singlefile_df = pd.concat(singlefile_timings, ignore_index=True)
multifile_df = pd.concat(multifile_timings, ignore_index=True)
pd.concat([singlefile_df, multifile_df], ignore_index=True).to_csv('data/multifile_times.csv', mode='w', index=False)