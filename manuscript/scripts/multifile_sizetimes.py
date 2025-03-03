import numpy as np
import pandas as pd
import sqlite3
import os
import glob
import random
import timeit
from mzsql import *

random.seed(123)
multifile_subset = random.sample(glob.glob("E:/mzsql/MTBLS10066/*.mzML"), 100)
multifile_subset = [f.replace("\\", "/") for f in multifile_subset]

file_sequence = [1, 3, 10, 30, 100]

for n_files in file_sequence:
    for db_sort in ["ordered", "unordered"]:
        ordering = "mz" if db_sort == "ordered" else None
        for db_name in ["duckdb", "sqlite", "parquet"]:
            db_ending = {"duckdb": ".duckdb", "sqlite": ".sqlite", "parquet": "_pqds"}
            file_ending = "_ord" if db_sort == "ordered" else ""
            file_ending += db_ending[db_name]
            for db_type in ["consolidated", "multifile"]:
                this_multifile_subset = multifile_subset[:n_files]
                if db_type == "multifile":
                    if n_files < max(file_sequence):
                        continue
                    for mzml_i in this_multifile_subset:
                        output_db = mzml_i.replace(".mzML", file_ending)
                        print(output_db)
                        combined_cmd = f"turn_mzml_{db_name}(mzml_i, output_db, ordered=ordering)"
                        #eval(combined_cmd)
                else:
                    output_db = f"E:/mzsql/MTBLS10066/consolidated_{n_files}{file_ending}"
                    print(output_db)
                    combined_cmd = f"turn_mzml_{db_name}(this_multifile_subset, output_db, ordered=ordering)"
                    #eval(combined_cmd)

conn = sqlite3.connect(f"E:/mzsql/MTBLS10066/consolidated_{max(file_sequence)}.sqlite")
top_int_df = pd.read_sql_query("SELECT * FROM MS1 ORDER BY int DESC LIMIT 30000", conn)
top_masses = []
top_rts = []
for i in range(10):
    top_rts.append(top_int_df["rt"].iloc[0])
    top_mz = top_int_df["mz"].iloc[0]
    top_masses.append(top_mz)
    mzmin, mzmax = pmppm(top_mz, 10)
    rtmin, rtmax = top_int_df["rt"].iloc[0]+(-1, 1)
    top_int_df = top_int_df[((top_int_df["mz"]<mzmin) | (top_int_df["mz"]>mzmax))]
conn.close()

multifile_timings = []
for mz_i in top_masses:
    for n_files in file_sequence:
        for db_sort in ["ordered", "unordered"]:
            ordering = "mz" if db_sort == "ordered" else None
            for db_name in ["duckdb", "sqlite", "parquet"]:
                db_ending = {"duckdb": ".duckdb", "sqlite": ".sqlite", "parquet": "_pqds"}
                file_ending = "_ord" if db_sort == "ordered" else ""
                file_ending += db_ending[db_name]
                for db_type in ["consolidated", "multifile"]:
                    this_multifile_subset = multifile_subset[:n_files]
                    if db_type == "multifile":
                        for mzml_i in this_multifile_subset:
                            output_db = mzml_i.replace(".mzML", file_ending)
                            rep_function = f"get_chrom_{db_name}('{output_db}', {mz_i}, 5)"
                            time_vals = timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
                            time_data = pd.DataFrame({"method":db_name, "mz_target":mz_i, "n_files":n_files, 
                                                      "type":db_type, "db_sort":db_sort, "time":time_vals})
                            multifile_timings.append(time_data)
                    else:
                        output_db = f"E:/mzsql/MTBLS10066/consolidated_{n_files}{file_ending}"
                        rep_function = f"get_chrom_{db_name}('{output_db}', {mz_i}, 5)"
                        time_vals = timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
                        time_data = pd.DataFrame({"method":db_name, "mz_target":mz_i, "n_files":n_files, 
                                                  "type":db_type, "db_sort":db_sort,
                                                  "time":time_vals})
                        multifile_timings.append(time_data)

def get_dbsize(basename, db_ending):
    if db_ending == "_pqds":
        return sum(
            os.path.getsize(os.path.join(root, f))
            for root, _, files in os.walk(basename)
            for f in files
        )
    return os.path.getsize(basename)

file_sizes = []
for n_files in file_sequence:
    for db_sort in ["ordered", "unordered"]:
        sort_ending = "_ord" if db_sort == "ordered" else ""
        for db_name in ["duckdb", "sqlite", "parquet"]:
            db_ending = {"duckdb": ".duckdb", "sqlite": ".sqlite", "parquet": "_pqds"}[db_name]
            for db_type in ["consolidated", "multifile"]:
                if db_type == "multifile":
                    appendage = sort_ending+db_ending
                    this_multifile_subset = [i.replace(".mzML", appendage) for i in multifile_subset[:n_files]]
                    db_size = sum([get_dbsize(file_i, db_ending) for file_i in this_multifile_subset])
                else:
                    db_size = get_dbsize(f"E:/mzsql/MTBLS10066/consolidated_{n_files}{sort_ending}{db_ending}", db_ending)
                size_data = pd.DataFrame({"method": db_name, "n_files": n_files, "type": db_type,
                                          "db_sort": db_sort, "size": [db_size]})
                file_sizes.append(size_data)

all_timings = pd.concat(multifile_timings, ignore_index=True)
all_sizes = pd.concat(file_sizes, ignore_index=True)
all_sizes.merge(all_timings).to_csv("data/multifile_timings.csv", index=False)

