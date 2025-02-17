
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mzsql import *
import sqlite3
import os
import random

random.seed(123)
files = ["E:/mzsql/MTBLS10066/"+f for f in os.listdir("E:/mzsql/MTBLS10066/") if f.endswith(".mzML")]
msfiles = random.sample(files, 10)

if(not(os.path.exists("E:/mzsql/MTBLS10066/20220923_LEAP-POS_QC04.sqlite"))):
    turn_mzml_duckdb(msfiles[0], "E:/mzsql/MTBLS10066/20220923_LEAP-POS_QC04.duckdb")
    turn_mzml_sqlite(msfiles[0], "E:/mzsql/MTBLS10066/20220923_LEAP-POS_QC04.sqlite")
    turn_mzml_parquet(msfiles[0], "E:/mzsql/MTBLS10066/20220923_LEAP-POS_QC04.parquet")

conn = sqlite3.connect("E:/mzsql/MTBLS10066/20220923_LEAP-POS_QC04.sqlite")
top_int_df = pd.read_sql_query("SELECT * FROM MS1 ORDER BY int DESC LIMIT 30000", conn)
top_masses = []
top_rts = []
for i in range(100):
    top_rts.append(top_int_df["rt"].iloc[0])
    top_mz = top_int_df["mz"].iloc[0]
    top_masses.append(top_mz)
    mzmin, mzmax = pmppm(top_mz, 10)
    top_int_df = top_int_df[((top_int_df["mz"]<mzmin) | (top_int_df["mz"]>mzmax))]

def get_multichrom_duckdb_loop(file, mzs, ppm):
    conn = duckdb.connect(file)
    all_chroms = []
    for mz_i in mzs:
        query = "SELECT mz, int, rt FROM MS1 WHERE mz BETWEEN ? AND ?"
        query_data = conn.execute(query, pmppm(mz_i, ppm)).fetchdf()
        query_data["mz_i"] = mz_i
        all_chroms.append(query_data)
    conn.close()
    return(pd.concat(all_chroms, ignore_index=True))

def get_multichrom_duckdb_join(file, mzs, ppm):
    mz_df = pd.DataFrame(mzs, columns=["mz_i"])
    mz_df["mzmin"], mz_df["mzmax"] = zip(*mz_df["mz_i"].apply(lambda mz_i: pmppm(mz_i, 5)))
    conn = duckdb.connect(file)
    query = """
        SELECT mz_df.mz_i, MS1.*
        FROM mz_df
        JOIN MS1 ON MS1.mz BETWEEN mz_df.mzmin AND mz_df.mzmax
    """
    result = conn.execute(query).df()
    conn.close()
    return(result)

def get_multichrom_duckdb_unified(file, mzs, ppm):
    conn = duckdb.connect(file)
    conditions = [
        f"(mz BETWEEN {pmppm(mz_i, ppm)[0]} AND {pmppm(mz_i, ppm)[1]})"
        for mz_i in mzs
    ]
    query = "SELECT * FROM MS1 WHERE " + " OR ".join(conditions)
    chrom_data = conn.execute(query).df()
    mz_ranges = [(mz, *pmppm(mz, 5)) for mz in mzs]
    chrom_data["mz_i"] = [next(mz for mz, mzmin, mzmax in mz_ranges if mzmin <= x <= mzmax) for x in chrom_data["mz"]]
    conn.close()
    return(chrom_data)

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

def get_multichrom_sqlite_join(file, mzs, ppm):
    mz_df = pd.DataFrame(mzs, columns=["mz_i"])
    mz_df["mzmin"], mz_df["mzmax"] = zip(*mz_df["mz_i"].apply(lambda mz_i: pmppm(mz_i, ppm)))
    conn = sqlite3.connect(file)
    mz_df.to_sql("mz_df", conn, if_exists="replace", index=False)
    query = """
        SELECT mz_df.mz_i, MS1.*
        FROM mz_df
        JOIN MS1 ON MS1.mz BETWEEN mz_df.mzmin AND mz_df.mzmax
    """
    result = pd.read_sql(query, conn)
    conn.cursor().execute("DROP TABLE mz_df")
    conn.close()
    return result

def get_multichrom_sqlite_unified(file, mzs, ppm):
    conn = sqlite3.connect(file)
    conditions = [
        f"(mz BETWEEN {pmppm(mz_i, ppm)[0]} AND {pmppm(mz_i, ppm)[1]})"
        for mz_i in mzs
    ]
    query = "SELECT * FROM MS1 WHERE " + " OR ".join(conditions)
    chrom_data = pd.read_sql(query, conn)
    conn.close()
    
    mz_ranges = [(mz, *pmppm(mz, 5)) for mz in mzs]
    chrom_data["mz_i"] = [next(mz for mz, mzmin, mzmax in mz_ranges if mzmin <= x <= mzmax) for x in chrom_data["mz"]]
    
    return chrom_data

def get_multichrom_parquet_loop(pqds_dir, mzs, ppm):
    scan_data =[]
    for mz_i in mzs:
        mz_min, mz_max = pmppm(mz_i, ppm)
        dataset = ds.dataset(f"{pqds_dir}/MS1", format="parquet")
        bet_df = dataset.to_table(filter=((ds.field("mslevel")=="MS1") & (ds.field('mz') >= mz_min) & (ds.field('mz') <= mz_max))).to_pandas()
        bet_df["mz_i"] = mz_i
        scan_data.append(bet_df)
    return pd.concat(scan_data, ignore_index=True)

def get_multichrom_parquet_unified(pqds_dir, mzs, ppm):
    filter_cond = None
    for mz_i in mzs:
        mz_min, mz_max = pmppm(mz_i, ppm)
        mz_condition = (ds.field('mz') >= mz_min) & (ds.field('mz') <= mz_max)
        # Combine conditions with OR
        if filter_cond is None:
            filter_cond = mz_condition
        else:
            filter_cond |= mz_condition
    dataset = ds.dataset(f"{pqds_dir}/MS1", format="parquet")
    chrom_data = dataset.to_table(filter=filter_cond).to_pandas()
    
    mz_ranges = [(mz, *pmppm(mz, 5)) for mz in mzs]
    chrom_data["mz_i"] = [next(mz for mz, mzmin, mzmax in mz_ranges if mzmin <= x <= mzmax) for x in chrom_data["mz"]]
    return chrom_data



init_df = pd.DataFrame(columns=["method", "time"])
init_df.to_csv('data/dbloopjoinuni_times.csv', index=False)

function_list = [
    ("get_multichrom_duckdb_loop", ".duckdb"),
    ("get_multichrom_duckdb_join", ".duckdb"),
    ("get_multichrom_duckdb_unified", ".duckdb"),
    ("get_multichrom_sqlite_loop", ".sqlite"),
    ("get_multichrom_sqlite_join", ".sqlite"),
    ("get_multichrom_sqlite_unified", ".sqlite"),
    ("get_multichrom_parquet_loop", "_pqds"),
    ("get_multichrom_parquet_unified", "_pqds"),
]


init_df = pd.DataFrame(columns=["method", "time"])
init_df.to_csv('data/dbloopjoinuni_times.csv', index=False)

random.seed(123)
for n_chroms in [1, 3, 10, 30, 100]:
    print(n_chroms)
    mz_list = random.sample(top_masses, n_chroms)
    for fun_i, file_ending in function_list:
        rep_function = f"{fun_i}('E:/mzsql/MTBLS10066/20220923_LEAP-POS_QC04{file_ending}', {mz_list}, 5)"
        time_vals = timeit.repeat(rep_function, globals=globals(), number=1, repeat=10)
        fun_name = fun_i.replace("get_multichrom_", "")
        time_df = pd.DataFrame({"n_chroms":n_chroms, "method":fun_name, "time":time_vals})
        time_df.to_csv('data/dbloopjoinuni_times.csv', mode='a', header=False, index=False)
