
import numpy as np
import pandas as pd
import sqlite3
import pyteomics.mzml
from .helpers import pmppm

def turn_mzml_sqlite(file, outfile, ordered = False ):
    #converts mzml file to sqlite database. Here for reference, and to show change for row id
    conn = sqlite3.connect(outfile)
    conn.execute("DROP TABLE IF EXISTS MS1")
    for spectrum in pyteomics.mzml.MzML(file):
        if spectrum['ms level'] == 1:
            #print(spectrum["index"])
            idx = int(spectrum['id'].split("scan=")[-1].split()[0])
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            df_scan = pd.DataFrame({'id':idx,'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(mz_vals)})
            df_scan.to_sql("MS1", conn, if_exists="append", index=False)
    if ordered:
        conn.execute("CREATE INDEX IF NOT EXISTS idx_mz ON MS1 (mz)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_int ON MS1 (int)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_rt ON MS1 (rt)")
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

def get_spec_sqlite(file, spectrum_idx):
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