import sqlite3
import pandas as pd
import pyteomics.mzml
import os
from .helpers import pmppm

def turn_mzml_sqlite(files, outfile, ordered=None):
    conn = sqlite3.connect(outfile)
    conn.execute("DROP TABLE IF EXISTS MS1")
    conn.execute("DROP TABLE IF EXISTS MS2")
    if isinstance(files, str):
        files = [files]
    for file in files:
        for spectrum in pyteomics.mzml.MzML(file):
            if spectrum['ms level'] == 1:
                idx = int(spectrum['id'].split("scan=")[-1].split()[0])
                mz_vals = spectrum['m/z array']
                int_vals = spectrum['intensity array']
                rt_val = spectrum['scanList']['scan'][0]['scan start time']
                df_scan = pd.DataFrame({'filename': os.path.basename(file), 'id': idx, 'mz': mz_vals, 
                                        'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
                df_scan.to_sql("MS1", conn, if_exists="append", index=False)
            if spectrum['ms level'] == 2:
                idx = int(spectrum['id'].split("scan=")[-1].split()[0])
                mz_vals = spectrum['m/z array']
                int_vals = spectrum['intensity array']
                rt_val = spectrum['scanList']['scan'][0]['scan start time']
                premz_val = spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
                df_scan = pd.DataFrame({'filename': os.path.basename(file), 'id': idx, 'premz': premz_val, 
                                        'fragmz': mz_vals, 'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
                df_scan.to_sql("MS2", conn, if_exists="append", index=False)
    
    if ordered is not None:
        index_name = f"idx_{ordered}"
        if ordered == "rt":
            conn.execute(f"CREATE INDEX IF NOT EXISTS {index_name} ON MS1 ({ordered})")
            conn.execute(f"CREATE INDEX IF NOT EXISTS {index_name} ON MS2 ({ordered})")
        if ordered == "mz":
            conn.execute(f"CREATE INDEX IF NOT EXISTS {index_name} ON MS1 ({ordered})")
        if ordered in ["fragmz", "premz"]:
            conn.execute(f"CREATE INDEX IF NOT EXISTS {index_name} ON MS2 ({ordered})")
    conn.close()

    return outfile


def get_chrom_sqlite(file, mz, ppm):
    """
    Extracts a chromatogram from an SQLite database based on an m/z value and tolerance.

    Parameters:
    - file (str): Path to the SQLite database file containing the MS1 data.
    - mz (float): The target m/z value for which to extract chromatographic data.
    - ppm (float): The parts-per-million (ppm) tolerance to calculate the m/z range.

    Returns:
    - pandas.DataFrame: A DataFrame containing rows from the MS1 table where the m/z value is within the specified range.

    Raises:
    - sqlite3.Error: If there is an issue with SQLite operations.
    - ValueError: If `mz` or `ppm` values are invalid.
    """
    if not os.path.exists(file):
        raise FileNotFoundError(f"Database '{file}' does not exist.")

    mz_min, mz_max = pmppm(mz, ppm)
    conn = sqlite3.connect(file)
    query = f"SELECT * FROM MS1 WHERE mz BETWEEN {mz_min} AND {mz_max}"
    query_data = pd.read_sql_query(query, conn)

    conn.close()
    return query_data

def get_spec_sqlite(file, spectrum_idx):
    """
    Retrieves a single spectrum from an SQLite database by spectrum ID.

    This function queries an SQLite database to extract all data from the `MS1` table
    corresponding to a specific spectrum ID.

    Parameters:
    - file (str): Path to the SQLite database file containing the MS1 data.
    - spectrum_idx (int): The ID of the spectrum to retrieve.

    Returns:
    - pandas.DataFrame: A DataFrame containing all rows from the MS1 table with the specified spectrum ID.

    Raises:
    - sqlite3.Error: If there is an issue with SQLite operations.
    - ValueError: If `spectrum_idx` is invalid.
    """
    if not os.path.exists(file):
        raise FileNotFoundError(f"Database '{file}' does not exist.")

    conn = sqlite3.connect(file)
    query = f"SELECT id, mz, int FROM MS1 WHERE id = {spectrum_idx} UNION ALL SELECT id, fragmz AS mz, int FROM MS2 WHERE id = {spectrum_idx}"
    spectrum_data = pd.read_sql_query(query, conn)
    conn.close()

    return spectrum_data

def get_rtrange_sqlite(file, rtstart, rtend):
    """
    Retrieves data from an SQLite database for a specific retention time (RT) range.

    This function queries an SQLite database to extract all rows from the `MS1` table
    where the retention time (RT) is within the specified range.

    Parameters:
    - file (str): Path to the SQLite database file containing the MS1 data.
    - rtstart (float): The starting retention time for the range (inclusive).
    - rtend (float): The ending retention time for the range (inclusive).

    Returns:
    - pandas.DataFrame: A DataFrame containing all rows from the MS1 table 
      where the retention time falls within the specified range.

    Raises:
    - sqlite3.Error: If there is an issue with SQLite operations.
    - ValueError: If `rtstart` or `rtend` are invalid.
    """
    if not os.path.exists(file):
        raise FileNotFoundError(f"Database '{file}' does not exist.")

    conn = sqlite3.connect(file)
    query = f"SELECT * FROM MS1 WHERE rt >= {rtstart} AND rt <= {rtend}"
    rt_range_data = pd.read_sql_query(query, conn)
    
    conn.close()
    
    return rt_range_data


def get_MS2fragmz_sqlite(file, fragment_mz, ppm_acc):
    conn = sqlite3.connect(file)
    mzmin, mzmax = pmppm(fragment_mz, ppm_acc)
    query = f"SELECT * FROM MS2 WHERE fragmz BETWEEN {mzmin} AND {mzmax}"
    query_data = pd.read_sql_query(query, conn)
    conn.close()
    return(query_data)

def get_MS2premz_sqlite(file, precursor_mz, ppm_acc):
    conn = sqlite3.connect(file)
    mzmin, mzmax = pmppm(precursor_mz, ppm_acc)
    query = f"SELECT * FROM MS2 WHERE premz BETWEEN {mzmin} AND {mzmax}"
    query_data = pd.read_sql_query(query, conn)
    conn.close()
    return(query_data)

def get_MS2nloss_sqlite(file, nloss_mz, ppm_acc):
    conn = sqlite3.connect(file)
    mzmin, mzmax = pmppm(nloss_mz, ppm_acc)
    query = f"SELECT * FROM MS2 WHERE premz-fragmz BETWEEN {mzmin} AND {mzmax}"
    query_data = pd.read_sql_query(query, conn)
    conn.close()
    return(query_data)