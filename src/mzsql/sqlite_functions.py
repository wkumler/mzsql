import sqlite3
import pandas as pd
import pyteomics.mzml
import os
from .helpers import pmppm

def turn_mzml_sqlite(file, outfile, ordered=None):
    """
    Converts an mzML file into an SQLite database.

    This function reads an mzML file, extracts MS1 spectra data, and stores it in an SQLite database.
    Optionally, it can create an index on a specified column ('mz', 'int', or 'rt') for faster querying.

    Parameters:
    - file (str): Path to the input mzML file.
    - outfile (str): Path to the output SQLite database file.
    - ordered (str, optional): Column name to create an index on. Must be one of 'mz', 'int', or 'rt'.
                               If None, no index is created. Default is None.

    Returns:
    - str: The path to the created SQLite database file.

    Raises:
    - sqlite3.Error: If there is an issue with SQLite operations.
    - KeyError: If the input mzML file does not contain expected fields.
    - ValueError: If the ordered column name is invalid.
    """
    if not os.path.exists(file):
        raise FileNotFoundError(f"Input mzML file '{file}' does not exist.")

    conn = sqlite3.connect(outfile)
    try:
        conn.execute("DROP TABLE IF EXISTS MS1")

        for spectrum in pyteomics.mzml.MzML(file):
            if spectrum['ms level'] == 1:
                try:
                    idx = int(spectrum['id'].split("scan=")[-1].split()[0])
                    mz_vals = spectrum['m/z array']
                    int_vals = spectrum['intensity array']
                    rt_val = spectrum['scanList']['scan'][0]['scan start time']
                except KeyError as e:
                    raise KeyError(f"Missing expected key in mzML spectrum data: {e}")

                df_scan = pd.DataFrame({'id': idx, 'mz': mz_vals, 'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
                df_scan.to_sql("MS1", conn, if_exists="append", index=False)

        if ordered in ['mz', 'int', 'rt']:
            index_name = f"idx_{ordered}"
            conn.execute(f"CREATE INDEX IF NOT EXISTS {index_name} ON MS1 ({ordered})")
        elif ordered is not None:
            raise ValueError(f"Invalid column for indexing: {ordered}. Must be one of 'mz', 'int', or 'rt'.")
    except sqlite3.Error as e:
        raise sqlite3.Error(f"SQLite error occurred: {e}")
    finally:
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
    query = f"SELECT * FROM MS1 WHERE id = {spectrum_idx}"
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
