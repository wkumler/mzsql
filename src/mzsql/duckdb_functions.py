
import numpy as np
import pandas as pd
import duckdb
import pyteomics.mzml
from .helpers import pmppm

def turn_mzml_duckdb(file, outfile, ordered=None):
    """
    Converts an mzML file into a DuckDB database table.
    
    Args:
        file (str): Path to the mzML file to be processed.
        outfile (str): Path to the output DuckDB file where the data will be stored.
        ordered (bool, optional): If True, the mz values are sorted before storing. Defaults to False.
    
    Returns:
        str: The path to the output DuckDB file.
    """
    conn = duckdb.connect(outfile)
    conn.execute("DROP TABLE IF EXISTS MS1")

    MS1_dfs = []
    MS2_dfs = []
    for spectrum in pyteomics.mzml.MzML(file):
        if spectrum['ms level'] == 1:
            idx = int(spectrum['id'].split("scan=")[-1].split()[0])
            mz_vals = spectrum['m/z array']
            int_vals = spectrum['intensity array']
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            df_scan = pd.DataFrame({'id': idx, 'mz': mz_vals, 'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
            MS1_dfs.append(df_scan)
        if spectrum["ms level"] == 2:
            idx = int(spectrum['id'].split("scan=")[-1].split()[0])
            premz_val = spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
            mz_vals = spectrum['m/z array']
            int_vals = spectrum['intensity array']
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            df_scan = pd.DataFrame({'id': idx, 'premz': premz_val, 'fragmz': mz_vals, 'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
            MS2_dfs.append(df_scan)

    all_MS1 = pd.concat(MS1_dfs, ignore_index=True)
    all_MS2 = pd.concat(MS2_dfs, ignore_index=True)
    if ordered is not None:
        all_MS1.sort_values(by=ordered, inplace=True)
        all_MS2.sort_values(by=ordered, inplace=True)

    conn.execute("CREATE TABLE MS1 AS SELECT * FROM all_MS1")
    conn.execute("CREATE TABLE MS2 AS SELECT * FROM all_MS2")

    conn.close()

    return outfile

def get_chrom_duckdb(file, mz, ppm):
    """
    Retrieves chromatogram data from the DuckDB database for a specified m/z range.

    Args:
        file (str): Path to the DuckDB file.
        mz (float): The target m/z value.
        ppm (float): The allowed mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing m/z, intensity, and retention time data.
    """
    conn = duckdb.connect(file)
    mz_min, mz_max = pmppm(mz, ppm)
    query = "SELECT mz, int, rt FROM MS1 WHERE mz BETWEEN ? AND ?"
    query_data = conn.execute(query, (mz_min, mz_max)).fetchdf()    
    conn.close()
    return query_data

def get_spec_duckdb(file, spectrum_idx):
    """
    Retrieves the spectrum data for a specific scan ID from the DuckDB database.

    Args:
        file (str): Path to the DuckDB file.
        spectrum_idx (int): The scan ID of the spectrum to retrieve.

    Returns:
        pd.DataFrame: A DataFrame containing all columns for the specified spectrum.
    """
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS1 WHERE id = ?"
    spectrum_data = conn.execute(query, (spectrum_idx,)).fetchdf()
    conn.close()
    return spectrum_data

def get_rtrange_duckdb(file, rtstart, rtend):
    """
    Retrieves spectrum data within a specified retention time range from the DuckDB database.

    Args:
        file (str): Path to the DuckDB file.
        rtstart (float): The start retention time value.
        rtend (float): The end retention time value.

    Returns:
        pd.DataFrame: A DataFrame containing all columns for spectra within the specified retention time range.
    """
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS1 WHERE rt BETWEEN ? AND ?"
    rt_range_data = conn.execute(query, (rtstart, rtend)).fetchdf()
    conn.close()
    return rt_range_data



def get_MS2scan_duckdb(file, spectrum_idx):
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS2 WHERE id = ?"
    spectrum_data = conn.execute(query, (spectrum_idx,)).fetchdf()
    conn.close()
    return(spectrum_data)

def get_MS2fragmz_duckdb(file, fragment_mz, ppm_acc):
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS2 WHERE fragmz BETWEEN ? AND ?"
    spectrum_data = conn.execute(query, pmppm(fragment_mz, ppm_acc)).fetchdf()
    conn.close()
    return(spectrum_data)

def get_MS2premz_duckdb(file, precursor_mz, ppm_acc):
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS2 WHERE premz BETWEEN ? AND ?"
    spectrum_data = conn.execute(query, pmppm(precursor_mz, ppm_acc)).fetchdf()
    conn.close()
    return(spectrum_data)

def get_MS2nloss_duckdb(file, nloss_mz, ppm_acc):
    conn = duckdb.connect(file)
    query = "SELECT * FROM MS2 WHERE (premz-fragmz) BETWEEN ? AND ?"
    spectrum_data = conn.execute(query, pmppm(nloss_mz, ppm_acc)).fetchdf()
    conn.close()
    return(spectrum_data)
