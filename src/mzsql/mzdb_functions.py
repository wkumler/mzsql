
import sqlite3
import pandas as pd
import numpy as np
import struct
from .helpers import pmppm


def unpack_raw_bb(raw_string):
    """
    Unpacks raw bounding box (BB) data from a binary string and converts it into a DataFrame.

    Args:
        raw_string (bytes): A raw binary string containing peak data for multiple scans.

    Returns:
        pd.DataFrame: A DataFrame with columns ['mz', 'int', 'scan_id'] representing 
                      m/z, intensity, and scan ID for each peak.
    """
    offset = 0
    scanids = []
    offsets = []
    nb_peaks_all = []
    while offset < len(raw_string):
        # Unpack two integers starting at a byte offset of zero
        scan_id, nb_peaks = struct.Struct('<2i').unpack_from(raw_string, offset)
        scanids.append(scan_id)
        nb_peaks_all.append(nb_peaks)
        offsets.append(offset)
        # New offset is 8 bytes (metadata) + 12 * number of peaks (8 for mz and 4 for int)
        offset+=8+12*nb_peaks
    
    scan_data = []
    for i in range(len(scanids)):
        # Encoded as little endian, then double (mz) and float (int) for each peak
        decomp_string = '<'+'df'*nb_peaks_all[i]
        # Unpack the data, remembering to skip the 8 bytes of metadata at the front!
        data_tuple = struct.Struct(decomp_string).unpack_from(raw_string, offsets[i]+8)
        mz_int_array = np.array(data_tuple).reshape(-1, 2)
        mz_int_scanid_array = np.column_stack((mz_int_array, np.full(mz_int_array.shape[0], scanids[i])))
        # Getting duplicate entries for some reason... using unique() to remove them for now
        scan_data.append(np.unique(mz_int_scanid_array, axis=0))
    
    return(pd.DataFrame(np.vstack(scan_data), columns=['mz', 'int', 'scan_id']))


def get_chrom_mzdb(file, mz, ppm):
    """
    Retrieves chromatogram data from an mzDB file for a specified m/z range.

    Args:
        file (str): Path to the mzDB file.
        mz (float): The target m/z value.
        ppm (float): The allowed mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing m/z, intensity, and retention time data 
                      for the specified m/z range.
    """
    connection = sqlite3.connect(file)
    cursor = connection.cursor()

    mzmin, mzmax = pmppm(mz, ppm)
    spec_bb_query = "SELECT MAX(begin_mz) FROM run_slice WHERE begin_mz < ?"
    # Assumes that the ppm range doesn't span a bounding box...
    bb_id_for_chrom = cursor.execute(spec_bb_query, (mzmin,)).fetchone()[0]

    # Get mapping between scan_id and rt
    scanid_rt_pd = pd.read_sql("SELECT id AS scan_id, time AS rt, ms_level FROM spectrum", connection)

    bb_query = """
    SELECT bounding_box.id, begin_mz, end_mz, first_spectrum_id, time AS first_spectrum_time, data
    FROM run_slice, bounding_box, spectrum
    WHERE bounding_box.run_slice_id = run_slice.id
    AND bounding_box.first_spectrum_id = spectrum.id
    AND begin_mz = ?
    """
    bb_dataframe = pd.read_sql(bb_query, connection, params=(bb_id_for_chrom,))
    unpacked_bb_list = [unpack_raw_bb(bb_data) for bb_data in bb_dataframe["data"]]
    bb_chrom = pd.concat(unpacked_bb_list).merge(scanid_rt_pd)
    bb_chrom["rt"] /= 60

    connection.close()
    
    return(bb_chrom[(mzmin < bb_chrom["mz"]) & (bb_chrom["mz"] < mzmax)])

def get_spec_mzdb(file, scan_num):
    """
    Retrieves spectrum data for a specific scan number from an mzDB file.

    Args:
        file (str): Path to the mzDB file.
        scan_num (int): The scan number to retrieve the spectrum for.

    Returns:
        pd.DataFrame: A DataFrame containing m/z and intensity values for the specified scan.
    """
    connection = sqlite3.connect(file)
    cursor = connection.cursor()
    spec_bb_query = "SELECT bb_first_spectrum_id FROM spectrum WHERE id = ?"
    bb_id_for_scan = cursor.execute(spec_bb_query, (scan_num,)).fetchone()[0]
    
    bb_query = """
    SELECT first_spectrum_id, data
    FROM bounding_box
    WHERE first_spectrum_id = ?
    """
    
    bb_dataframe = pd.read_sql(bb_query, connection, params=(bb_id_for_scan,))
    unpacked_bb_list = [unpack_raw_bb(bb_data) for bb_data in bb_dataframe["data"]]
    bb_spec = pd.concat(unpacked_bb_list)

    connection.close()
    return(bb_spec[bb_spec["scan_id"]==scan_num])

def get_rtrange_mzdb(file, rtstart, rtend):
    """
    Retrieves spectrum data within a specified retention time range from an mzDB file.

    Args:
        file (str): Path to the mzDB file.
        rtstart (float): The start retention time value in minutes.
        rtend (float): The end retention time value in minutes.

    Returns:
        pd.DataFrame: A DataFrame containing m/z, intensity, and retention time values for 
                      spectra within the specified retention time range.
    """
    connection = sqlite3.connect(file)
    scanid_rt_pd = pd.read_sql("SELECT id AS scan_id, time AS rt, ms_level FROM spectrum", connection)

    bb_query = """
    SELECT DISTINCT bounding_box.data
    FROM bounding_box, spectrum
    WHERE bounding_box.first_spectrum_id = spectrum.bb_first_spectrum_id
    AND spectrum.time BETWEEN ? AND ?
    """
    bb_dataframe = pd.read_sql(bb_query, connection, params=(rtstart*60, rtend*60))
    unpacked_bb_list = [unpack_raw_bb(bb_data) for bb_data in bb_dataframe["data"]]
    rtrange_data = pd.concat(unpacked_bb_list).merge(scanid_rt_pd)
    rtrange_data["rt"] /= 60

    connection.close()
    return(rtrange_data)