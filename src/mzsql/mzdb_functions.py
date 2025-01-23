
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
    AND spectrum.ms_level = 1
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
    AND spectrum.ms_level = 1
    AND spectrum.time BETWEEN ? AND ?
    """
    bb_dataframe = pd.read_sql(bb_query, connection, params=(rtstart*60, rtend*60))
    unpacked_bb_list = [unpack_raw_bb(bb_data) for bb_data in bb_dataframe["data"]]
    rtrange_data = pd.concat(unpacked_bb_list).merge(scanid_rt_pd)
    rtrange_data["rt"] /= 60

    connection.close()
    return(rtrange_data)


def parse_mzDB_premz_string(precursor_xml):
    start_index = precursor_xml.find('name="selected ion m/z"')
    value_index = precursor_xml.find('value="', start_index) + len('value="')
    end_index = precursor_xml.find('"', value_index)
    premz_val = float(precursor_xml[value_index:end_index])
    return(premz_val)

def get_MS2scan_mzdb(mzdb_file, spectrum_idx):
    connection = sqlite3.connect(mzdb_file)
    spec_bb_query = "SELECT bb_first_spectrum_id, precursor_list FROM spectrum WHERE id = ?"
    spectrum_metadata = pd.read_sql(spec_bb_query, connection, params=(spectrum_idx+1,))
    bb_id_for_scan = spectrum_metadata["bb_first_spectrum_id"][0]
    precursor_xml = spectrum_metadata["precursor_list"][0]

    premz_val = parse_mzDB_premz_string(precursor_xml)
    
    bb_query = """
    SELECT first_spectrum_id, data
    FROM bounding_box
    WHERE first_spectrum_id = ?
    """
    
    bb_dataframe = pd.read_sql(bb_query, connection, params=(str(bb_id_for_scan),))
    unpacked_bb_list = [unpack_raw_bb(bb_data) for bb_data in bb_dataframe["data"]]
    bb_spec = pd.concat(unpacked_bb_list)
    bb_spec.columns = ["fragmz", "int", "scan_id"]
    scan_data = bb_spec[bb_spec["scan_id"]==spectrum_idx+1]
    scan_data["premz"] = premz_val
    
    connection.close()
    return(scan_data)

def get_MS2premz_mzdb(mzdb_file, precursor_mz, ppm_acc):
    connection = sqlite3.connect(mzdb_file)
    mzmin, mzmax = pmppm(precursor_mz, ppm_acc)
    spec_bb_query = "SELECT bb_first_spectrum_id, precursor_list FROM spectrum WHERE ms_level = 2"
    spectrum_metadata = pd.read_sql(spec_bb_query, connection)
    premz_vals = np.array([parse_mzDB_premz_string(x) for x in spectrum_metadata["precursor_list"]])
    chosen_scans = spectrum_metadata[(mzmin < premz_vals) & (premz_vals < mzmax)].astype("str")
    
    bb_query = """
    SELECT first_spectrum_id, data
    FROM bounding_box
    WHERE first_spectrum_id IN ({})
    """.format(','.join(chosen_scans["bb_first_spectrum_id"]))
    
    bb_dataframe = pd.read_sql(bb_query, connection)
    unpacked_bb_list = [unpack_raw_bb(bb_data) for bb_data in bb_dataframe["data"]]
    bb_spec = pd.concat(unpacked_bb_list)
    bb_spec.columns = ["fragmz", "int", "scan_id"]
    
    connection.close()
    return(bb_spec)

def get_MS2fragmz_mzdb(mzdb_file, fragment_mz, ppm_acc):
    connection = sqlite3.connect(mzdb_file)
    spec_bb_query = "SELECT bb_first_spectrum_id, precursor_list FROM spectrum WHERE ms_level = 2"
    spectrum_metadata = pd.read_sql(spec_bb_query, connection)
    ms2_query = '''
    SELECT data, precursor_list
    FROM bounding_box, spectrum
    WHERE ms_level = 2
    AND bounding_box.first_spectrum_id = spectrum.id
    '''
    ms2_bbs = pd.read_sql(ms2_query, connection)
    mzmin, mzmax = pmppm(fragment_mz, ppm_acc)
    scan_dfs = []
    premz_vals = np.array([parse_mzDB_premz_string(x) for x in ms2_bbs["precursor_list"]])
    for index, bb_data in enumerate(ms2_bbs["data"]):
        scan_data = unpack_raw_bb(bb_data)
        bet_idxs = (mzmin < scan_data["mz"]) & (scan_data["mz"] < mzmax)
        frag_data = scan_data[bet_idxs].copy()
        frag_data.columns = ["fragmz", "int", "scan_id"]
        frag_data["premz"] = premz_vals[index]
        scan_dfs.append(frag_data)
    spectrum_data = pd.concat(scan_dfs)
    return(spectrum_data)

def get_MS2nloss_mzdb(mzdb_file, neutral_loss, ppm_acc):
    connection = sqlite3.connect(mzdb_file)
    spec_bb_query = "SELECT bb_first_spectrum_id, precursor_list FROM spectrum WHERE ms_level = 2"
    spectrum_metadata = pd.read_sql(spec_bb_query, connection)
    ms2_query = '''
    SELECT data, precursor_list
    FROM bounding_box, spectrum
    WHERE ms_level = 2
    AND bounding_box.first_spectrum_id = spectrum.id
    '''
    ms2_bbs = pd.read_sql(ms2_query, connection)
    mzmin, mzmax = pmppm(neutral_loss, ppm_acc)
    scan_dfs = []
    premz_vals = np.array([parse_mzDB_premz_string(x) for x in ms2_bbs["precursor_list"]])
    for index, bb_data in enumerate(ms2_bbs["data"]):
        scan_data = unpack_raw_bb(bb_data)
        bet_idxs = (mzmin < premz_vals[index]-scan_data["mz"]) & (premz_vals[index]-scan_data["mz"] < mzmax)
        frag_data = scan_data[bet_idxs].copy()
        frag_data.columns = ["fragmz", "int", "scan_id"]
        frag_data["premz"] = premz_vals[index]
        scan_dfs.append(frag_data)
    spectrum_data = pd.concat(scan_dfs)
    return(spectrum_data)

