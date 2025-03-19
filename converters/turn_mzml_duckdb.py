import duckdb
import pandas as pd
import pyteomics.mzml
import os
def pmppm(mz, ppm):
    return(mz*(1-ppm/1000000), mz*(1+ppm/1000000))

def turn_mzml_duckdb(files, outfile, ordered=None):
    conn = duckdb.connect(outfile)
    conn.execute("DROP TABLE IF EXISTS MS1")
    conn.execute("DROP TABLE IF EXISTS MS2")
    conn.execute("DROP TABLE IF EXISTS file_info")
    conn.execute("DROP TABLE IF EXISTS scan_info")
    empty_ms1 = pd.DataFrame(columns=["filename", "id", "mz", "int", "rt"]).astype(
        {"filename": "string", "id": "int64", "mz": "float64", "int": "float64", "rt": "float64"}
    )
    conn.execute("CREATE TABLE MS1 AS SELECT * FROM empty_ms1")
    empty_ms2 = pd.DataFrame(columns=["filename", "id", "premz", "fragmz", "int", "rt"]).astype(
        {"filename": "string", "id": "int64", "premz": "float64", "fragmz": "float64", "int": "float64", "rt": "float64"}
    )
    conn.execute("CREATE TABLE MS2 AS SELECT * FROM empty_ms2")
    empty_file_info = pd.DataFrame(columns=['filename', 'n_scans', 'rt_start', 'rt_end']).astype(
        {'filename':'string', 'n_scans':'int64', 'rt_start':'float64', 'rt_end':'float64'}
    )
    conn.execute("CREATE TABLE file_info AS SELECT * FROM empty_file_info")
    empty_scan_info = pd.DataFrame(columns=[
      'filename', 'rt', 'scan_num', 'min_mz', 'max_mz', 'ms_level', 'centroided', 'polarity'
    ]).astype(
        {'filename':"string", 'rt':"float64",
         'scan_num':"int64", 'min_mz':"float64", 
         'max_mz':"float64",
         'ms_level':"int64",
         'centroided':"bool",
         'polarity':"string"}
    )
    conn.execute("CREATE TABLE scan_info AS SELECT * FROM empty_scan_info")

    
    
    if isinstance(files, str):
        files = [files]
    for file in files:
        open_file = pyteomics.mzml.MzML(file)
        file_info = []
        n_scans = len(open_file)
        min_rt = open_file[0]["scanList"]["scan"][0]["scan start time"]
        max_rt = open_file[n_scans-1]["scanList"]["scan"][0]["scan start time"]
        file_info.append(pd.DataFrame({'filename':[os.path.basename(file)], 'n_scans':n_scans,
                                       'rt_start':min_rt, 'rt_end':max_rt}))
        # Might need to switch to pymzml for more metadata...
    
        scan_info = []
        MS1_dfs = []
        MS2_dfs = []
        for spectrum in open_file:
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            idx = int(spectrum['id'].split("scan=")[-1].split()[0])
            scan_info.append(pd.DataFrame({'filename':[os.path.basename(file)], 'rt':rt_val,
                                           'scan_num':idx, 'min_mz':spectrum['lowest observed m/z'], 
                                           'max_mz':spectrum['highest observed m/z'],
                                           'ms_level':spectrum["ms level"],
                                           'centroided':'centroid spectrum' in spectrum,
                                           'polarity':("negative", "positive")['positive scan' in spectrum]}))
            
            if spectrum['ms level'] == 1:
                mz_vals = spectrum['m/z array']
                int_vals = spectrum['intensity array']
                df_scan = pd.DataFrame({'filename': os.path.basename(file), 'id': idx, 'mz': mz_vals, 
                                        'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
                MS1_dfs.append(df_scan)
            if spectrum['ms level'] == 2:
                mz_vals = spectrum['m/z array']
                int_vals = spectrum['intensity array']
                premz_val = spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
                df_scan = pd.DataFrame({'filename': os.path.basename(file), 'id': idx, 'premz': premz_val, 
                                        'fragmz': mz_vals, 'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
                MS2_dfs.append(df_scan)
        all_MS1 = pd.concat(MS1_dfs, ignore_index=True)
        all_MS2 = pd.concat(MS2_dfs, ignore_index=True)
        all_scan_info = pd.concat(scan_info, ignore_index=True)
        all_file_info = pd.concat(file_info, ignore_index=True)
        if ordered is not None:
            if ordered == "rt":
                all_MS1.sort_values("rt", inplace=True)
                all_MS2.sort_values("rt", inplace=True)
            if ordered == "mz":
                all_MS1.sort_values("mz", inplace=True)
            if ordered in ["fragmz", "premz"]:
                all_MS1.sort_values(ordered, inplace=True)
        conn.execute("INSERT INTO MS1 SELECT * FROM all_MS1")
        conn.execute("INSERT INTO MS2 SELECT * FROM all_MS2")
        conn.execute("INSERT INTO file_info SELECT * FROM all_file_info")
        conn.execute("INSERT INTO scan_info SELECT * FROM all_scan_info")
    conn.close()
    
    return outfile


# Example usage
import glob
ms_files = glob.glob("demo_data/*.mzML")
turn_mzml_duckdb(ms_files, "demo_data/msdata.duckdb")

# Example chromatogram extraction
conn = duckdb.connect("demo_data/msdata.duckdb")
chr_statement = "SELECT * FROM MS1 WHERE mz BETWEEN ? AND ?"
chrom = conn.execute(chr_statement, pmppm(118.0865, 10)).fetchdf()
conn.close()
import matplotlib.pyplot as plt
plt.plot(chrom["rt"], chrom["int"])
plt.show()
