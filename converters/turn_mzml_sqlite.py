import sqlite3
import pandas as pd
import pyteomics.mzml
import os
def pmppm(mz, ppm):
    return(mz*(1-ppm/1000000), mz*(1+ppm/1000000))

def turn_mzml_sqlite(files, outfile, ordered=None):
    conn = sqlite3.connect(outfile)
    conn.execute("DROP TABLE IF EXISTS MS1")
    conn.execute("DROP TABLE IF EXISTS MS2")
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
                                           'centroided':'centroid spectrum' in demo_lib,
                                           'polarity':("negative", "positive")['positive scan' in demo_lib]}))
            
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
        pd.concat(file_info, ignore_index=True).to_sql("file_info", conn, if_exists="append", index=False)
        pd.concat(scan_info, ignore_index=True).to_sql("scan_info", conn, if_exists="append", index=False)
        pd.concat(MS1_dfs, ignore_index=True).to_sql("MS1", conn, if_exists="append", index=False)
        pd.concat(MS2_dfs, ignore_index=True).to_sql("MS2", conn, if_exists="append", index=False)
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



import glob
ms_files = glob.glob("demo_data/*.mzML")
turn_mzml_sqlite(ms_files, "demo_data/msdata.sqlite")

# Example chromatogram extraction
conn = sqlite3.connect("demo_data/msdata.sqlite")
chr_statement = "SELECT * FROM MS1 WHERE mz BETWEEN ? AND ? ORDER BY filename, rt"
chrom = pd.read_sql(chr_statement, conn, params=pmppm(118.0865, 10))
conn.close()
import matplotlib.pyplot as plt
plt.plot(chrom["rt"], chrom["int"])
plt.show()
