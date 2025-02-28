
import numpy as np
import pandas as pd
import pyteomics.mzml
from .helpers import pmppm
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc
import pyarrow.dataset as ds
import os

def turn_mzml_parquet(files, outdir, ordered=None):
    if isinstance(files, str):
        files = [files]
    for file in files:
        MS1_dfs = []
        MS2_dfs = []
        for spectrum in pyteomics.mzml.MzML(file):
            idx = int(spectrum['id'].split("scan=")[-1].split()[0])
            mz_vals = spectrum['m/z array']
            int_vals = spectrum['intensity array']
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            if spectrum['ms level'] == 1:
                df_scan = pd.DataFrame({'id': idx, 'mslevel':"MS1", 'mz': mz_vals, 'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
                MS1_dfs.append(df_scan)
            if spectrum["ms level"] == 2:
                premz_val = spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
                df_scan = pd.DataFrame({'id': idx, 'mslevel':"MS2", 'premz': premz_val, 'fragmz': mz_vals, 'int': int_vals, 'rt': [rt_val] * len(mz_vals)})
                MS2_dfs.append(df_scan)
    
        all_MS1 = pd.concat(MS1_dfs, ignore_index=True)
        all_MS2 = pd.concat(MS2_dfs, ignore_index=True)
        if ordered is not None:
            if ordered == "rt":
                all_MS1.sort_values(by=ordered, inplace=True)
                all_MS2.sort_values(by=ordered, inplace=True)
            if ordered == "mz":
                all_MS1.sort_values(by=ordered, inplace=True)
            if ordered == "fragmz":
                all_MS2.sort_values(by=ordered, inplace=True)
            if ordered == "premz":
                all_MS2.sort_values(by=ordered, inplace=True)
    
        table_MS1 = pa.Table.from_pandas(all_MS1)
        table_MS2 = pa.Table.from_pandas(all_MS2)
    
        basename = os.path.splitext(os.path.basename(file))[0]
        os.makedirs(outdir, exist_ok=True)
        os.makedirs(f"{outdir}/MS1", exist_ok=True)
        os.makedirs(f"{outdir}/MS2", exist_ok=True)
        pq.write_table(table_MS1, f"{outdir}/MS1/{basename}.parquet")
        pq.write_table(table_MS2, f"{outdir}/MS2/{basename}.parquet")

    return outdir



def get_chrom_parquet(pqds_dir, mz, ppm):
    mz_min, mz_max = pmppm(mz, ppm)
    dataset = ds.dataset(f"{pqds_dir}/MS1", format="parquet")
    bet_df = dataset.to_table(filter=((ds.field("mslevel")=="MS1") & (ds.field('mz') >= mz_min) & (ds.field('mz') <= mz_max))).to_pandas()
    return bet_df

def get_spec_parquet(pqds_dir, spectrum_idx):
    ms1_ds = ds.dataset(f"{pqds_dir}/MS1", format="parquet")
    ms1_df = ms1_ds.to_table(filter=((ds.field('mslevel')=="MS1") & (ds.field('id') == spectrum_idx)), columns=["mz", "int"]).to_pandas()
    ms2_ds = ds.dataset(f"{pqds_dir}/MS2", format="parquet")
    ms2_df = ms2_ds.to_table(filter=((ds.field('mslevel')=="MS2") & (ds.field('id') == spectrum_idx)), columns=["fragmz", "int"]).to_pandas()
    ms2_df.columns = ["mz", "int"]
    if(ms1_df.shape[0]>0):
        spec_df = ms1_df
    else:
        spec_df = ms2_df
    return spec_df

def get_rtrange_parquet(pqds_dir, rtstart, rtend):
    dataset = ds.dataset(f"{pqds_dir}/MS1", format="parquet")
    rtrange_df = dataset.to_table(filter=((ds.field("mslevel")=="MS1") & (ds.field('rt') >= rtstart) & (ds.field('rt') <= rtend))).to_pandas()
    return rtrange_df



def get_MS2fragmz_parquet(pqds_dir, fragment_mz, ppm_acc):
    mz_min, mz_max = pmppm(fragment_mz, ppm_acc)
    dataset = ds.dataset(f"{pqds_dir}/MS2", format="parquet")
    bet_df = dataset.to_table(filter=((ds.field("mslevel")=="MS2") & (ds.field('fragmz') >= mz_min) & (ds.field('fragmz') <= mz_max))).to_pandas()
    return bet_df

def get_MS2premz_parquet(pqds_dir, precursor_mz, ppm_acc):
    mz_min, mz_max = pmppm(precursor_mz, ppm_acc)
    dataset = ds.dataset(f"{pqds_dir}/MS2", format="parquet")
    bet_df = dataset.to_table(filter=((ds.field("mslevel")=="MS2") & (ds.field('premz') >= mz_min) & (ds.field('premz') <= mz_max))).to_pandas()
    return bet_df

def get_MS2nloss_parquet(pqds_dir, nloss_mz, ppm_acc):
    mz_min, mz_max = pmppm(nloss_mz, ppm_acc)
    dataset = ds.dataset(f"{pqds_dir}/MS2", format="parquet")
    bet_df = dataset.to_table(filter=((ds.field("mslevel")=="MS2") & 
                                      (ds.field('premz')-ds.field('fragmz') >= mz_min) & 
                                      (ds.field('premz')-ds.field('fragmz') <= mz_max))).to_pandas()
    return bet_df
