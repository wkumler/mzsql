
import numpy as np
import pandas as pd
import pyteomics.mzml
from .helpers import pmppm
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc
import pyarrow.dataset as ds
import os

def turn_mzml_parquet(file, outdir, ordered=None):
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

    table_MS1 = pa.Table.from_pandas(all_MS1)
    table_MS2 = pa.Table.from_pandas(all_MS2)

    basename = os.path.splitext(os.path.basename(file))[0]
    os.makedirs(outdir)
    os.makedirs(f"{outdir}/mslevel=MS1", exist_ok=True)
    os.makedirs(f"{outdir}/mslevel=MS2", exist_ok=True)
    pq.write_table(table_MS1, f"{outdir}/mslevel=MS1/{basename}.parquet")
    pq.write_table(table_MS2, f"{outdir}/mslevel=MS2/{basename}.parquet")

    return outdir



def get_chrom_parquet(pqds_dir, mz, ppm):
    mz_min, mz_max = pmppm(mz, ppm)
    dataset = ds.dataset(pqds_dir, format="parquet", partitioning="hive")
    bet_df = dataset.to_table(filter=((ds.field("mslevel")=="MS1") & (ds.field('mz') >= mz_min) & (ds.field('mz') <= mz_max))).to_pandas()
    return bet_df

def get_spec_parquet(pqds_dir, spectrum_idx):
    dataset = ds.dataset(pqds_dir, format="parquet", partitioning="hive")
    spec_df = dataset.to_table(filter=((ds.field("mslevel")=="MS1") & (ds.field('id') == spectrum_idx))).to_pandas()
    return spec_df

def get_rtrange_parquet(pqds_dir, rtstart, rtend):
    dataset = ds.dataset(pqds_dir, format="parquet", partitioning="hive")
    rtrange_df = dataset.to_table(filter=((ds.field("mslevel")=="MS1") & (ds.field('rt') >= rtstart) & (ds.field('rt') <= rtend))).to_pandas()
    return rtrange_df


def get_MS2scan_parquet(pqds_dir, spectrum_idx):
    dataset = ds.dataset(pqds_dir, format="parquet", partitioning="hive")
    spec_df = dataset.to_table(filter=((ds.field("mslevel")=="MS1") & (ds.field('id') == spectrum_idx))).to_pandas()
    return spec_df

def get_MS2fragmz_parquet(pqds_dir, fragment_mz, ppm_acc):
    mz_min, mz_max = pmppm(fragment_mz, ppm_acc)
    dataset = ds.dataset(pqds_dir, format="parquet", partitioning="hive")
    bet_df = dataset.to_table(filter=((ds.field("mslevel")=="MS2") & (ds.field('fragmz') >= mz_min) & (ds.field('fragmz') <= mz_max))).to_pandas()
    return bet_df

def get_MS2premz_parquet(pqds_dir, precursor_mz, ppm_acc):
    mz_min, mz_max = pmppm(precursor_mz, ppm_acc)
    dataset = ds.dataset(pqds_dir, format="parquet", partitioning="hive")
    bet_df = dataset.to_table(filter=((ds.field("mslevel")=="MS2") & (ds.field('premz') >= mz_min) & (ds.field('premz') <= mz_max))).to_pandas()
    return bet_df

def get_MS2nloss_parquet(pqds_dir, nloss_mz, ppm_acc):
    mz_min, mz_max = pmppm(nloss_mz, ppm_acc)
    dataset = ds.dataset(pqds_dir, format="parquet", partitioning="hive")
    bet_df = dataset.to_table(filter=((ds.field("mslevel")=="MS2") & 
                                      (ds.field('premz')-ds.field('fragmz') >= mz_min) & 
                                      (ds.field('premz')-ds.field('fragmz') <= mz_max))).to_pandas()
    return bet_df
