
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyteomics import mzml

def pmppm(mass, ppm=4):
    return(mass*(1-ppm/1000000), mass*(1+ppm/1000000))


# mzML things
def get_chrom_mzml(file, mz, ppm):
    scan_dfs = []
    for spectrum in mzml.MzML(file):
        if spectrum['ms level'] == 1:
            rt_val = spectrum['scanList']['scan'][0]['scan start time']
            mz_vals=spectrum['m/z array']
            int_vals = spectrum['intensity array']

            mzmin, mzmax = pmppm(mz, ppm)
            bet_idxs = (spectrum["m/z array"]>mzmin) & (spectrum["m/z array"]<mzmax)
            if(sum(bet_idxs)>0):
                df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
                scan_dfs.append(df_scan)    
    return(pd.concat(scan_dfs, ignore_index=True))

def get_spectrum_mzml(file, spectrum_idx):
    raise Exception("mzML spectrum extraction yet implemented")

def get_rtrange_mzml(file, rtstart, rtend):
    raise Exception("mzML rtrange extraction yet implemented")

