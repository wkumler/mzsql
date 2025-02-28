
import pandas as pd
import numpy as np
import timeit
import glob
import random
import pyopenms
import pyteomics.mzml
import pyteomics.mzmlb
import pymzml
import h5py
import matplotlib.pyplot as plt
from mzsql import *

random.seed(123)
basename=random.sample(glob.glob("E:/mzsql/MTBLS10066/*.mzML"), 1)[0].replace(".mzML", "").replace("\\", "/")
print(basename)

def get_randscan_mzml_pyopenms(basename, scan_num):
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(f"{basename}.mzML", exp)
    spec_data = exp[scan_num].get_peaks()
    return(pd.DataFrame({"mz":spec_data[0], "int":spec_data[1]}))

def get_randscan_mzml_pyteomics(basename, scan_num):
    file_data = pyteomics.mzml.MzML(f"{basename}.mzML")
    return(pd.DataFrame({"mz":file_data[scan_num]['m/z array'], "int":file_data[scan_num]['intensity array']}))

def get_randscan_mzml_pymzml(basename, scan_num):
    run = pymzml.run.Reader(f"{basename}.mzML", build_index_from_scratch=True)
    spec1_data = run[scan_num].peaks("raw")
    return(pd.DataFrame({"mz":spec1_data[:,0], "int":spec1_data[:,1]}))

def get_randscan_mzmlb(basename, scan_num):
    file_data = pyteomics.mzmlb.MzMLb(f"{basename}.mzMLb")
    return(pd.DataFrame({"mz":file_data[scan_num]['m/z array'], "int":file_data[scan_num]['intensity array']}))

def get_randscan_mza(basename, scan_num):
    mza = h5py.File(f"{basename}.mza", 'r')
    intensities = mza[f"Arrays_intensity/{scan_num}"][:]
    mz = mza[f"Arrays_mz/{scan_num}"][:]
    mza.close()
    spec_df = pd.DataFrame({"mz":mz, "int":intensities})
    return(spec_df)

def get_randscan_mz5(basename, scan_num):
    mz5_file = h5py.File(f"{basename}.mz5", 'r')
    scan_idxs = np.concatenate(([0], mz5_file["SpectrumIndex"][...]))
    lower_bound = scan_idxs[scan_num-1]
    upper_bound = scan_idxs[scan_num]
    mz_vals = np.cumsum(mz5_file["SpectrumMZ"][lower_bound:upper_bound])
    int_vals = mz5_file["SpectrumIntensity"][lower_bound:upper_bound]
    return(pd.DataFrame({"mz":mz_vals, "int":int_vals}))



def min_randscan_mzml_pyopenms(scan_num):
    spec_data = pyop_exp[scan_num].get_peaks()
    return(pd.DataFrame({"mz":spec_data[0], "int":spec_data[1]}))

def min_randscan_mzml_pyteomics(scan_num):
    return(pd.DataFrame({"mz":pyteo_data[scan_num]['m/z array'], "int":pyteo_data[scan_num]['intensity array']}))

def min_randscan_mzml_pymzml(scan_num):
    spec1_data = pymzml_run[scan_num].peaks("raw")
    return(pd.DataFrame({"mz":spec1_data[:,0], "int":spec1_data[:,1]}))

def min_randscan_mzmlb(scan_num):
    return(pd.DataFrame({"mz":mzmlb_data[scan_num]['m/z array'], "int":mzmlb_data[scan_num]['intensity array']}))

def min_randscan_mza(scan_num):
    intensities = mza_file["Arrays_intensity/"+str(scan_num)][:]
    mz = mza_file["Arrays_mz/"+str(scan_num)][:]
    spec_df = pd.DataFrame({"mz":mz, "int":intensities})
    return(spec_df)

def min_randscan_mz5(scan_num):
    scan_idxs = np.concatenate(([0], mz5_file["SpectrumIndex"][...]))
    lower_bound = scan_idxs[scan_num-1]
    upper_bound = scan_idxs[scan_num]
    mz_vals = np.cumsum(mz5_file["SpectrumMZ"][lower_bound:upper_bound])
    int_vals = mz5_file["SpectrumIntensity"][lower_bound:upper_bound]
    return(pd.DataFrame({"mz":mz_vals, "int":int_vals}))

# Manually specify the same random scans here from the singlefile_timing.py script
rand_MS1_scans = (9044, 187, 7532, 2396, 10409, 87, 5259, 69, 5157, 7988)
rand_MS2_scans = (3778, 3815, 5685, 6859, 4223, 3372, 6948, 6354, 3140, 912)

all_timings = []
for scan_num in rand_MS1_scans:
    print(scan_num)
    pyteo=timeit.repeat(f'get_randscan_mzml_pyteomics("{basename}", {scan_num})', globals=globals(), number=1, repeat=3)
    pyopen=timeit.repeat(f'get_randscan_mzml_pyopenms("{basename}", {scan_num})', globals=globals(), number=1, repeat=3)
    pymzml_data=timeit.repeat(f'get_randscan_mzml_pymzml("{basename}", {scan_num})', globals=globals(), number=1, repeat=3)
    mzmlb=timeit.repeat(f'get_randscan_mzmlb("{basename}", {scan_num})', globals=globals(), number=1, repeat=3)
    mza=timeit.repeat(f'get_randscan_mza("{basename}", {scan_num})', globals=globals(), number=1, repeat=3)
    mz5=timeit.repeat(f'get_randscan_mz5("{basename}", {scan_num})', globals=globals(), number=1, repeat=3)
    time_info = pd.DataFrame({"id":scan_num, "pyteomics": pyteo, "pyopenms": pyopen, "pymzml": pymzml_data, 
                              "mzMLb":mzmlb, "MZA": mza, "mz5": mz5, "fun_type": "inclusive"})
    all_timings.append(time_info)

for scan_num in rand_MS1_scans:
    print(scan_num)
    pyop_exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(f"{basename}.mzML", pyop_exp)
    pyopen=timeit.repeat(f'min_randscan_mzml_pyopenms({scan_num})', globals=globals(), number=1, repeat=3)

    pyteo_data = pyteomics.mzml.MzML(f"{basename}.mzML")
    pyteo=timeit.repeat(f"min_randscan_mzml_pyteomics({scan_num})", globals=globals(), number=1, repeat=3)

    pymzml_run = pymzml.run.Reader(f"{basename}.mzML", build_index_from_scratch=True)
    pymzml_data=timeit.repeat(f'min_randscan_mzml_pymzml({scan_num})', globals=globals(), number=1, repeat=3)

    mzmlb_data = pyteomics.mzmlb.MzMLb(f"{basename}.mzMLb")
    mzmlb=timeit.repeat(f"min_randscan_mzmlb({scan_num})", globals=globals(), number=1, repeat=3)
    
    mza_file = h5py.File(f"{basename}.mza", 'r')
    mza=timeit.repeat(f'min_randscan_mza({scan_num})', globals=globals(), number=1, repeat=3)

    mz5_file = h5py.File(f"{basename}.mz5", 'r')
    mz5=timeit.repeat(f'min_randscan_mz5({scan_num})', globals=globals(), number=1, repeat=3)

    time_info = pd.DataFrame({"id":scan_num, "pyteomics": pyteo, "pyopenms": pyopen, "pymzml": pymzml_data, 
                              "mzMLb":mzmlb, "MZA": mza, "mz5": mz5, "fun_type": "minimal"})
    all_timings.append(time_info)


pd.concat(all_timings).to_csv("data/singlefile_supp.csv", index=False)


# Rerun timings for pyopenms using the precompiled methods
pyop_exp = pyopenms.MSExperiment()
pyopenms.MzMLFile().load(f"{basename}.mzML", pyop_exp)

def get_chrom_pyop_precomp(mz, ppm):
    mzmin, mzmax = pmppm(mz, ppm)
    scan_dfs = []
    for spectrum in pyop_exp:
        if(spectrum.getMSLevel()==1):
            rt_val = spectrum.getRT()
            mz_vals, int_vals = spectrum.get_peaks()
            bet_idxs = (mzmin < mz_vals) & (mz_vals < mzmax)
            if(sum(bet_idxs)>0):
                df_scan = pd.DataFrame({'mz':mz_vals[bet_idxs], 'int':int_vals[bet_idxs], 'rt':[rt_val]*sum(bet_idxs)})
                scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
    
def get_chrom_pyop_precomp_2DPeak(mz, ppm):
    mzmin, mzmax = pmppm(mz, ppm)
    chrom_data=pyop_exp.get2DPeakDataLong(min_mz=mzmin, max_mz=mzmax, min_rt=0, max_rt=1e7)
    return(pd.DataFrame({"rt":chrom_data[0], "mz":chrom_data[1], "int":chrom_data[2]}))
    
def get_spec_pyop_precomp(scan_num):
    i = 0
    while i <= pyop_exp.size():
        if(int(pyop_exp[i].getNativeID().split("scan=")[-1].split()[0]) == scan_num):
            spec_data = pyop_exp[i].get_peaks()
            return(pd.DataFrame({"mz":spec_data[0], "int":spec_data[1]}))
        else:
            i += 1
    raise Exception(f"No scan number {scan_num} found")

def get_rtrange_pyop_precomp(rtstart, rtend):
    scan_dfs = []
    for spectrum in pyop_exp:
        if(spectrum.getMSLevel()==1):
            rt_val = spectrum.getRT()
            if(rtstart*60 < rt_val < rtend*60):
                mz_vals, int_vals = spectrum.get_peaks()
                df_scan = pd.DataFrame({'mz':mz_vals, 'int':int_vals, 'rt':[rt_val]*len(int_vals)})
                scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))
    
def get_rtrange_pyop_precomp_2DPeak(rtstart, rtend):
    rtrange_data=pyop_exp.get2DPeakDataLong(min_mz=0, max_mz=1e6, min_rt=rtstart*60, max_rt=rtend*60)
    return(pd.DataFrame({"rt":rtrange_data[0], "mz":rtrange_data[1], "int":rtrange_data[2]}))

def get_MS2premz_pyop_precomp(precursor_mz, ppm_acc):
    mzmin, mzmax = pmppm(precursor_mz, ppm_acc)
    scan_dfs = []
    for spectrum in pyop_exp:
        if(spectrum.getMSLevel()==2):
            premz_val = spectrum.getPrecursors()[0].getMZ()
            if(mzmin < premz_val < mzmax):
                rt_val = spectrum.getRT()/60
                spec_data = spectrum.get_peaks()
                df_scan = pd.DataFrame({"rt":rt_val, "premz":premz_val, "fragmz":spec_data[0], "int":spec_data[1]})
                scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))

def get_MS2fragmz_pyop_precomp(fragment_mz, ppm_acc):
    mzmin, mzmax = pmppm(fragment_mz, ppm_acc)
    scan_dfs = []
    for spectrum in pyop_exp:
        if(spectrum.getMSLevel()==2):
            premz_val = spectrum.getPrecursors()[0].getMZ()
            rt_val = spectrum.getRT()
            mz_vals, int_vals = spectrum.get_peaks()
            bet_idxs = (mzmin < mz_vals) & (mz_vals < mzmax)
            if(sum(bet_idxs)>0):
                df_scan = pd.DataFrame({"rt":rt_val, "premz":premz_val, "fragmz":mz_vals[bet_idxs], "int":int_vals[bet_idxs]})
                scan_dfs.append(df_scan)
    return(pd.concat(scan_dfs, ignore_index=True))



precomp_timings = [
    timeit.repeat(f"get_spec_pyop_precomp(9044)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_chrom_pyop_precomp(829.79730224, 5)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_chrom_pyop_precomp_2DPeak(829.79730224, 5)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_rtrange_pyop_precomp(11, 12)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_rtrange_pyop_precomp_2DPeak(11, 12)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_MS2premz_pyop_precomp(752.60633, 5)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_MS2fragmz_pyop_precomp(184.07304, 5)", globals=globals(), number=1, repeat=3)
]

init_timings = [
    timeit.repeat(f"get_spec_mzml_pyopenms('{basename}.mzML', 9044)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_chrom_mzml_pyopenms('{basename}.mzML', 829.79730224, 5)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_chrom_mzml_pyopenms_2DPeak('{basename}.mzML', 829.79730224, 5)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_rtrange_mzml_pyopenms('{basename}.mzML', 11, 12)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_rtrange_mzml_pyopenms_2DPeak('{basename}.mzML', 11, 12)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_MS2premz_mzml_pyopenms('{basename}.mzML', 752.60633, 5)", globals=globals(), number=1, repeat=3),
    timeit.repeat(f"get_MS2fragmz_mzml_pyopenms('{basename}.mzML', 184.07304, 5)", globals=globals(), number=1, repeat=3)
]

init_times = pd.DataFrame(init_timings, columns=["rep1", "rep2", "rep3"])
init_times["method"]=["spec", "chrom", "chrom_2D", "rtrange", "rtrange_2D", "premz", "fragmz"]
init_times["type"] = "init"
precomp_times = pd.DataFrame(precomp_timings, columns=["rep1", "rep2", "rep3"])
precomp_times["method"]=["spec", "chrom", "chrom_2D", "rtrange", "rtrange_2D", "premz", "fragmz"]
precomp_times["type"] = "precomp"

pd.concat([init_times, precomp_times]).to_csv("data/pyopenms_precomp.csv", index=False)