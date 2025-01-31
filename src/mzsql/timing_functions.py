import time
import timeit
from mzsql import *

def time_chrom(fun_suffix, file_ending, target_mz, ppm_acc, verbose=True):
    if(fun_suffix in ["mztree", "mzMD"]):
        rep_function = f"get_chrom_{fun_suffix}('{file_ending}', {target_mz}, {ppm_acc})"
    else:
        rep_function = f"get_chrom_{fun_suffix}('E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01{file_ending}', {target_mz}, {ppm_acc})"
    if(verbose):
        start_time = time.time()
        print(f"Running {fun_suffix}", end="\r")
    time_vals = timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
    if(verbose):
        elapsed = time.time() - start_time
        print(f"Running {fun_suffix}... ({elapsed:.2f}s)")
    return(time_vals)

def time_spec(fun_suffix, file_ending, spec_id, verbose=True):
    rep_function = f"get_spec_{fun_suffix}('E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01{file_ending}', {spec_id})"
    if(verbose):
        print(f"Running {fun_suffix}", end="\r")
        start_time = time.time()
    time_vals = timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
    if(verbose):
        elapsed = time.time() - start_time
        print(f"Running {fun_suffix}... ({elapsed:.2f}s)")
    return(time_vals)

def time_rtrange(fun_suffix, file_ending, rt, verbose=True):
    if(fun_suffix in ["mztree", "mzMD"]):
        rep_function = f"get_rtrange_{fun_suffix}('{file_ending}', {rt-0.2}, {rt+0.2})"
    else:
        rep_function = f"get_rtrange_{fun_suffix}('E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01{file_ending}', {rt-0.2}, {rt+0.2})"
    if(verbose):
        print(f"Running {fun_suffix}", end="\r")
        start_time = time.time()
    time_vals = timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
    if(verbose):
        elapsed = time.time() - start_time
        print(f"Running {fun_suffix}... ({elapsed:.2f}s)")
    return(time_vals)

def time_premz(fun_suffix, file_ending, target_premz, ppm_acc, verbose=True):
    rep_function = f"get_MS2premz_{fun_suffix}('E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01{file_ending}', {target_premz}, {ppm_acc})"
    if(verbose):
        start_time = time.time()
        print(f"Running {fun_suffix}", end="\r")
    time_vals = timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
    if(verbose):
        elapsed = time.time() - start_time
        print(f"Running {fun_suffix}... ({elapsed:.2f}s)")
    return(time_vals)

def time_fragmz(fun_suffix, file_ending, target_fragmz, ppm_acc, verbose=True):
    rep_function = f"get_MS2fragmz_{fun_suffix}('E:/mzsql/MTBLS10066/20220921_LEAP-POS_BL01{file_ending}', {target_fragmz}, {ppm_acc})"
    if(verbose):
        start_time = time.time()
        print(f"Running {fun_suffix}", end="\r")
    time_vals = timeit.repeat(rep_function, globals=globals(), number=1, repeat=1)
    if(verbose):
        elapsed = time.time() - start_time
        print(f"Running {fun_suffix}... ({elapsed:.2f}s)")
    return(time_vals)
