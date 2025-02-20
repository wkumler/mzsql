import requests
from requests.exceptions import ConnectionError
import pandas as pd
from .helpers import pmppm

def get_chrom_mzMD(url_port, mz, ppm):
    """
    Retrieves chromatogram data from a server for a specific m/z range.

    Args:
        url_port (str): The URL and port of the server, e.g., 'http://localhost:5000'.
        mz (float): The target m/z value.
        ppm (float): The mass tolerance in parts per million (ppm).

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified m/z range.
    """
    mz_min, mz_max = pmppm(mz, ppm)
    request_string = f"{url_port}/api/v2/getpoints?mzmin={mz_min}&mzmax={mz_max}&rtmin=0&rtmax=1000000&n=0&m=0"
    
    try:
        response = requests.get(request_string)
        response.raise_for_status()
        if response.status_code == 204:
            print("No content returned. Have you loaded a file?")
    except ConnectionError as e:
        if "Connection refused" in str(e):
            print("Connection refused. Have you started the server?")
        else:
            print(f"Connection error: {e}")

    chrom_data = pd.DataFrame(response.json(), columns=["pointId", "traceId", "mz", "rt", "intensity"])
    chrom_data = chrom_data[["rt", "mz", "intensity"]].sort_values(by="rt")
    chrom_data.columns = ["rt", "mz", "int"]
    return chrom_data

def get_rtrange_mzMD(url_port, rtstart, rtend):
    """
    Retrieves chromatogram data from a server for a specific retention time range.

    Args:
        url_port (str): The URL and port of the server, e.g., 'http://localhost:5000'.
        rtstart (float): The start retention time value.
        rtend (float): The end retention time value.

    Returns:
        pd.DataFrame: A DataFrame containing retention time (rt), m/z, and intensity data 
                      within the specified retention time range.
    """
    request_string = f"{url_port}/api/v2/getpoints?mzmin=0&mzmax=1000000&rtmin={rtstart}&rtmax={rtend}&n=0&m=0"
    
    try:
        response = requests.get(request_string)
        response.raise_for_status()
        if response.status_code == 204:
            print("No content returned. Have you loaded a file?")
    except ConnectionError as e:
        if "Connection refused" in str(e):
            print("Connection refused. Have you started the server?")
        else:
            print(f"Connection error: {e}")

    chrom_data = pd.DataFrame(response.json(), columns=["pointId", "traceId", "mz", "rt", "intensity"])
    chrom_data = chrom_data[["rt", "mz", "intensity"]].sort_values(by="rt")
    chrom_data.columns = ["rt", "mz", "int"]
    return chrom_data