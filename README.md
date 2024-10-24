## Data and code demonstrating the efficacy of databases relative to existing mass-spectrometry data formats

This repository pioneers the use of SQL databases for the efficient storage and access of raw
mass-spectrometry (MS) data. Here, "raw" MS data refers to the data encoded in vendor or .mzML files and contains
coordinates of retention time, *m/z* ratio, and intensity.

Existing data formats fail to provide intuitive, rapid, and programmatic search of raw MS data and require learning
the quirks and conceits of idiosyncratic file formats. Databases have
consistent and transferrable syntax for queries, indexing and binary search for rapid data extraction,
benefit from multi-file aggregations of data, can store processed data alongside, and are established as an industry standard.

### Consistent and transferrable syntax

SQL is widely known and understood. Transforming MS data into a format that can be queried using SQL allows
developers to focus on building downstream products instead of comprehending the particulars of a given
file type. SQL also often results in queries much closer to natural language - compare 

```SELECT * FROM MS1 WHERE rt BETWEEN 5 AND 6```

to

```import matplotlib.pyplot as plt
from pyteomics import mzml
import numpy as np

mz_values = []
retention_times = []
int_values = []
with mzml.MzML(file) as reader:
    for spectrum in reader:
        if spectrum['ms level'] == 1:
            mz_values.extend(spectrum['m/z array'])
            mz_values.extend(spectrum['intensity array'])
            rt = spectrum['scanList']['scan'][0]['scan start time']
            retention_times.extend([rt] * len(spectrum['m/z array']))
```

### Indexing and binary search for rapid data extraction

Databases can be highly optimized for the extraction of data subsets. For MS data, this means that
it's just as easy to find a single scan as it is to extract an entire chromatogram. Many existing
data formats have indexes built in that identify the start or end byte of a given scan, but still
require every scan to be extracted when the information *within* the scan is of interest (as when
subsetting by *m/z* or intensity).

### Benefit from multi-file data aggregations

Existing MS data formats preserve the idea of "one sample = one file". With databases, the filename or
identifier can be simply encoded as an additional column and the entire dataset can be concatenated
together, allowing for rapid searches across every file simultaneously. Without this, MS datasets
increase linearly in time with the number of files in the dataset as each file is looped over
one at a time.

### Can store processed data alongside

After chromatographic peakpicking or compound identification has been performed, the resulting
peak list (start_rt, end_rt, mean_mz, area, etc.) can be stored as an additional table in the
database. Not only does this make it easy to share all the data (often a single .sql file or 
database connection instead of multiple mzMLs and CSVs), but it also allows for powerful
join operations between the raw and processed data that enable rapid extraction of the
individual data points that compose a chromatographic peak. Additionally, metadata or corrected
retention times similarly benefit from simple and powerful joins that minimize memory usage.

### Established as an industry standard

Databases are being actively developed by corporations much larger than any mass-spectrometry lab. This
means that the improvements in access speed and compression efficiency can be passed along to the MS
data within them without any extra effort. Additionally, databases are supported by an enormous community
of experts rather than individual developers who may leave the project and result in a deprecated file type.

## Goals

1. Determine when and how databases offer a performance improvement over existing MS file types
2. Explore multiple database types and provide use-case recommendations for each
3. Provide code that allows for the elegant conversion of existing MS data into database

## Existing formats for comparison:

| Format | Publication | Year | Interface | In msconvert? | m/z values searchable? | Written/examples in | Notes |
| --- | --- | --- | --- | --- | --- | --- | ---
| mzML | [Link](https://www.mcponline.org/article/S1535-9476(20)31387-6/fulltext) | 2010 | XML | Yes | No | Many | |
| mzRTree | [Link](https://dx.doi.org/10.1016/j.jprot.2010.02.006) | 2010 | Tree | No | ? | Java | |
| YAFMS | [Link](https://dx.doi.org/10.1016/j.jasms.2010.06.014) | 2010 | SQL | No | No | C# | Deprecated? |
| mz5 | [Link](https://dx.doi.org/10.1074/mcp.O111.011379) | 2012 | HDF5 | Yes | ? | ? | |
| mzDB | [Link](https://dx.doi.org/10.1074/mcp.O114.039115) | 2015 | SQL | No | Yes? | Java, C++ | |
| Indexed mzML | [Link](https://dx.doi.org/10.1371/journal.pone.0125108) | 2015 | XML | Yes | No | C++, Python | |
| MzTree | [Link](https://dx.doi.org/10.1371/journal.pone.0188059) | 2017 | SQLite, Tree | No | Yes? | Java | |
| Toffee | [Link](https://dx.doi.org/10.1038/s41598-020-65015-y) | 2020 | HDF5(?) | No | Yes? | Python | Only for TOF data? |
| mzMLb | [Link](https://dx.doi.org/10.1021/acs.jproteome.0c00192) | 2021 | XML | Yes | No | Python | |
| Aird | [Link](https://dx.doi.org/10.1186/s12859-021-04490-0) | 2022 | Unique | No | ? | C# | |
| mzMD | [Link](https://dx.doi.org/10.1093/bioinformatics/btac098) | 2022 | SQL | No | ? | Java | |
| MZA | [Link](https://dx.doi.org/10.1021/acs.jproteome.2c00313) | 2023 | HDF5 | No | Yes? | Python | Separate Python package [here](https://dx.doi.org/10.1021/acs.analchem.3c01653)

Not mentioned above: mzData, mzXML, Shaduf, imzML, others?

