## Data and code demonstrating the efficacy of databases relative to existing mass-spectrometry data formats

This repository pioneers the use of SQL databases for the efficient storage and access of raw
mass-spectrometry (MS) data. Here, "raw" MS data refers to the data encoded in vendor or .mzML files and contains
coordinates of retention time, *m/z* ratio, and intensity.

Existing data formats fail to provide intuitive, rapid, and programmatic search of raw MS data and require learning
the quirks and conceits of idiosyncratic file formats. Databases have
consistent and transferrable syntax for queries, indexing and binary search for rapid data extraction,
benefit from multi-file aggregations of data, and are established as an industry standard.

### Consistent and transferrable syntax

SQL is widely known and understood. Transforming MS data into a format that can be queried using SQL allows
developers to focus on building downstream products instead of comprehending the particulars of a given
file type.

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

### Established as an industry standard

Databases are being actively developed by corporations much larger than any mass-spectrometry lab. This
means that the improvements in access speed and compression efficiency can be passed along to the MS
data within them without any extra effort. Additionally, databases are supported by an enormous community
of experts rather than individual developers who may leave the project and result in a deprecated file type.

## Goals

1. Determine when and how databases offer a performance improvement over existing MS file types
2. Explore multiple database types and provide use-case recommendations for each
3. Provide code that allows for the elegant conversion of existing MS data into database
