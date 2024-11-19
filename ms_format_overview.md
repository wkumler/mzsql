## Overview of MS data formats and notes about their status

### mzML

Lots of different ways to get data out of mzMLs, including 3 dedicated Python libraries I've discovered:

  - pyteomics
  - pymzml
  - pyopenms: Wrapper around a bunch of C code to do things
      - Has multiple methods of accessing mz/rt chunks

I can't seem to get the indexed version working, at least in such a way that speeds up spectrum access. It feels like something along the lines of the below code should work (for pyteomics), but returns strange errors. 

```python
file_data=mzml.PreIndexedMzML("demo_data/180205_Poo_TruePoo_Full1_idx.mzML")
file_data=file_data.build_byte_index()
```

I also can't tell if something similar is implemented in the other libraries. Also, it may be that the other libraries just handle indexing automagically.

### mzMLb

Seems to only be implemented in pyteomics? Double check this later.

### mzDB

Creation is pretty straightforward via the raw2mzDB.exe tool but we can't really figure out how to access the data. Package code exists (the original developers have a Python (and R) port of some Rust code but have responded saying that the [rt/mz range extraction isn't yet supported](https://github.com/mzdb/mzdb-rs/issues/3) and hasn't been [since Johannes asked about it in 2018?](https://github.com/mzdb/rmzdb/issues/3). Developer seems responsive but just overwhelmed. R-specific rmzdb package (not to be confused with the rmzdb port of the Rust code) hasn't been touched in 5 years, the Rust/Python/R hasn't been updated in one and a half.

Someone else [wrote a Python package for mzDB access](https://github.com/jerkos/pymzdb) (annoyingly, named the same as the dev Rust port) which has all the associated functions but it's 9 years old and I haven't gotten the chance to give it a try yet and the README is sparse.

Dev version of pymzdb: 8 commits total, 2 years ago

Other version of pymzdb: 10 commits total, 9 years ago

### MzTree

No idea where to get started with this one because the README is just documentation of a web server API that I can't figure out how to query. Admittedly, haven't tried installing Java and following the instructions yet but I just have no idea what the inputs should be. mzML? .raw files? Some other format?
In the sole Github issue that exists for this, the developer says "MZTree is not a format like, say, mzML or .raw. It is a storage and retrieval system." Frankly, I don't understand what this means - the data's gotta be stored somewhere and the article certainly makes it sound like it's a method for storing the data for rapid access.
Pythonic access to this will likely be via local server requests but it at least looks like that should be doable.

8 commits in the repo, last one 7 years ago.

### mz5

Successfully created mz5 file via msconvert and wrote Python code for accessing mz5 files but it's pretty basic. Could maybe be sped up/accelerated with the [pymz5](https://github.com/jmchilton/pymz5/tree/master) library but that library is pretty old now and didn't seem to take off at all (zero issues, very few interactions). Doesn't provide any info about how to use the library either beyond install even though it seems like a ton of work was put into it (782 commits).
I still need to try this library and see whether it's functional but I'm skeptical - 11 years ago was still Python 2, right?

782 commits in the repo, last one 11 years ago

### MZA

Files were created successfully via the command-line tool (although the README was only recently updated to tell the user how to download it from the releases) but are inaccessible via the mzapy package because it currently only supports TOF files. Branches for Orbitrap are under development but can't be currently installed. Wrote my own Python access code but it just seems to be a badly-organized hdf5 file format since I can't take advantage of their custom functions / structures.

20 commits in the repo, last one 2 days ago (after I opened an issue for it, prior to which it had been 2 years)

### Aird

Haven't started yet. Lots of development, lots of support and interest. Seems to have built-in Python support, last updated 2 years ago. Unclear whether it's faster or just more compressed - they don't really show a speed vs format graph AFAIK. For whatever reason, Python isn't recognizing the install but I think that's on me. Conversion was successful and I've uploaded the .aird and .json files to demo_data.

466 commits in the repo, last one 7 months ago (2 years for the python-specific stuff)



## Incomparables

### YAFMS

Seems fully deprecated - the link in the original manuscript now just redirects to PNNL's general software page and I can't find any residuals of the code via quick Googling.

### mzMD

Feels like a strange port of MzTree, which makes sense given that that's what it's derived from. Their Github is just a pared-down version of MzTree's, no install or startup instructions (maybe the exact same as MzTree's?). Kinda surprised that the Bioinformatics manuscript made it past review and editing because it's filled with grammatical errors. Also uses a super weird definition of "peak" as "tuple of three parameters that are mass-to-charge (m/z) value, retention time and intensity".

3 commits in the repo (basically just "add files via upload"), last one 4 years ago.

### mzRTree

Haven't started yet. Seems like an older file type and only really compares to GUI formats. No idea how to convert the data into this format or access it - there's no Github I've been able to find and the only code associated I have found [is in Matlab](https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fwww.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fe5af2033-4a80-11e4-9553-005056977bd0%2F65b8e141-bbf2-ecf1-80a0-1929890c734e%2Ffiles%2F3DSpectra%2FmzRTreeCreation.m&embed=web). mzDB somehow made it work (and found an example file?) but I have no idea how they did this.

### Toffee

Only works with TOF data, as far as I can tell. That disqualifies it from comparison in my opinion.

