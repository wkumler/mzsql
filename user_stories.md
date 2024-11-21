User Stories (For Class Purposes)

1. Trish is a student, using mass spectrometry data for her capstone project. She is wants to use this data for _____. She is looking for an easy way to sort through large quantites of data from a variety of file formats. She specifically does not want to deal with keeping track off all the 
differences between file formats. She has a limited knowledge of programming, with only a couple of classes years ago.

2. Karl is a researcher in medicine. He uses mass spectrometry data to find new drugs. He wants a fast and consistant method of filtering data, with an emphasis on forward compatibility. Karl is interested in making it easier to work with others, who work with different file formats. He has extensive knowledge both in programming and in the use of existing file formats.

3. Laura is a professor of an environmental omics lab. She's accumulated years of MS data and wants to scan through all of it quickly and concisely but the current storage format is slow and memory-intensive to scan through. She's got some exposure to R and Python from her students but doesn't want to learn the details of any particular language.

4. Ethan is brand new to mass spectrometry and can't figure out how the things he's learning in MS courses is encoded in the existing mzML file types. He wants to see the variance in m/z values for a known compound over the course of a mass spec run and make chromatograms but doesn't know where the mz/rt/int matrix is or how to read it into Python so he can use his prior experience with the language to make pretty plots.

5. Justine is a software developer for mass spectrometry instruments. She wants to develop visuals for her clients that show the differences between profile mode data and multiple centroiding algorithms, but the existing software doesn't provide a clear way of plotting m/z across retention time.

6. Lana is a reviewer for a paper that claims to have discovered a new compound in the environment. She would like to be able to quickly and easily view a chromatogram for the compound alongside a plot of the MS2 fragments that were collected, ideally without needing to download many packages or lots of data.

7. Dale is a new graduate student who is trying to figure out how to interface with a new MS2 library. His data's coming off the instrument in proprietary formats but he can't figure out how to quickly and easily extract the fragments associated with the features he got from MS-DIAL so that they can be correctly formatted for upload to the web server and matching against the library.

8. Mary is a lab manager who's concerned about the new MS instrument that's arriving shortly because it's known for generating really large data files. She would like a robust file format that's smaller on disk than mzML but accessible in a way that vendor file types aren't, ideally without requiring loading the entire massive file all at once.


Use Cases
1.
User: inputs file for conversion.
System: Converts file to new sql dataframe style.
Implied: Ensures that the file is not corrupt, and follows the expected format 
Implied: has a user interface to load/interact files

2.
User: inputs request for specific data section. I.e. outputs above a certain magnitude only.
System: Computes and returns the selection of outputs that meet the requirement.
Implied: 
Implied: 




Component specification

  - Database!
    - MS1 table
    - MS2 table
  - Manuscript draft?
  - Vignettes?

