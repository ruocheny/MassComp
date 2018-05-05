# MassComp

## About MassComp
MassComp is a loseless compressor for mass spectrometry data. It compresses the mass-to-charge ratio and intensity pairs in mzXML files efficiently by calculating the hexadecimal difference of consecutive m/z values, and by searching for parts of the intensity values that match previous ones. The remaining parts of the mzXML (e.g., metadata associated to the experiments) is compressed with the general compression algorithm 7zip. 


## Getting Started
Download the full project.

## Run MassComp
### Linux system
To compile:

`g++ -o masscomp MassComp.cpp tinyxml2.cpp`

To compress:

`./masscomp -c fileOri.mzXML fileMasscomp`

To decompress

`./masscomp -d fileMasscomp fileDecomp.mzXML` for single precision
`./masscomp -d fileMasscomp fileDecomp.mzXML -64` for double precision

### Windows system
Current implementation of the code can be run by visual studio on windows system (we will provide a Linux and Mac implementation shortly).

Here's an example of this. Folder 'MSV000080896' is downloaded from MassIVE with id MSV000080896 and contains two mzXML files.

Run the executbale file MassComp in the project.
With the hint "please input the path of files to be compressing:", input the folder path "\MSV000080896\peak\Data_mzXML" to start compressing. The output folder need an external compression with 7zip.

With the hint "please input the path of files to be decompressing:", input the folder path "\output\MSV000080896\peak\Data_mzXML" to start decompressing. Before start the decompressing, the zip file needs an external decompression with 7zip.


## Datasets
Datasets of mass spectrometry data can be downloaded from MassIVE https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp

## Authors
MassComp was created by Ruochen Yang and Idoia Ochoa at University of Illinois at Urbana-Champaign.

## Contact
If you have any problem, please email me at rcyang624@126.com
