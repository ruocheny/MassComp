# MassComp

## About MassComp
MassComp is a loseless compresspr for mass spectrometry data. It compresses the mass-to-charge ratio and intensity pairs in mzXML files efficiently by calculating the hexadecimal difference of consecutive m/z values, and by searching for parts of the intensity values that match previous ones. 


## Getting Started
Download the full project.

## Run MassComp
The code can be run by visual studio 2012 on windows system.

Here's an example of this. Folder 'MSV000080896' is downloaded from MassIVE with id MSV000080896 and contains two mzXML files.

Run the executbale file MassComp in the project.
With the hint "please input the path of files to be compressing:", input the folder path "\MSV000080896\peak\Data_mzXML" to start compressing.

With the hint "please input the path of files to be decompressing:", input the folder path "\output\MSV000080896\peak\Data_mzXML" to start decompressing.


## Datasets
Datasets of mass spectrometry data can be downloaded from MassIVE https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp

## Contact
If you have any problem, please email me at rcyang624@126.com
