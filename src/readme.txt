//compile
g++ -o masscomp MassComp.cpp tinyxml2.cpp

// to compress
./masscomp -c fileOri.mzXML fileMasscomp

// to decompress
./masscomp -d fileMasscomp fileDecomp.mzXML

// to compare the m/z-int pairs
./masscomp -cmp fileOri.mzXML fileDecomp.mzXML
