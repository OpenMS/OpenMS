#!/bin/sh
python pyTOPP/PeakPickerHiRes.py -write_ini pp.ini
python pyTOPP/PeakPickerHiRes.py -ini pp.ini -in example1.mzXML -out ex1picked.mzML
python pyTOPP/PeakPickerHiRes.py -ini pp.ini -in example2.mzXML -out ex2picked.mzML

python pyTOPP/FeatureFinderCentroided.py -write_ini fd.ini
python pyTOPP/FeatureFinderCentroided.py -ini fd.ini -in example1.mzXML -out feat1.featureXML
python pyTOPP/FeatureFinderCentroided.py -ini fd.ini -in example2.mzXML -out feat2.featureXML

python pyTOPP/MapAlignerPoseClustering.py -write_ini align.ini
python pyTOPP/MapAlignerPoseClustering.py -ini align.ini -in feat1.featureXML,feat2.featureXML -out align1.featureXML,align2.featureXML -trafo_out trafo1.transformationXML,trafo2.transformationXML
