#!/bin/sh

echo
echo ============================================================
echo PeakPickerHiRes.py
echo ============================================================
echo
python pyTOPP/PeakPickerHiRes.py -write_ini pp.ini
python pyTOPP/PeakPickerHiRes.py -ini pp.ini -in ../share/OpenMS/examples/peakpicker_tutorial_1_baseline_and_noise_filtered.mzML  -out picked_1.mzML
python pyTOPP/PeakPickerHiRes.py -ini pp.ini -in ../share/OpenMS/examples/peakpicker_tutorial_2.mzML -out picked_2.mzML
python pyTOPP/PeakPickerHiRes.py -ini pp.ini -in ../share/OpenMS/examples/LCMS-centroided.mzML -out picked_3.mzML

echo
ls -lh picked_?.mzML
echo


echo
echo ============================================================
echo FeatureFinderCentroided.py
echo ============================================================
echo
python pyTOPP/FeatureFinderCentroided.py -write_ini fd.ini
python pyTOPP/FeatureFinderCentroided.py -ini fd.ini -in ../share/OpenMS/examples/LCMS-centroided.mzML -out feat_1.featureXML
python pyTOPP/FeatureFinderCentroided.py -ini fd.ini -in ../share/OpenMS/examples/LCMS-centroided.mzML -out feat_2.featureXML
python pyTOPP/FeatureFinderCentroided.py -ini fd.ini -in picked_3.mzML -out feat_3.featureXML

echo
ls -lh feat_?.featureXML
echo

echo
echo ============================================================
echo MapAlignerPoseClustering.py
echo ============================================================
echo
python pyTOPP/MapAlignerPoseClustering.py -write_ini align.ini
python pyTOPP/MapAlignerPoseClustering.py -ini align.ini -in feat_1.featureXML,feat_2.featureXML -out align_1.featureXML,align_2.featureXML -trafo_out trafo_1.transformationXML,trafo_2.transformationXML

#rm pp.ini
#rm picked_1.mzML
#rm picked_2.mzML

#rm fd.ini
#rm feat1.featureXML
#rm feat2.featureXML

#rm align.ini
#rm align1.featureXML
#rm align2.featureXML
#rm trafo1.featureXML
#rm trafo2.featureXML
