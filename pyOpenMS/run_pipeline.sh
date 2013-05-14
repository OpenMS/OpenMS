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

echo
ls -lh align_?.featureXML trafo_?.transformationXML
echo

echo
echo ============================================================
echo FeatureLinkerUnlabeledQT
echo ============================================================
echo
python pyTOPP/FeatureLinkerUnlabeledQT.py -write_ini linker.ini
python pyTOPP/FeatureLinkerUnlabeledQT.py -ini linker.ini -in align_1.featureXML,align_2.featureXML -out cons.consensusXML 

echo
ls -lh cons.consensusXML
echo

echo
echo ============================================================
echo IDMapper
echo ============================================================
echo
python pyTOPP/IDMapper.py -write_ini mapper.ini
python pyTOPP/IDMapper.py -ini mapper.ini -in align_1.featureXML -out cons_mapped.consensusXML -id /usr/share/OpenMS/examples/BSA/BSA1_OMSSA.idXML
python pyTOPP/IDMapper.py -ini mapper.ini -in cons.consensusXML -out cons_mapped.consensusXML -id /usr/share/OpenMS/examples/BSA/BSA1_OMSSA.idXML

echo
ls -lh cons_mapped.consensusXML
echo

rm pp.ini
rm picked_1.mzML
rm picked_2.mzML

rm fd.ini
rm feat_1.featureXML
rm feat_2.featureXML
rm feat_3.featureXML

rm align.ini
rm align_1.featureXML
rm align_2.featureXML
rm trafo_1.transformationXML
rm trafo_2.transformationXML

rm linker.ini
rm cons.consensusXML
