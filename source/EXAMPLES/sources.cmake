### the directory name
set(directory source/EXAMPLES)

### list all filenames of the directory here
set(sources_list
Tutorial_AASequence.C
Tutorial_Clustering.C
Tutorial_ComparatorUtils.C
Tutorial_DPosition.C
Tutorial_DRange.C
Tutorial_Element.C
Tutorial_EmpiricalFormula.C
Tutorial_FeatureFinder.C
Tutorial_FeatureMap.C
Tutorial_FileIO.C
Tutorial_GaussFilter.C
Tutorial_InternalCalibration.C
Tutorial_Labeled.C
Tutorial_MSExperiment.C
Tutorial_MSSpectrum.C
Tutorial_MapAlignment.C
Tutorial_MetaInfo.C
Tutorial_Param.C
Tutorial_ParamEditor.C
Tutorial_PeakPickerCWT.C
Tutorial_RangeManager.C
Tutorial_Residue.C
Tutorial_SavitzkyGolayFilter.C
Tutorial_Spectrum1D.C
Tutorial_TOFCalibration.C
Tutorial_TheoreticalSpectrumGenerator.C
Tutorial_TopHatFilter.C
Tutorial_Unlabeled.C
Tutorial_typeAsString.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### source group definition
source_group("Source Files\\EXAMPLES" FILES ${sources})

