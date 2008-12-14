### the directory name
set(directory source/EXAMPLES)

### list all filenames of the directory here
set(executables_list
Tutorial_AASequence
Tutorial_Clustering
Tutorial_ComparatorUtils
Tutorial_DPosition
Tutorial_DRange
Tutorial_Element
Tutorial_EmpiricalFormula
Tutorial_FeatureFinder
Tutorial_FeatureMap
Tutorial_FileIO
Tutorial_GaussFilter
Tutorial_InternalCalibration
Tutorial_Labeled
Tutorial_MSExperiment
Tutorial_MSSpectrum
Tutorial_MapAlignment
Tutorial_MetaInfo
Tutorial_Param
Tutorial_ParamEditor
Tutorial_PeakPickerCWT
Tutorial_RangeManager
Tutorial_Residue
Tutorial_SavitzkyGolayFilter
Tutorial_Spectrum1D
Tutorial_TOFCalibration
Tutorial_TheoreticalSpectrumGenerator
Tutorial_TopHatFilter
Tutorial_Unlabeled
Tutorial_typeAsString
)

### add path to the filenames
#set(executables)
#foreach(i ${executables_list})
#	list(APPEND executables ${directory}/${i})
#endforeach(i)

### pass source file list to the upper instance
set(EXAMPLES_executables ${EXAMPLES_executables} ${executables_list})

