# --------------------------------------------------------------------------
# the directory name
set(directory source/EXAMPLES)

# --------------------------------------------------------------------------
# list all filenames of the directory here
set(executables_list
Tutorial_AASequence
Tutorial_Clustering
Tutorial_DPosition
Tutorial_DRange
Tutorial_Element
Tutorial_EmpiricalFormula
Tutorial_Enzyme
Tutorial_FeatureFinder
Tutorial_FeatureMap
Tutorial_FileIO
Tutorial_FileIO_Consumer
Tutorial_FileIO_mzML
Tutorial_GaussFilter
Tutorial_IdentificationClasses
Tutorial_Labeled
Tutorial_Logger
Tutorial_MapAlignment
Tutorial_MetaInfo
Tutorial_MorphologicalFilter
Tutorial_MSChromatogram
Tutorial_MSExperiment
Tutorial_MSSpectrum
Tutorial_Param
Tutorial_Precursor
Tutorial_RangeManager
Tutorial_Residue
Tutorial_ResidueModification
Tutorial_SavitzkyGolayFilter
Tutorial_Template           ## only used in developer_tutorial.doxygen, but still testing here if it compiles
Tutorial_TheoreticalSpectrumGenerator
Tutorial_TOPP
Tutorial_Unlabeled
)

# --------------------------------------------------------------------------
# pass source file list to the upper instance
set(EXAMPLES_executables ${EXAMPLES_executables} ${executables_list})

# --------------------------------------------------------------------------
set(executables_list
Tutorial_GUI_ListEditor
Tutorial_GUI_ParamEditor
Tutorial_GUI_Plot1D
)

# pass source file list to the upper instance
set(GUI_EXAMPLES_executables ${GUI_EXAMPLES_executables} ${executables_list})

# --------------------------------------------------------------------------
# add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${executables_list})
	list(APPEND sources_VS "${i}.cpp")
endforeach(i)
source_group("EXAMPLES" FILES ${sources_VS})
