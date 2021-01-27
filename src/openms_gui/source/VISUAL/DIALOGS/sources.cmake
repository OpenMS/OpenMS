### the directory name
set(directory source/VISUAL/DIALOGS)

### list all filenames of the directory here
set(sources_list
DataFilterDialog.cpp
DataFilterDialog.ui
FeatureEditDialog.cpp
FeatureEditDialog.ui
HistogramDialog.cpp
LayerStatisticsDialog.cpp
LayerStatisticsDialog.ui
ListFilterDialog.cpp
ListFilterDialog.ui
PythonModuleRequirement.cpp
PythonModuleRequirement.ui
PythonSelector.cpp
PythonSelector.ui
SaveImageDialog.cpp
Plot1DGoToDialog.cpp
Plot1DGoToDialog.ui
Plot1DPrefDialog.cpp
Plot1DPrefDialog.ui
Plot2DGoToDialog.cpp
Plot2DGoToDialog.ui
Plot2DPrefDialog.cpp
Plot2DPrefDialog.ui
Plot3DPrefDialog.cpp
Plot3DPrefDialog.ui
SpectrumAlignmentDialog.cpp
SpectrumAlignmentDialog.ui
SwathTabWidget.cpp
SwathTabWidget.ui
TOPPASIOMappingDialog.cpp
TOPPASIOMappingDialog.ui
TOPPASInputFileDialog.cpp
TOPPASInputFileDialog.ui
TOPPASInputFilesDialog.cpp
TOPPASInputFilesDialog.ui
TOPPASOutputFilesDialog.cpp
TOPPASOutputFilesDialog.ui
TOPPASToolConfigDialog.cpp
TOPPASVertexNameDialog.cpp
TOPPASVertexNameDialog.ui
TOPPViewOpenDialog.cpp
TOPPViewOpenDialog.ui
TOPPViewPrefDialog.cpp
TOPPViewPrefDialog.ui
TheoreticalSpectrumGenerationDialog.cpp
TheoreticalSpectrumGenerationDialog.ui
ToolsDialog.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMSVisual_sources ${OpenMSVisual_sources} ${sources})

### source group definition
source_group("Source Files\\VISUAL\\DIALOGS" FILES ${sources})

