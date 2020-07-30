### the directory name
set(directory include/OpenMS/VISUAL/DIALOGS)

### list all header files of the directory here
set(sources_list_h
DataFilterDialog.h
FeatureEditDialog.h
HistogramDialog.h
LayerStatisticsDialog.h
PythonModuleRequirement.h
PythonSelector.h
SaveImageDialog.h
Spectrum1DGoToDialog.h
Spectrum1DPrefDialog.h
Spectrum2DGoToDialog.h
Spectrum2DPrefDialog.h
Spectrum3DPrefDialog.h
SpectrumAlignmentDialog.h
SwathTabWidget.h
TOPPASIOMappingDialog.h
TOPPASInputFileDialog.h
TOPPASInputFilesDialog.h
TOPPASOutputFilesDialog.h
TOPPASToolConfigDialog.h
TOPPASVertexNameDialog.h
TOPPViewOpenDialog.h
TOPPViewPrefDialog.h
TheoreticalSpectrumGenerationDialog.h
ToolsDialog.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### treat as source files, for autoMOC'ing instead of manually calling QT5_WRAP_CPP()
set(OpenMSVisual_sources ${OpenMSVisual_sources} ${sources_h})
### pass header file list to the upper instance
set(OpenMSVisual_sources_h ${OpenMSVisual_sources_h} ${sources_h})

### header group definition for IDE's
source_group("Header Files\\OpenMS\\VISUAL\\DIALOGS" FILES ${sources_h})

