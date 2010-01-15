### the directory name
set(directory source/VISUAL/DIALOGS)

### list all filenames of the directory here
set(sources_list
DemoDialog.C
DBOpenDialog.C
DataFilterDialog.C
FeatureEditDialog.C
HistogramDialog.C
LayerStatisticsDialog.C
SaveImageDialog.C
Spectrum1DGoToDialog.C
Spectrum1DPrefDialog.C
Spectrum2DGoToDialog.C
Spectrum2DPrefDialog.C
Spectrum3DPrefDialog.C
TOPPViewOpenDialog.C
TOPPViewPrefDialog.C
ToolsDialog.C
TheoreticalSpectrumGenerationDialog.C
SpectrumAlignmentDialog.C
TOPPASInputFileDialog.C
TOPPASInputFilesDialog.C
TOPPASOutputFilesDialog.C
TOPPASToolConfigDialog.C
TOPPASIOMappingDialog.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\VISUAL\\DIALOGS" FILES ${sources})

