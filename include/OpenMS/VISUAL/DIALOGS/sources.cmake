### the directory name
set(directory include/OpenMS/VISUAL/DIALOGS)

### list all MOC filenames of the directory here
set(sources_list
DemoDialog.h
DBOpenDialog.h
DataFilterDialog.h
FeatureEditDialog.h
HistogramDialog.h
LayerStatisticsDialog.h
SaveImageDialog.h
Spectrum1DGoToDialog.h
Spectrum1DPrefDialog.h
Spectrum2DGoToDialog.h
Spectrum2DPrefDialog.h
Spectrum3DPrefDialog.h
TOPPViewOpenDialog.h
TOPPViewPrefDialog.h
ToolsDialog.h
SpectrumAlignmentDialog.h
TheoreticalSpectrumGenerationDialog.h
TOPPASInputFileDialog.h
TOPPASInputFilesDialog.h
TOPPASOutputFilesDialog.h
TOPPASToolConfigDialog.h
TOPPASIOMappingDialog.h
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### Apply MOC compiler
QT4_WRAP_CPP(mocced_sources ${sources})

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${mocced_sources})

source_group("Source Files\\OpenMS\\VISUAL\\DIALOGS" FILES ${mocced_sources})

### list all header files of the directory here
set(sources_list_h
DemoDialog.h
DBOpenDialog.h
DataFilterDialog.h
FeatureEditDialog.h
HistogramDialog.h
LayerStatisticsDialog.h
SaveImageDialog.h
Spectrum1DGoToDialog.h
Spectrum1DPrefDialog.h
Spectrum2DGoToDialog.h
Spectrum2DPrefDialog.h
Spectrum3DPrefDialog.h
TOPPViewOpenDialog.h
TOPPViewPrefDialog.h
ToolsDialog.h
SpectrumAlignmentDialog.h
TheoreticalSpectrumGenerationDialog.h
TOPPASInputFileDialog.h
TOPPASInputFilesDialog.h
TOPPASOutputFilesDialog.h
TOPPASToolConfigDialog.h
TOPPASIOMappingDialog.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\VISUAL\\DIALOGS" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

