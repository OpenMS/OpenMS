### the directory name
set(directory include/OpenMS/VISUAL/DIALOGS/UIC)

### list all filenames of the directory here
set(sources_list
DataFilterDialog.ui
FeatureEditDialog.ui
LayerStatisticsDialog.ui
Spectrum1DGoToDialog.ui
Spectrum1DPrefDialog.ui
Spectrum2DGoToDialog.ui
Spectrum2DPrefDialog.ui
Spectrum3DPrefDialog.ui
TOPPViewOpenDialog.ui
TOPPViewPrefDialog.ui
SpectrumAlignmentDialog.ui
TheoreticalSpectrumGenerationDialog.ui
TOPPASInputFileDialog.ui
TOPPASInputFilesDialog.ui
TOPPASOutputFilesDialog.ui
TOPPASIOMappingDialog.ui
TOPPASVertexNameDialog.ui
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### Apply UIC compiler
QT4_WRAP_UI_OWN(uiced_sources ${sources})

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${uiced_sources})

