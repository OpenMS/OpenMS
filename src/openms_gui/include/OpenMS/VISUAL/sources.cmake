### the directory name
set(directory include/OpenMS/VISUAL)

### list all header files of the directory here
set(sources_list_h
AxisTickCalculator.h
AxisPainter.h
AxisWidget.h
ColorSelector.h
EnhancedTabBar.h
EnhancedTabBarWidgetInterface.h
EnhancedWorkspace.h
FilterList.h
GUIProgressLoggerImpl.h
HistogramWidget.h
InputFile.h
InputFileList.h
LayerListView.h
LayerData.h
ListEditor.h
LogWindow.h
MetaDataBrowser.h
MultiGradient.h
MultiGradientSelector.h
OutputDirectory.h
ParamEditor.h
RecentFilesMenu.h
SpectraViewWidget.h
SpectraIdentificationViewWidget.h
SpectraSelectionTabs.h
Spectrum1DCanvas.h
Spectrum1DWidget.h
Spectrum2DCanvas.h
Spectrum2DWidget.h
Spectrum3DCanvas.h
Spectrum3DOpenGLCanvas.h
Spectrum3DWidget.h
SpectrumCanvas.h
SpectrumWidget.h
SwathLibraryStats.h
TOPPASEdge.h
TOPPASInputFileListVertex.h
TOPPASMergerVertex.h
TOPPASOutputFileListVertex.h
TOPPASResource.h
TOPPASResources.h
TOPPASScene.h
TOPPASSplitterVertex.h
TOPPASToolVertex.h
TOPPASTreeView.h
TOPPASVertex.h
TOPPASWidget.h
TOPPViewIdentificationViewBehavior.h
TOPPViewMenu.h
#TOPPViewPreferences.h
TOPPViewSpectraViewBehavior.h
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
source_group("Header Files\\OpenMS\\VISUAL" FILES ${sources_h})

