### the directory name
set(directory source/VISUAL)

### list all filenames of the directory here
set(sources_list
AxisPainter.cpp
AxisTickCalculator.cpp
AxisWidget.cpp
ColorSelector.cpp
EnhancedTabBar.cpp
EnhancedTabBarWidgetInterface.cpp
EnhancedWorkspace.cpp
FilterList.cpp
FilterList.ui
GUIProgressLoggerImpl.cpp
HistogramWidget.cpp
InputFile.cpp
InputFile.ui
InputFileList.cpp
InputFileList.ui
LayerListView.cpp
LayerData.cpp
ListEditor.cpp
LogWindow.cpp
MetaDataBrowser.cpp
MultiGradient.cpp
MultiGradientSelector.cpp
OutputDirectory.cpp
OutputDirectory.ui
ParamEditor.cpp
ParamEditor.ui
RecentFilesMenu.cpp
SpectraIdentificationViewWidget.cpp
SpectraViewWidget.cpp
SpectraSelectionTabs.cpp
Spectrum1DCanvas.cpp
Spectrum1DWidget.cpp
Spectrum2DCanvas.cpp
Spectrum2DWidget.cpp
Spectrum3DCanvas.cpp
Spectrum3DOpenGLCanvas.cpp
Spectrum3DWidget.cpp
SpectrumCanvas.cpp
SpectrumWidget.cpp
SwathLibraryStats.cpp
SwathLibraryStats.ui
TOPPASEdge.cpp
TOPPASInputFileListVertex.cpp
TOPPASMergerVertex.cpp
TOPPASOutputFileListVertex.cpp
TOPPASResource.cpp
TOPPASResources.cpp
TOPPASScene.cpp
TOPPASSplitterVertex.cpp
TOPPASToolVertex.cpp
TOPPASTreeView.cpp
TOPPASVertex.cpp
TOPPASWidget.cpp
TOPPViewIdentificationViewBehavior.cpp
TOPPViewMenu.cpp
#TOPPViewPreferences.cpp
TOPPViewSpectraViewBehavior.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMSVisual_sources ${OpenMSVisual_sources} ${sources})

### source group definition
source_group("Source Files\\VISUAL" FILES ${sources})


### icons
# add   : icons are added to source/VISUAL/ICONS/resources.qrc
# remove: after removing an icon, you have to rerun 'cmake' to fix the dependencies
QT5_ADD_RESOURCES(qt_resource_file source/VISUAL/ICONS/resources.qrc)
set(OpenMSVisual_sources ${OpenMSVisual_sources} ${qt_resource_file})
set_property(SOURCE ${qt_resource_file} PROPERTY SKIP_AUTOGEN ON)
