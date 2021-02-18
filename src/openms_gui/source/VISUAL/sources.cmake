### the directory name
set(directory source/VISUAL)

### list all filenames of the directory here
set(sources_list
AxisPainter.cpp
AxisTickCalculator.cpp
AxisWidget.cpp
ColorSelector.cpp
DataSelectionTabs.cpp
DIATreeTab.cpp
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
Plot1DCanvas.cpp
Plot1DWidget.cpp
Plot2DCanvas.cpp
Plot2DWidget.cpp
Plot3DCanvas.cpp
Plot3DOpenGLCanvas.cpp
Plot3DWidget.cpp
PlotCanvas.cpp
PlotWidget.cpp
RecentFilesMenu.cpp
SpectraIDViewTab.cpp
SpectraTreeTab.cpp
SwathLibraryStats.cpp
SwathLibraryStats.ui
TableView.cpp
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
TOPPViewMenu.cpp
TreeView.cpp
TVDIATreeTabController.cpp
TVIdentificationViewController.cpp
TVSpectraViewController.cpp
TVControllerBase.cpp
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
