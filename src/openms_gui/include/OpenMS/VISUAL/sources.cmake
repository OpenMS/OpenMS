### the directory name
set(directory include/OpenMS/VISUAL)

### list all header files of the directory here
set(sources_list_h
AxisTickCalculator.h
AxisPainter.h
AxisWidget.h
ColorSelector.h
DataSelectionTabs.h
DIATreeTab.h
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
Plot1DCanvas.h
Plot1DWidget.h
Plot2DCanvas.h
Plot2DWidget.h
Plot3DCanvas.h
Plot3DOpenGLCanvas.h
Plot3DWidget.h
PlotCanvas.h
PlotWidget.h
RecentFilesMenu.h
SpectraTreeTab.h
SpectraIDViewTab.h
SwathLibraryStats.h
TableView.h
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
TOPPViewMenu.h
TreeView.h
TVControllerBase.h
TVDIATreeTabController.h
TVIdentificationViewController.h
TVSpectraViewController.h
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

