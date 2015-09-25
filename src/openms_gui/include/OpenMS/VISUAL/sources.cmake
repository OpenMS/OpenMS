### the directory name
set(directory include/OpenMS/VISUAL)

### list all MOC filenames of the directory here
set(sources_list
AxisWidget.h
ColorSelector.h
EnhancedTabBar.h
EnhancedWorkspace.h
HistogramWidget.h
ListEditor.h
MetaDataBrowser.h
MultiGradientSelector.h
ParamEditor.h
SpectraIdentificationViewWidget.h
SpectraViewWidget.h
Spectrum1DCanvas.h
Spectrum1DWidget.h
Spectrum2DCanvas.h
Spectrum2DWidget.h
Spectrum3DCanvas.h
Spectrum3DOpenGLCanvas.h
Spectrum3DWidget.h
SpectrumCanvas.h
SpectrumWidget.h
TOPPASEdge.h
TOPPASInputFileListVertex.h
TOPPASLogWindow.h
TOPPASMergerVertex.h
TOPPASOutputFileListVertex.h
TOPPASResource.h
TOPPASResources.h
TOPPASScene.h
TOPPASSplitterVertex.h
TOPPASTabBar.h
TOPPASToolVertex.h
TOPPASTreeView.h
TOPPASVertex.h
TOPPASWidget.h
TOPPViewBehaviorInterface.h
TOPPViewIdentificationViewBehavior.h
TOPPViewSpectraViewBehavior.h
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### Apply MOC compiler
QT4_WRAP_CPP(mocced_sources ${sources} OPTIONS ${BOOST_MOC_ARGS})

### pass source file list to the upper instance
set(OpenMSVisual_sources ${OpenMSVisual_sources} ${mocced_sources})

source_group("Source Files\\OpenMS\\VISUAL" FILES ${mocced_sources})

### list all header files of the directory here
set(sources_list_h
AxisTickCalculator.h
AxisPainter.h
AxisWidget.h
ColorSelector.h
EnhancedTabBar.h
GUIProgressLoggerImpl.h
HistogramWidget.h
LayerData.h
MetaDataBrowser.h
MultiGradient.h
MultiGradientSelector.h
ParamEditor.h
SpectraViewWidget.h
SpectraIdentificationViewWidget.h
Spectrum1DCanvas.h
Spectrum1DWidget.h
Spectrum2DCanvas.h
Spectrum2DWidget.h
Spectrum3DCanvas.h
Spectrum3DOpenGLCanvas.h
Spectrum3DWidget.h
SpectrumCanvas.h
SpectrumWidget.h
ListEditor.h
TOPPASWidget.h
TOPPASScene.h
TOPPASVertex.h
TOPPASToolVertex.h
TOPPASInputFileListVertex.h
TOPPASOutputFileListVertex.h
TOPPASMergerVertex.h
TOPPASEdge.h
TOPPASTabBar.h
TOPPASTreeView.h
TOPPASResource.h
TOPPASResources.h
TOPPViewBehaviorInterface.h
TOPPViewIdentificationViewBehavior.h
TOPPViewSpectraViewBehavior.h
EnhancedWorkspace.h
EnhancedTabBarWidgetInterface.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\VISUAL" FILES ${sources_h})

set(OpenMSVisual_sources_h ${OpenMSVisual_sources_h} ${sources_h})
