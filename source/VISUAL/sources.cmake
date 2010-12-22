### the directory name
set(directory source/VISUAL)

### list all filenames of the directory here
set(sources_list
Annotation1DDistanceItem.C
Annotation1DItem.C
Annotation1DPeakItem.C
Annotation1DTextItem.C
Annotations1DContainer.C
AxisTickCalculator.C
AxisWidget.C
ColorSelector.C
EnhancedTabBar.C
HistogramWidget.C
LayerData.C
MetaDataBrowser.C
MultiGradient.C
MultiGradientSelector.C
ParamEditor.C
SpectraViewWidget.C
SpectraIdentificationViewWidget.C
Spectrum1DCanvas.C
Spectrum1DWidget.C
Spectrum2DCanvas.C
Spectrum2DWidget.C
Spectrum3DCanvas.C
Spectrum3DOpenGLCanvas.C
Spectrum3DWidget.C
SpectrumCanvas.C
SpectrumWidget.C
ListEditor.C
TOPPASWidget.C
TOPPASScene.C
TOPPASVertex.C
TOPPASToolVertex.C
TOPPASInputFileListVertex.C
TOPPASOutputFileListVertex.C
TOPPASMergerVertex.C
TOPPASEdge.C
TOPPASTabBar.C
TOPPASTreeView.C
TOPPASResource.C
TOPPASResources.C
TOPPViewBehaviorInterface.C
TOPPViewIdentificationViewBehavior.C
TOPPViewSpectraViewBehavior.C
EnhancedWorkspace.C
EnhancedTabBarWidgetInterface.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\VISUAL" FILES ${sources})


### icons
# add   : icons are added to source/VISUAL/ICONS/resources.qrc
# remove: after removing an icon, you have to rerun 'cmake' to fix the dependencies
QT4_ADD_RESOURCES(qt_resource_file source/VISUAL/ICONS/resources.qrc)
set(OpenMS_sources ${OpenMS_sources} ${qt_resource_file})
