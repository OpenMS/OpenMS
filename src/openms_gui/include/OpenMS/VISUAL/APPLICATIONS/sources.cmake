### the directory name
set(directory include/OpenMS/VISUAL/APPLICATIONS)

### list all header files of the directory here
set(sources_list_h
INIFileEditorWindow.h
SwathWizardBase.h
TOPPViewBase.h
TOPPASBase.h
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
source_group("Header Files\\OpenMS\\VISUAL\\APPLICATIONS" FILES ${sources_h})

