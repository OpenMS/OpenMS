### the directory name
set(directory include/OpenMS/VISUAL/UIC)

### list all filenames of the directory here
set(sources_list
ParamEditor.ui
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### Apply UIC compiler
QT4_WRAP_UI_OWN(uiced_sources ${sources})

### pass source file list to the upper instance
set(OpenMSVisual_sources ${OpenMSVisual_sources} ${uiced_sources})

