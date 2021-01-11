### the directory name
set(directory source/VISUAL/APPLICATIONS/GUITOOLS)

### list GUI executables here
set(GUI_executables
INIFileEditor
SwathWizard
TOPPAS
TOPPView
)


### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${GUI_executables})
	list(APPEND sources_VS "${i}.cpp")
endforeach(i)
source_group("" FILES ${sources_VS})
