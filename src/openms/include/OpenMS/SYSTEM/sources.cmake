### the directory name
set(directory include/OpenMS/SYSTEM)

### list all MOC filenames of the directory here
set(sources_list
FileWatcher.h
NetworkGetRequest.h
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### treat as source files, for autoMOC'ing instead of manually calling QT5_WRAP_CPP()
set(OpenMS_sources ${OpenMS_sources} ${sources})
source_group("Source Files\\OpenMS\\SYSTEM" FILES ${sources})

### list all header files of the directory here
set(sources_list_h
ExternalProcess.h
File.h
FileWatcher.h
JavaInfo.h
NetworkGetRequest.h
PythonInfo.h
RWrapper.h
StopWatch.h
SysInfo.h
UpdateCheck.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\SYSTEM" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

