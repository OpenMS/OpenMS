### the directory name
set(directory include/OpenMS/SYSTEM)

### list all MOC filenames of the directory here
set(sources_list
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
BuildInfo.h
CurlInit.h
ExternalProcess.h
File.h
JavaInfo.h
Network.h
NetworkGetRequest.h
PathUtils.h
PythonInfo.h
RWrapper.h
StopWatch.h
SysInfo.h
SystemSettings.h
TempFiles.h
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

### Private (non-installed) header: SIMDe.h pulls in <simde/x86/ssse3.h> and, on
### MSVC, defines operators on simde__m128i. Its own comment says to include it
### from .cpp files only, and nothing but libOpenMS sources does. Keeping it off
### OpenMS_sources_h is what lets SIMDe be a PRIVATE link dependency: no SIMDe
### type or include appears in any installed header.
set(private_headers_list_h
SIMDe.h
)
set(private_sources_h)
foreach(i ${private_headers_list_h})
	list(APPEND private_sources_h ${directory}/${i})
endforeach(i)
source_group("Header Files\\OpenMS\\SYSTEM" FILES ${private_sources_h})
set(OpenMS_private_headers ${OpenMS_private_headers} ${private_sources_h})
