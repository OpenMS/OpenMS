### the directory name
set(directory source/CONCEPT)

### list all filenames of the directory here
set(sources_list
ClassTest.cpp
Constants.cpp
Exception.cpp
FuzzyStringComparator.cpp
GlobalExceptionHandler.cpp
Init.cpp
LogConfigHandler.cpp
LogStream.cpp
PrecisionWrapper.cpp
ProgressLogger.cpp
SingletonRegistry.cpp
StreamHandler.cpp
TypeAsString.cpp
Types.cpp
UniqueIdGenerator.cpp
UniqueIdIndexer.cpp
UniqueIdInterface.cpp
VersionInfo.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\CONCEPT" FILES ${sources})

