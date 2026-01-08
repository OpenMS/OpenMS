### the directory name
set(directory include/OpenMS/CONCEPT)

### list all header files of the directory here
set(sources_list_h
ClassTest.h
Colorizer.h
CommonEnums.h
Constants.h
EnumHelpers.h
Exception.h
FuzzyStringComparator.h
GlobalExceptionHandler.h
Helpers.h
Init.h
LogConfigHandler.h
LogStream.h
Macros.h
MacrosTest.h
PrecisionWrapper.h
ProgressLogger.h
Qt5Port.h
RAIICleanup.h
StreamHandler.h
Types.h
UniqueIdGenerator.h
UniqueIdIndexer.h
UniqueIdInterface.h
VersionInfo.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\CONCEPT" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})
