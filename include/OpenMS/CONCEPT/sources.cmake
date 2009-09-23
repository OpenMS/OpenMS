### the directory name
set(directory include/OpenMS/CONCEPT)

### list all header files of the directory here
set(sources_list_h
ClassTest.h
Constants.h
Exception.h
Factory.h
FactoryBase.h
FuzzyStringComparator.h
Macros.h
ProgressLogger.h
SingletonRegistry.h
Types.h
VersionInfo.h
LogStream.h
UniqueIdGenerator.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\CONCEPT" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

