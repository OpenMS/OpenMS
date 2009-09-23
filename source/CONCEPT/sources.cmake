### the directory name
set(directory source/CONCEPT)

### list all filenames of the directory here
set(sources_list
ClassTest.C
Constants.C
Exception.C
FuzzyStringComparator.C
ProgressLogger.C
SingletonRegistry.C
VersionInfo.C
LogStream.C
UniqueIdGenerator.C
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

