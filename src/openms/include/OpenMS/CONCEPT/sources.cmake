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
GlobalExceptionHandler.h
LogStream.h
LogConfigHandler.h
Macros.h
ProgressLogger.h
SingletonRegistry.h
StreamHandler.h
Types.h
UniqueIdGenerator.h
UniqueIdInterface.h
UniqueIdIndexer.h
VersionInfo.h
BinaryComposeFunctionAdapter.h
UnaryComposeFunctionAdapter.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\CONCEPT" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

