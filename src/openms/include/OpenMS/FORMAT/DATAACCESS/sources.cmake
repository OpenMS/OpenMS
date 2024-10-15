### the directory name
set(directory include/OpenMS/FORMAT/DATAACCESS)

### list all header files of the directory here
set(sources_list_h
  MSDataAggregatingConsumer.h
  MSDataCachedConsumer.h
  MSDataChainingConsumer.h
  MSDataStoringConsumer.h
  MSDataSqlConsumer.h
  MSDataTransformingConsumer.h
  MSDataWritingConsumer.h
  NoopMSDataConsumer.h
  SiriusFragmentAnnotation.h
  SwathFileConsumer.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\FORMAT\\DATAACCESS" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

