### the directory name
set(directory include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS)

### list all header files of the directory here
set(sources_list_h
DataAccessHelper.h
MRMFeatureAccessOpenMS.h
SimpleOpenMSSpectraAccessFactory.h
SpectrumAccessOpenMS.h
SpectrumAccessOpenMSCached.h
SpectrumAccessOpenMSInMemory.h
SpectrumAccessSqMass.h
SpectrumAccessTransforming.h
SpectrumAccessQuadMZTransforming.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\OPENSWATH" FILES ${sources_h})
set_source_files_properties(${directory}/sources.cmake PROPERTIES HEADER_FILE_ONLY TRUE)

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

