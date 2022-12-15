set(OpenMSML_sources  CACHE INTERNAL "This variable should hold all OpenMS sources at the end of the config step" )

## added to OpenMSML library: ${OpenMSML_sources}
include(source/ML/sources.cmake)

set(OpenMSML_sources_h  CACHE INTERNAL "This variable should hold all OpenMS sources at the end of the config step" )

## added to OpenMSML library: ${OpenMSML_sources_h}
include(include/OpenMS/ML/sources.cmake)

## merge all headers to sources (for source group view in VS)
list(APPEND OpenMSML_sources ${OpenMSML_sources_h})
