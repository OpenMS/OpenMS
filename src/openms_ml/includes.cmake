set(OpenMSML_sources  CACHE INTERNAL "This variable should hold all OpenMS sources at the end of the config step" )

## added to OpenMSVisual library: ${OpenMSML_sources}
# TODO include(source/VISUAL/ANNOTATION/sources.cmake)

set(OpenMSML_sources_h  CACHE INTERNAL "This variable should hold all OpenMS sources at the end of the config step" )

## added to OpenMSVisual library: ${OpenMSML_sources}
# TODO include(include/OpenMS/VISUAL/ANNOTATION/sources.cmake)

## merge all headers to sources (for source group view in VS)
list(APPEND OpenMSML_sources ${OpenMSML_sources_h})
