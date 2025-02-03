set(OpenMS_IO_sources  CACHE INTERNAL "This variable should hold all OpenMS_IO sources at the end of the config step" )
set(OpenMS_IO_sources_h  CACHE INTERNAL "This variable should hold all OpenMS_IO headers at the end of the config step" )

## ATTENTION: The order of includes should be similar to the inclusion hierarchy
include(source/FORMAT/DATAACCESS/sources.cmake)
include(source/FORMAT/HANDLERS/sources.cmake)
include(source/FORMAT/MSNUMPRESS/sources.cmake)
include(source/FORMAT/VALIDATORS/sources.cmake)
include(source/FORMAT/OPTIONS/sources.cmake)
include(source/FORMAT/sources.cmake)

## add configured config.h&Co to source group
source_group("Header Files\\OpenMS_IO" FILES ${OpenMS_IO_configured_headers})
## merge all headers to sources (for source group view in VS)
list(APPEND OpenMS_IO_sources ${OpenMS_IO_sources_h} ${OpenMS_IO_configured_headers})

# Remove duplicates
list(REMOVE_DUPLICATES OpenMS_IO_sources)