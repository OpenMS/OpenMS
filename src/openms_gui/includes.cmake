set(OpenMSVisual_sources  CACHE INTERNAL "This variable should hold all OpenMS sources at the end of the config step" )

## added to OpenMSVisual library: ${OpenMSVisual_sources}
include(source/VISUAL/ANNOTATION/sources.cmake)
include(source/VISUAL/APPLICATIONS/sources.cmake)
include(source/VISUAL/APPLICATIONS/MISC/sources.cmake)
include(source/VISUAL/DIALOGS/sources.cmake)
include(source/VISUAL/INTERFACES/sources.cmake)
include(source/VISUAL/MISC/sources.cmake)
include(source/VISUAL/VISUALIZER/sources.cmake)
include(source/VISUAL/VISITORS/sources.cmake)
include(source/VISUAL/sources.cmake)
#include(include/OpenMS/VISUAL/UIC/sources.cmake) ## uic are "sources" of OpenMS because they add .ui depedencies to the lib
#include(include/OpenMS/VISUAL/DIALOGS/UIC/sources.cmake) ## uic are "sources" of OpenMS because they add .ui depedencies to the lib

set(OpenMSVisual_sources_h  CACHE INTERNAL "This variable should hold all OpenMS sources at the end of the config step" )

## added to OpenMSVisual library: ${OpenMSVisual_sources}
include(include/OpenMS/VISUAL/ANNOTATION/sources.cmake)
include(include/OpenMS/VISUAL/APPLICATIONS/sources.cmake)
include(include/OpenMS/VISUAL/INTERFACES/sources.cmake)
include(include/OpenMS/VISUAL/MISC/sources.cmake)
include(include/OpenMS/VISUAL/VISITORS/sources.cmake)
include(include/OpenMS/VISUAL/APPLICATIONS/MISC/sources.cmake)  ## MOC sources are included here
include(include/OpenMS/VISUAL/sources.cmake)					 					## and here ...
include(include/OpenMS/VISUAL/DIALOGS/sources.cmake)   					## and here ...
include(include/OpenMS/VISUAL/VISUALIZER/sources.cmake)					## and here ...

## merge all headers to sources (for source group view in VS)
list(APPEND OpenMSVisual_sources ${OpenMSVisual_sources_h})

# TODO track why the duplicate warnings are thrown for all (!) MOC sources
# Macro problem?
#list(REMOVE_DUPLICATES OpenMSVisual_sources)

