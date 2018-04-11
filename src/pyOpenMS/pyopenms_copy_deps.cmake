INCLUDE(GetPrerequisites)
get_prerequisites(${DEPS} DEPENDENCIES 1 1 "" "${LOOKUP_DIRS}")
#message ("Found deps for ${DEPS}: "  ${DEPENDENCIES})

get_filename_component(BIN_DIR ${DEPS} DIRECTORY)
set(LOOKUP_DIRS "${BIN_DIR};${LOOKUP_DIRS}")

#message ("Searching in: " ${LOOKUP_DIRS})
foreach(DEPENDENCY_FILE ${DEPENDENCIES})
  gp_resolve_item(${DEPS} "${DEPENDENCY_FILE}" "" "${LOOKUP_DIRS}" resolved_file)
  #message(resolved_file='${resolved_file}')
  file(COPY ${resolved_file} DESTINATION ${TARGET})
endforeach()