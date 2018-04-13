INCLUDE(GetPrerequisites)
get_prerequisites(${DEPS} DEPENDENCIES 1 1 "" "${LOOKUP_DIRS}")

get_filename_component(BIN_DIR ${DEPS} DIRECTORY)
set(LOOKUP_DIRS "${BIN_DIR};${LOOKUP_DIRS}")

foreach(DEPENDENCY_FILE ${DEPENDENCIES})
  gp_resolve_item(${DEPS} "${DEPENDENCY_FILE}" "" "${LOOKUP_DIRS}" resolved_file)
  
  get_filename_component(resolved_file_nosymlink ${resolved_file} REALPATH)
  #message("realpath of resolved dependency " ${resolved_file} " is " ${resolved_file_nosymlink})
  file(COPY ${resolved_file_nosymlink} DESTINATION ${TARGET})
endforeach()
