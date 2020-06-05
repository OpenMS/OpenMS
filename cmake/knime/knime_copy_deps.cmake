INCLUDE(GetPrerequisites)

get_filename_component(BIN_DIR ${DEPS} DIRECTORY)
set(LOOKUP_DIRS "${BIN_DIR};${LOOKUP_DIRS}")

get_prerequisites(${DEPS} DEPENDENCIES 0 1 "${BIN_DIR}" "${LOOKUP_DIRS}" "\$ORIGIN/../lib")

foreach(DEPENDENCY_FILE ${DEPENDENCIES})
  ## Skip things like glibc, gomp, stdcpp
  if (DEPENDENCY_FILE MATCHES "libz" OR DEPENDENCY_FILE MATCHES "libbz" OR NOT DEPENDENCY_FILE MATCHES "^/lib/.*" AND NOT DEPENDENCY_FILE MATCHES "libstdc\\+\\+")
  gp_resolve_item(${DEPS} "${DEPENDENCY_FILE}" "" "${LOOKUP_DIRS}" resolved_file)
  
  get_filename_component(resolved_file_nosymlink ${resolved_file} REALPATH)
  get_filename_component(resolved_file_nosymlink_name ${resolved_file_nosymlink} NAME)
  get_filename_component(resolved_filename ${resolved_file} NAME)
  message("realpath of resolved dependency " ${resolved_file} " is " ${resolved_file_nosymlink})
  file(COPY ${resolved_file_nosymlink} DESTINATION ${TARGET})
  file(RENAME ${TARGET}/${resolved_file_nosymlink_name} ${TARGET}/${resolved_filename})
  else()
  message("Skipped system dependency " ${DEPENDENCY_FILE})
  endif()
endforeach()
