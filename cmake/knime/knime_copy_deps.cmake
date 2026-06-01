INCLUDE(GetPrerequisites)

get_filename_component(BIN_DIR ${DEPS} DIRECTORY)
set(LOOKUP_DIRS "${BIN_DIR};${LOOKUP_DIRS}")

get_prerequisites(${DEPS} DEPENDENCIES 0 1 "${BIN_DIR}" "${LOOKUP_DIRS}" "\$ORIGIN/../lib")

foreach(DEPENDENCY_FILE ${DEPENDENCIES})
  ## Skip things like glibc
  if (DEPENDENCY_FILE MATCHES ".*/libc\\..*"
   OR DEPENDENCY_FILE MATCHES ".*/libm\\..*"
   OR DEPENDENCY_FILE MATCHES ".*/libdl\\..*"
   OR DEPENDENCY_FILE MATCHES ".*/libpthread\\..*"
   OR DEPENDENCY_FILE MATCHES ".*/ld-linux-.*"
   OR DEPENDENCY_FILE MATCHES ".*/linux-vdso.*"
  )
    message("Skipped system dependency " ${DEPENDENCY_FILE})
  else()
    gp_resolve_item(${DEPS} "${DEPENDENCY_FILE}" "" "${LOOKUP_DIRS}" resolved_file)
    get_filename_component(resolved_file_nosymlink ${resolved_file} REALPATH)
    get_filename_component(resolved_file_nosymlink_name ${resolved_file_nosymlink} NAME)
    get_filename_component(resolved_filename ${resolved_file} NAME)
    message("realpath of resolved dependency " ${resolved_file} " is " ${resolved_file_nosymlink})
    file(COPY ${resolved_file_nosymlink} DESTINATION ${TARGET})
    file(RENAME ${TARGET}/${resolved_file_nosymlink_name} ${TARGET}/${resolved_filename})
  endif()
endforeach()
