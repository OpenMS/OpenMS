#------------------------------------------------------------------------------
## fills ${varname} with the names of Debug and Release libraries (which usually only differ on MSVC)
## @param varname Name of the variable which will hold the result string (e.g. "optimized myLib.so debug myLibDebug.so")
## @param libnames   List of library names which are searched (release libs)
## @param libnames_d List of library names which are searched (debug libs)
## @param human_libname Name of the library (for display only)
macro (OPENMS_CHECKLIB varname libnames libnames_d human_libname)
	# force find_library to run again
	set(${varname}_OPT "${varname}_OPT-NOTFOUND" CACHE FILEPATH "Cleared." FORCE)
	find_library(${varname}_OPT NAMES ${libnames} PATHS ${CONTRIB_LIB_DIR} DOC "${human_libname} library dir" NO_DEFAULT_PATH)
	if ("${varname}_OPT" STREQUAL "${varname}_OPT-NOTFOUND")
		message(FATAL_ERROR "Unable to find ${human_libname} library! Searched names are: [${libnames}]\nPlease make sure it is part of the contrib (which we assume to be in either of these directories: ${CONTRIB_LIB_DIR}). Set custom contrib paths using the CMAKE_FIND_ROOT_PATH variable in CMake.")
	else()
		message(STATUS "Found ${human_libname} library (Release) at: " ${${varname}_OPT})
	endif()
	# force find_library to run again
	set(${varname}_DBG "${varname}_DBG-NOTFOUND" CACHE FILEPATH "Cleared." FORCE)
	find_library(${varname}_DBG NAMES ${libnames_d} PATHS ${CONTRIB_LIB_DIR} DOC "${human_libname} (Debug) library dir" NO_DEFAULT_PATH)
	if ("${varname}_DBG" STREQUAL "${varname}_DBG-NOTFOUND")
		message(FATAL_ERROR "Unable to find ${human_libname} (Debug) library! Searched names are: [${libnames}]\nPlease make sure it is part of the contrib (which we assume to be in either of these directories: ${CONTRIB_LIB_DIR}). Set custom contrib paths using the CMAKE_FIND_ROOT_PATH variable in CMake.")
	else()
		message(STATUS "Found ${human_libname} library (Debug) at: " ${${varname}_DBG})
	endif()
	## combine result and include "optimized" and "debug" keywords which are essential for target_link_libraries()
	set(${varname} optimized ${${varname}_OPT} debug ${${varname}_DBG})
endmacro (OPENMS_CHECKLIB)

#------------------------------------------------------------------------------
## Copy the dll produced by the given target to the test/doc binary path.
## @param targetname The target to modify.
## @note This macro will do nothing with non MSVC generators.
macro(copy_dll_to_extern_bin targetname)
  if(MSVC)
    get_target_property(WIN32_DLLLOCATION ${targetname} LOCATION)
    get_filename_component(WIN32_DLLPATH ${WIN32_DLLLOCATION} PATH)

    ## copy OpenMS.dll to test executables dir "$(TargetFileName)" is a placeholder filled by VS at runtime
    file(TO_NATIVE_PATH "${WIN32_DLLPATH}/$(TargetFileName)" DLL_SOURCE)

    file(TO_NATIVE_PATH "${OPENMS_HOST_BINARY_DIRECTORY}/src/tests/class_tests/bin/$(ConfigurationName)/$(TargetFileName)" DLL_TEST_TARGET)
    file(TO_NATIVE_PATH "${OPENMS_HOST_BINARY_DIRECTORY}/src/tests/class_tests/bin/$(ConfigurationName)" DLL_TEST_TARGET_PATH)

    file(TO_NATIVE_PATH "${OPENMS_HOST_BINARY_DIRECTORY}/doc/doxygen/parameters/$(ConfigurationName)/$(TargetFileName)" DLL_DOC_TARGET)
    file(TO_NATIVE_PATH "${OPENMS_HOST_BINARY_DIRECTORY}/doc/doxygen/parameters/$(ConfigurationName)" DLL_DOC_TARGET_PATH)


    add_custom_command(TARGET ${targetname}
                      POST_BUILD
                      COMMAND ${CMAKE_COMMAND} -E make_directory "${DLL_TEST_TARGET_PATH}"
                      COMMAND ${CMAKE_COMMAND} -E copy ${DLL_SOURCE} ${DLL_TEST_TARGET}
                      COMMAND ${CMAKE_COMMAND} -E make_directory "${DLL_DOC_TARGET_PATH}"
                      COMMAND ${CMAKE_COMMAND} -E copy ${DLL_SOURCE} ${DLL_DOC_TARGET})
  endif(MSVC)
endmacro()

#------------------------------------------------------------------------------
## export a single option indicating if boost static libs should be preferred
option(BOOST_USE_STATIC "Use Boost static libraries." ON)

#------------------------------------------------------------------------------
## Wraps the common find boost code into a single call
## @param .. simply add all required components to the call
## @note This macro will define BOOST_MOC_ARGS that should be added to all moc
##       calls (see https://bugreports.qt-project.org/browse/QTBUG-22829)
macro(find_boost)
  set(Boost_USE_STATIC_LIBS ${BOOST_USE_STATIC})
  set(Boost_USE_MULTITHREADED  ON)
  set(Boost_USE_STATIC_RUNTIME OFF)
  add_definitions(/DBOOST_ALL_NO_LIB) ## disable auto-linking of boost libs (boost tends to guess wrong lib names)
  set(Boost_COMPILER "")

  # help boost finding it's packages
  set(Boost_ADDITIONAL_VERSIONS "1.47.0" "1.48.0" "1.49.0" "1.50.0" "1.51.0" "1.52.0" "1.53.0" "1.54.0")

  # 1st attempt does not explicitly requires boost to enable second check (see below)
  find_package(Boost 1.42.0 COMPONENTS ${ARGN})

  set(BOOST_MOC_ARGS "")

  # see: https://bugreports.qt-project.org/browse/QTBUG-22829
  # Confirmed only on mac os x and leads to problems on win32 and lnx
  # so we handle it for now only on mac os x and boost versions > 1.52
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" OR ${Boost_MINOR_VERSION} GREATER "52")
  	set(BOOST_MOC_ARGS "-DBOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION")
  endif()
endmacro(find_boost)

#------------------------------------------------------------------------------
## Unity Build of a set of cpp files
## i.e., make one large compilation unit which usually compiles a lot faster
##
## pros: 
##    - compiles a lot faster 
##    - reveals double definitions of internal classes in different cpp files
##    - reveals weird variable names (like 'IN'), which are defined as macro in certain headers which are now part of the compilation unit
## cons:
##    - not desirable for developing, since a small change will trigger recreation of the whole lib
function(convert_to_unity_build UB_SUFFIX SOURCE_FILES_NAME)
   set(files ${${SOURCE_FILES_NAME}})
   # generate a unique file name for the unity build translation unit
   set(unit_build_file ${CMAKE_CURRENT_BINARY_DIR}/ub_${UB_SUFFIX}.cpp)
   # exclude all translation units from compilation
   set_source_files_properties(${files} PROPERTIES HEADER_FILE_ONLY true)
   # open the unity build file
   FILE(WRITE ${unit_build_file} "// Unity Build generated by CMake\n")
   # add include statement for each translation unit
   foreach(source_file ${files} )
     # we have headers in there as well, which should not be included explicitly
     if (${source_file} MATCHES "\\.cpp|\\.cxx") # cxx for moc's; 
       if (IS_ABSOLUTE ${source_file}) 
         FILE( APPEND ${unit_build_file} "#include<${source_file}>\n")
       else() 
         FILE( APPEND ${unit_build_file} "#include<${CMAKE_CURRENT_SOURCE_DIR}/${source_file}>\n")
       endif() 
     endif()
   endforeach(source_file)
   # add unity build aggregate as source file 
   set(${SOURCE_FILES_NAME} ${${SOURCE_FILES_NAME}} ${unit_build_file} PARENT_SCOPE)
endfunction(convert_to_unity_build)
#------------------------------------------------------------------------------
## Checks if the user supplied package type is valid and aborts if not
## @param package_type The given package type
macro(is_valid_package package_type)
  list(FIND VALID_PACKAGE_TYPES ${package_type} list_pos)
  if( ${list_pos} EQUAL -1 )
  	message(STATUS "The PACKAGE_TYPE ${package_type} is invalid")
  	message(STATUS "Valid PACKAGE_TYPEs are:")
  	foreach( _vpt ${VALID_PACKAGE_TYPES} )
  		message(STATUS " * ${_vpt}")
  	endforeach()
  	message(FATAL_ERROR "Aborting ...")
  endif()
endmacro()

#------------------------------------------------------------------------------
## Wrap add_library() and also deal with unity build at the same time
macro(oms_add_library libname sources_variable_name)
	if (ENABLE_UNITYBUILD)
		MESSAGE(STATUS "Unity Build for ${libname} lib: Enabled (ENABLE_UNITYBUILD=ON)")
		convert_to_unity_build(${libname} ${sources_variable_name})
	else()
		MESSAGE(STATUS "Unity Build for ${libname} lib: Disabled (ENABLE_UNITYBUILD=OFF)")
	endif()
	add_library(${libname} ${${sources_variable_name}})
endmacro()
