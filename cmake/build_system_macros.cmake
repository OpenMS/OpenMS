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
