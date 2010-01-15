# Locate Intel Threading Building Blocks include paths and libraries
# TBB can be found at http://www.threadingbuildingblocks.org/ 
# Originally Written by Hannes Hofmann, hannes.hofmann _at_ informatik.uni-erlangen.de
# Modified by Andreas Bertsch

# This module requires
# MT_TBB_INCLUDE_DIR
# TBB_LIB_DIR

# This module defines
# TBB_INCLUDE_DIRS, where to find task_scheduler_init.h, etc.
# TBB_LIBRARY_DIRS, where to find libtbb, libtbbmalloc
# TBB_INSTALL_DIR, the base TBB install directory
# TBB_LIBRARIES, the libraries to link against to use TBB.
# TBB_DEBUG_LIBRARIES, the libraries to link against to use TBB with debug symbols.
# TBB_FOUND, If false, don't try to use TBB.

if (WIN32)
  set(_TBB_LIB_NAME "tbb")
  set(_TBB_LIB_MALLOC_NAME "${_TBB_LIB_NAME}malloc")
  set(_TBB_LIB_DEBUG_NAME "${_TBB_LIB_NAME}_debug")
  set(_TBB_LIB_MALLOC_DEBUG_NAME "${_TBB_LIB_MALLOC_NAME}_debug")
endif (WIN32)

if (UNIX)
    if (APPLE)
        # MAC
        set(_TBB_DEFAULT_INSTALL_DIR "/Library/Frameworks/Intel_TBB.framework/Versions")
        # libs: libtbb.dylib, libtbbmalloc.dylib, *_debug
        set(_TBB_LIB_NAME "tbb")
        set(_TBB_LIB_MALLOC_NAME "${_TBB_LIB_NAME}malloc")
        set(_TBB_LIB_DEBUG_NAME "${_TBB_LIB_NAME}_debug")
        set(_TBB_LIB_MALLOC_DEBUG_NAME "${_TBB_LIB_MALLOC_NAME}_debug")
        # has only one flavor: ia32/cc4.0.1_os10.4.9
        set(_TBB_COMPILER "cc4.0.1_os10.4.9")
        set(_TBB_ARCHITECTURE "ia32")
    else (APPLE)
	# LINUX
        set(_TBB_LIB_NAME "tbb")
        set(_TBB_LIB_MALLOC_NAME "${_TBB_LIB_NAME}malloc")
        set(_TBB_LIB_DEBUG_NAME "${_TBB_LIB_NAME}_debug")
        set(_TBB_LIB_MALLOC_DEBUG_NAME "${_TBB_LIB_MALLOC_NAME}_debug")
    endif (APPLE)
endif (UNIX)

#-- Clear the public variables
set (TBB_FOUND "NO")

#-- Look for libraries
find_library(TBB_LIBRARY        ${_TBB_LIB_NAME}        ${MT_TBB_LIBRARY_DIR} NO_DEFAULT_PATH)
find_library(TBB_MALLOC_LIBRARY ${_TBB_LIB_MALLOC_NAME} ${MT_TBB_LIBRARY_DIR} NO_DEFAULT_PATH)
mark_as_advanced(TBB_LIBRARY TBB_MALLOC_LIBRARY)

find_library(TBB_LIBRARY_DEBUG        ${_TBB_LIB_DEBUG_NAME}        ${MT_TBB_LIBRARY_DIR} NO_DEFAULT_PATH)
find_library(TBB_MALLOC_LIBRARY_DEBUG ${_TBB_LIB_MALLOC_DEBUG_NAME} ${MT_TBB_LIBRARY_DIR} NO_DEFAULT_PATH)
mark_as_advanced(TBB_LIBRARY_DEBUG TBB_MALLOC_LIBRARY_DEBUG)

if (MT_TBB_INCLUDE_DIR)

	# sample tbb source code to test includes
	FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testTBBInclude.C
	"#include <tbb/parallel_for.h>
	int main() 
	{
		return 0;	
	}")

	try_compile(TBB_COMPILE_SUCCESS_HEADERS 
							${CMAKE_BINARY_DIR} 
							${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testTBBInclude.C
							CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MT_TBB_INCLUDE_DIR}")
	if (NOT TBB_COMPILE_SUCCESS_HEADERS)
		message (STATUS "Could NOT find TBB headers.")
	else()
		message (STATUS "Intel Threading Building Blocks includes found in "${MT_TBB_INCLUDE_DIR})
		if (MT_TBB_LIBRARY_DIR)
			# sample tbb source code to test library linking
			FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testTBBLinking.C
		  "#include <tbb/task_scheduler_init.h>
			using namespace tbb;
		  int main()
		  {
				task_scheduler_init init;
		    return 0;
		  }")
			#set(LIST TBB_LINK_LIBS ${TBB_LIBRARY} ${TBB_MALLOC_LIBRARY})
			try_compile(TBB_COMPILE_SUCCESS_LIBRARIES
              ${CMAKE_BINARY_DIR}
              ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testTBBLinking.C
              CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MT_TBB_INCLUDE_DIR}" 
							CMAKE_FLAGS "-DLINK_DIRECTORIES=${MT_TBB_LIBRARY_DIR}"
							CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${_TBB_LIB_NAME};${_TBB_LIB_MALLOC_NAME}")
			
			if (TBB_COMPILE_SUCCESS_LIBRARIES)
				set (TBB_LIBRARIES ${TBB_LIBRARY} ${TBB_MALLOC_LIBRARY})
				message (STATUS "Intel Threading Building Block libraries found in "${MT_TBB_LIBRARY_DIR})
				set (TBB_FOUND "YES")
			endif()
		endif()
	endif()
endif()

if (NOT TBB_FOUND)
    message("ERROR: Intel TBB NOT found!")
    message(STATUS "Looked for Threading Building Blocks includes in ${MT_TBB_INCLUDE_DIR}")
		message(STATUS "Looked for Threading Building Blocks libraries in ${MT_TBB_LIBRARY_DIR}")

    # do only throw fatal, if this pkg is REQUIRED
    if (TBB_FIND_REQUIRED)
        message(FATAL_ERROR "Could NOT find TBB library.")
    endif (TBB_FIND_REQUIRED)
endif (NOT TBB_FOUND)

