# CMake installation file
# Created 2008/10 during OpenMS retreat by CB && AB

PROJECT("OpenMS")

# TODO LIST:
# - install with prefix (chris) => enable this in the test (marc)
# - test for contrib headers (CONTRIB_HEADER_DIR ccmake variable - test one header for each required library)
# - parallel build of test executables

########################################################
### entries meant to be configured using CMake cache ###
### - do NOT hardcode them here!										 ###
### - edit them within CMakeCache.txt using ccmake	 ###
########################################################
# OPENMS_BUILD_MODE
# CONTRIB_CUSTOM_DIR
# OPENMS_DEBUG
# STL_DEBUG
# OPENMS_ENABLE_EXTERNALMEMORY
# QT_DB_PLUGIN

########################################################
###    manual entries (edit this for new release)    ###
########################################################

set(CF_PACKAGE_VERSION "1.3" CACHE INTERNAL "OpenMS VERSION" FORCE)

########################################################
###    end manual																	   ###
########################################################

# Heart of the BUILD system : only edit when you know what you´re doing (we don´t)
# quick manual for most commands: http://www.cmake.org/cmake/help/cmake2.6docs.html
# useful predefined variables: http://www.paraview.org/Wiki/CMake_Useful_Variables

cmake_minimum_required(VERSION 2.6)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if (CMAKE_SIZEOF_VOID_P MATCHES "8")
	set(OPENMS_64BIT_ARCHITECTURE 1 CACHE INTERNAL "Architecture-bits")
else()
	set(OPENMS_64BIT_ARCHITECTURE 0 CACHE INTERNAL "Architecture-bits")
endif()


########################################################
###    compiler flags																 ###
########################################################

if (CMAKE_COMPILER_IS_GNUCXX)
	add_definitions(-Wall -Wextra -pedantic -Wno-non-virtual-dtor -Wno-long-long)
elseif (MSVC)
	#add_definitions(/Wall) ## disable for now.. its just too much!

	## disable stupid warnings. enable if you want to have fun with MSVC 8
	## warning! 4275,4251 is for dllexport/import warnings! this could be vital. recheck!
	#add_definitions(/wd4267 /wd4290 /wd4996 /wd4521 /wd4522 /wd4275 /wd4251)
	#add_definitions(/wd4267)

	##don't warn about unchecked std::copy()
	add_definitions(/D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_DEPRECATE /D_SECURE_SCL=0)

	## disable min and max macros by default
	## see CGAL and BOOST configs for more elaborate explanations
	add_definitions(/DNOMINMAX)

	## xerces bug workaround
	add_definitions(/DOPENMS_XERCESDLL)
	
	# compile such that GSL is in DLL mode
	add_definitions(/DGSL_DLL)

	# minimal code rebuild
	#add_definitions(/Gm /Zi)

	#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT.lib,libcpmt.lib")
	#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT.lib,libcpmt.lib")
else ()
	set(CMAKE_COMPILER_IS_INTELCXX true CACHE INTERNAL "Is Intel C++ compiler (icpc)")
endif()

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)
if (NOT WIN32) # we only want fPIC on non-windows systems (fPIC is implicitly true there)

	CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC)
	if (WITH_FPIC)
		add_definitions(-fPIC)
	endif()
endif()

# -Wconversion flag for GCC
set(CXX_WARN_CONVERSION OFF CACHE BOOL "Enables warnings for type conversion problems (GCC only)")
if (CXX_WARN_CONVERSION)
	if (CMAKE_COMPILER_IS_GNUCXX)
		add_definitions(-Wconversion)
	endif()
endif()
Message(STATUS "Compiler checks for conversion: ${CXX_WARN_CONVERSION}")


#default build type (for Makefile generators): Release
set(OPENMS_BUILD_TYPE "Release" CACHE STRING "Build type: 'Debug' or 'Release'")
if (OPENMS_BUILD_TYPE STREQUAL "Debug")
	Message(STATUS "OpenMS build type: ${OPENMS_BUILD_TYPE} (compiler debug info and OpenMS debug mode enabled)")
	set(CMAKE_BUILD_TYPE "Debug" CACHE INTERNAL "This must be equal to OPENMS_BUILD_TYPE!" FORCE)
elseif (OPENMS_BUILD_TYPE STREQUAL "Release")
	Message(STATUS "OpenMS build type: ${OPENMS_BUILD_TYPE} (optimization enabled)")
	set(CMAKE_BUILD_TYPE "Release" CACHE INTERNAL "This must be equal to OPENMS_BUILD_TYPE!" FORCE)
else()
	Message(FATAL_ERROR "OPENMS_BUILD_TYPE must be 'Debug' or 'Release'!")
endif()

# user OPENMS DEBUG flag
set(OPENMS_ASSERTIONS OFF CACHE BOOL "Enables debug messages (precondition and postconditions are enabled, a bit slower) - this is NOT changing any compiler flags!")
if (OPENMS_BUILD_TYPE STREQUAL "Debug")
	set(OPENMS_ASSERTIONS ON CACHE BOOL "Enables debug messages (precondition and postconditions are enabled, a bit slower) - this is NOT changing any compiler flags!" FORCE)
endif()
Message(STATUS "OpenMS assertion mode: ${OPENMS_ASSERTIONS}")

## STL-DEBUG (only for GCC and in debug mode)
set(STL_DEBUG OFF CACHE BOOL "[GCC only] Enable STL-DEBUG mode (very slow).")
if (STL_DEBUG)
  if (CMAKE_COMPILER_IS_GNUCXX)
		if (OPENMS_BUILD_TYPE STREQUAL "Debug")
			# add compiler flag
    	add_definitions(/D_GLIBCXX_DEBUG)
    	Message(STATUS "STL debug mode: ${STL_DEBUG}")
	  else()
	    Message(WARNING "STL debug mode is supported for OpenMS debug mode only")
	  endif()
  else()
    Message(WARNING "STL debug mode is supported for compiler GCC only")
  endif()
else()
	Message(STATUS "STL debug mode: ${STL_DEBUG}")
endif()

########################################################
###    external libs																 ###
########################################################

## contrib

## disable AndiMS for 64-bit platforms and Visual Studio (until fixed)
if (CMAKE_SIZEOF_VOID_P MATCHES "8" OR MSVC)
	set(USE_ANDIMS 0 CACHE INTERNAL "Support AndiMS?!")
else()
  set(USE_ANDIMS 1 CACHE INTERNAL "Support AndiMS?!")
endif()

set(CONTRIB_CUSTOM_DIR CACHE PATH "User defined location of contrib dir. If left empty we assume the contrib to be in OpenMS/contrib!")
set(CONTRIB_DIR ${PROJECT_SOURCE_DIR}/contrib/ CACHE INTERNAL "Final contrib path after looking at custom_contrib_path. defaults to OpenMS/contrib")

#message("CUSTOM contrib is: " ${CONTRIB_CUSTOM_DIR})

if (CONTRIB_CUSTOM_DIR STREQUAL "")
	message(STATUS "Using DEFAULT setting for contrib directory: ${CONTRIB_DIR}")
  message(STATUS "If you want to set your own path to contrib use -D CONTRIB_CUSTOM_DIR:PATH=\"/my/path/to/contrib\"")
else()
	set(CONTRIB_DIR ${CONTRIB_CUSTOM_DIR})
	message(STATUS "Using user defined contrib directory: ${CONTRIB_DIR}")
endif()

FIND_PATH(CONTRIB_CONFIRM_DIR README_contrib1.3.txt PATHS ${CONTRIB_DIR} NO_DEFAULT_PATH)

if (CONTRIB_CONFIRM_DIR STREQUAL "CONTRIB_CONFIRM_DIR-NOTFOUND")
	message(FATAL_ERROR "The path to the Contrib directory (${CONTRIB_DIR}) is not correct or you need to update your contrib!")
endif()

set(CONTRIB_INCLUDE_DIR ${CONTRIB_DIR}/include/ CACHE INTERNAL "contrib include dir")
set(CONTRIB_LIB_DIR ${CONTRIB_DIR}/lib/ CACHE INTERNAL "contrib lib dir")

message(STATUS "Contrib lib dir: ${CONTRIB_LIB_DIR}")

###########################################################
###							 find libs (for linking)								###
###########################################################
## on windows we need the *.lib versions (dlls alone won't do for linking)
## Also we need to be very careful to link against the correct debug/release version in MSVC!
## TODO: check if we also need to distinguish between compiler versions

## fills ${varname} with the names of Debug and Release libraries (which usually only differ on MSVC)
## @param varname Name of the variable which will hold the result string (e.g. "optimized myLib.so debug myLibDebug.so")
## @param libnames   List of library names which are searched (release libs)
## @param libnames_d List of library names which are searched (debug libs)
## @param human_libname Name of the library (for display only) 
MACRO (OPENMS_CHECKLIB varname libnames libnames_d human_libname)
	FIND_LIBRARY(${varname}_OPT NAMES ${libnames} PATHS ${CONTRIB_LIB_DIR} DOC "${human_libname} library dir" NO_DEFAULT_PATH)
	if (${varname}_OPT STREQUAL "${varname}_OPT-NOTFOUND")
		MESSAGE(FATAL_ERROR "Unable to find ${human_libname} library! Searched names are: [${libnames}] Please make sure it is part of the contrib (which we assume to be at: ${CONTRIB_DIR}")
	else()
		MESSAGE(STATUS "Found ${human_libname} library at: " ${${varname}_OPT})
	endif()
	FIND_LIBRARY(${varname}_DBG NAMES ${libnames_d} PATHS ${CONTRIB_LIB_DIR} DOC "${human_libname} (Debug) library dir" NO_DEFAULT_PATH)
	if (${varname}_DBG STREQUAL "${varname}_DBG-NOTFOUND")
		MESSAGE(FATAL_ERROR "Unable to find ${human_libname} (Debug) library! Searched names are: [${libnames}] Please make sure it is part of the contrib (which we assume to be at: ${CONTRIB_DIR}")
	else()
		MESSAGE(STATUS "Found ${human_libname} library at: " ${${varname}_DBG})
	endif()
	## combine result and include "optimized" and "debug" keywords which are essential for target_link_libraries()
	set(${varname} optimized ${${varname}_OPT} debug ${${varname}_DBG})
ENDMACRO (OPENMS_CHECKLIB)

OPENMS_CHECKLIB(CONTRIB_XERCESC "xerces-c_3;xerces-c_static_3;libxerces-c;xerces-c" "xerces-c_3D;xerces-c_static_3D;libxerces-c;xerces-c" "xerces_c")

OPENMS_CHECKLIB(CONTRIB_GSL "libgsl_dll;gsl" "libgsl_dll_d;gsl" "GSL")

OPENMS_CHECKLIB(CONTRIB_GSLCBLAS "libgslcblas_dll;gslcblas" "libgslcblas_dll_d;gslcblas" "GSL-CBLAS")

OPENMS_CHECKLIB(CONTRIB_CGAL "cgal-vc80-mt;cgal-vc90-mt;CGAL" "cgal-vc80-mt-gd;cgal-vc90-mt-gd;CGAL" "CGAL")

# TODO windoes????
if (NOT MSVC AND USE_ANDIMS)
	OPENMS_CHECKLIB(CONTRIB_ANDIMS "ms10lib-TODO_CHRIS_JA_DU;ms10lib.a" "ms10lib-TODO_CHRIS_JA_DU;ms10lib.a" "ms10lib")

	OPENMS_CHECKLIB(CONTRIB_NETCDF_CPP "netCDF-TODO_CHRIS_JA_DU;libnetcdf_c++.a" "netCDF-TODO_CHRIS_JA_DU;libnetcdf_c++.a" "netCDF_c++")
	OPENMS_CHECKLIB(CONTRIB_NETCDF_C "netCDF-TODO_CHRIS_JA_DU;libnetcdf.a" "netCDF-TODO_CHRIS_JA_DU;libnetcdf.a" "netCDF")
	set(CONTRIB_NETCDF ${CONTRIB_NETCDF_CPP} ${CONTRIB_NETCDF_C})
endif()

# libSVM:

if (MSVC) # we created a lib here, because we cannot just add an object..
	OPENMS_CHECKLIB(SVM_LIBRARY "svm" "svmd" "LibSVM")
	#set(SVM_LIBRARY optimized ${CONTRIB_LIB_DIR}svm.obj debug ${CONTRIB_LIB_DIR}svm_d.obj CACHE INTERNAL "object file of svm library")
else()
	# there is no FIND_OBJECT.. so just hard-code it (FIND_LIBRARY won't work)
	set(SVM_LIBRARY general ${CONTRIB_LIB_DIR}svm.o CACHE INTERNAL "object file of svm library") # do not distinguish between optimized and debug
endif()


INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/include/) ## for configured files, i.e. config.h
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include/)
INCLUDE_DIRECTORIES(${CONTRIB_INCLUDE_DIR})

## QT
SET(QT_USE_QTXML 1)
SET(QT_USE_QTNETWORK 1)
SET(QT_USE_QTSQL 1)
SET(QT_USE_QTOPENGL 1)
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})
INCLUDE(UseQt4)

Message(STATUS "QT qmake at ${QT_QMAKE_EXECUTABLE}")
Message(STATUS "QT moc at ${QT_MOC_EXECUTABLE}")
Message(STATUS "QT uic at ${QT_UIC_EXECUTABLE}")

########################################################
###    useful programms															 ###
########################################################

find_program(SVNVERSION_EXECUTABLE
             svnversion
             PATHS "c:/programme/subversion/bin" "c:/program files/subversion/bin"   ### additional search paths (along with $PATH)
             DOC "svnversion executable which helps in determining the svn revision when building TOPP tools")

if (SVNVERSION_EXECUTABLE STREQUAL "SVNVERSION_EXECUTABLE-NOTFOUND")
	message(STATUS "The programm svnversion could not be found. SVN-revision information will not be available! Add the location of svnversion(.exe) to your PATH environment variable if you require SVN-revision.")
	set(OPENMS_HAS_SVNVERSION OFF CACHE INTERNAL "SVNVersion(.exe) present?")
else()
	set(OPENMS_HAS_SVNVERSION ON CACHE INTERNAL "SVNVersion(.exe) present?")
endif()

########################################################
###    configure config.h														 ###
########################################################

## define some directories
if (DO_INSTALL)
	set(CF_OPENMS_DATA_PATH ${CMAKE_INSTALL_PREFIX}/share/OpenMS-${CF_PACKAGE_VERSION} CACHE INTERNAL "Path to the shared documents of OpenMS.")
else()
	set(CF_OPENMS_DATA_PATH ${PROJECT_SOURCE_DIR}/share/OpenMS CACHE INTERNAL "Path to the shared documents of OpenMS.")
endif()

set(CF_OPENMS_TEST_DATA_PATH ${PROJECT_SOURCE_DIR}/source/TEST/data/ CACHE INTERNAL "Path to the test data")

# check for Microsoft Visual Studio compiler
if (MSVC)
	set(OPENMS_COMPILER_MSVC "1" CACHE INTERNAL "Do we use Microsoft Compiler?")
endif()
# check for G++
if (CMAKE_COMPILER_IS_GNUCXX)
	set(OPENMS_COMPILER_GXX "1" CACHE INTERNAL "Do we use G++ Compiler?")
endif()

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(OPENMS_BIG_ENDIAN)

# check 32/64 bit architecture (defined above!)
if (NOT DEFINED OPENMS_64BIT_ARCHITECTURE)
	message(FATAL_ERROR "Cmake script was re-ordered and is now invalid! Please make sure that OPENMS_64BIT_ARCHITECTURE is defined when config.h.in is configured!")
endif()


include(CheckTypeSize) #: Check sizeof a type
CHECK_TYPE_SIZE("unsigned char" SIZE_UCHAR)
CHECK_TYPE_SIZE("unsigned short" SIZE_USHORT)
CHECK_TYPE_SIZE("unsigned int" SIZE_UINT)
CHECK_TYPE_SIZE("unsigned long" SIZE_ULONG)
CHECK_TYPE_SIZE("unsigned long long" SIZE_ULONGLONG)
CHECK_TYPE_SIZE("short" SIZE_SHORT)
CHECK_TYPE_SIZE("int" SIZE_INT)
CHECK_TYPE_SIZE("long" SIZE_LONG)
CHECK_TYPE_SIZE("long long" SIZE_LONGLONG)

CHECK_TYPE_SIZE("int32_t" SIZE_INT32)
if (HAVE_SIZE_INT32)
	set(CF_OPENMS_INT32_TYPE int32_t)
else()
	# search for another Int32 type
	if (SIZE_INT MATCHES "4")
		set(CF_OPENMS_INT32_TYPE int)
	elseif (SIZE_SHORT MATCHES "4")
		set(CF_OPENMS_INT32_TYPE short)
	elseif (SIZE_LONG MATCHES "4")
		set(CF_OPENMS_INT32_TYPE long)
	else()
		message(FATAL_ERROR "Cannot find signed 32bit integer type. Please contact the developers!")
	endif()
endif()

CHECK_TYPE_SIZE("int64_t" SIZE_INT64)
if (HAVE_SIZE_INT64)
	set(CF_OPENMS_INT64_TYPE int64_t)
else()
	# search for another Int64 type
	if (SIZE_INT MATCHES "8")
		set(CF_OPENMS_INT64_TYPE int)
	elseif (SIZE_LONG MATCHES "8")
		set(CF_OPENMS_INT64_TYPE long)
	elseif (SIZE_LONGLONG MATCHES "8")
		set(CF_OPENMS_INT64_TYPE "long long")
	else()
		message(FATAL_ERROR "Cannot find signed 64bit integer type. Please contact the developers!")
	endif()
endif()


#### this is currently disabled due to too many problems
#### e.g. calling setValue("myParam",114) will not work, because it setValue would require a setValue(String, int) signature,
####      but only has setValue(String, long long). Casting everything to (Int) (which is long long on 64bit VS) is unintuitiv!)
## set OpenMS Types (Int) and (UInt) to the size of size_t and not simply as "int" or "unsigned int"
## - under unix this seems equivalent, but MSVS defined "int" always as 32bit, even on 64bit platforms.

#if (OPENMS_64BIT_ARCHITECTURE)
#	set(CF_OPENMS_SIZE_T_SIGNED ${CF_OPENMS_INT64_TYPE})
#else()
#	set(CF_OPENMS_SIZE_T_SIGNED ${CF_OPENMS_INT32_TYPE})
#endif()



CHECK_TYPE_SIZE("uint8_t" SIZE_UINT8)
if (HAVE_SIZE_UINT8)
	set(CF_OPENMS_BYTE_TYPE uint8_t)
else()
	# search for another uint8 type
	if (SIZE_UCHAR MATCHES "1")
		set(CF_OPENMS_BYTE_TYPE "unsigned char")
	elseif (SIZE_USHORT MATCHES "1")
		set(CF_OPENMS_BYTE_TYPE "unsigned short")
	else()
		message(FATAL_ERROR "Cannot find unsigned 8bit integer (byte) type. Please contact the developers!")
	endif()
endif()


CHECK_TYPE_SIZE("uint64_t" SIZE_UINT64)
if (HAVE_SIZE_UINT64)
	set(CF_OPENMS_UINT64_TYPE uint64_t)
else()
	# search for another uint64 type
	if (SIZE_ULONG MATCHES "8")
		set(CF_OPENMS_UINT64_TYPE "unsigned long")
	elseif (SIZE_ULONGLONG MATCHES "8")
		set(CF_OPENMS_UINT64_TYPE "unsigned long long")
	else()
		message(FATAL_ERROR "Cannot find uint64 type. Please contact the developers!")
	endif()
endif()

## system headers:
include(CheckIncludeFileCXX) # Check if the include file exists.

CHECK_INCLUDE_FILE_CXX("unistd.h" OPENMS_HAS_UNISTD_H)
CHECK_INCLUDE_FILE_CXX("process.h" OPENMS_HAS_PROCESS_H)

CHECK_INCLUDE_FILE_CXX("time.h" OPENMS_HAS_TIME_H)
CHECK_INCLUDE_FILE_CXX("sys/types.h" OPENMS_HAS_SYS_TYPES_H)
CHECK_INCLUDE_FILE_CXX("sys/times.h" OPENMS_HAS_SYS_TIMES_H)
CHECK_INCLUDE_FILE_CXX("sys/time.h"  OPENMS_HAS_SYS_TIME_H)

# user flag
set(OPENMS_ENABLE_EXTERNALMEMORY OFF CACHE BOOL "User switch to enable external memory support (if you run out of RAM+swap). Warning: performance penalty!")

include(CheckFunctionExists)
## in MinGW we have the signal.h header, but no kill() as in Linux, so we need to check for the kill() function
CHECK_FUNCTION_EXISTS("kill" OPENMS_HAS_KILL)
CHECK_FUNCTION_EXISTS("sysconf" OPENMS_HAS_SYSCONF)

## ANDIMS flag USE_ANDIMS defined above!
if (NOT DEFINED USE_ANDIMS)
	message(FATAL_ERROR "Cmake script was re-ordered and is now invalid! Please make sure that USE_ANDIMS is defined when config.h.in is configured!")
endif()

# user flag with default "QMYSQL"
set(QT_DB_PLUGIN "QMYSQL" CACHE STRING "User switch to change the Qt database plugin.")


########################################################
###    configure Doxyfile														 ###
########################################################

## CF_PACKAGE_VERSION (already defined above in manual section)

set(CF_OPENMS_PATH ${PROJECT_BINARY_DIR})

## replace any variables in Doxyfile.in with current values
configure_file(${PROJECT_SOURCE_DIR}/source/config/Doxyfile.in ${PROJECT_BINARY_DIR}/doc/doxygen/Doxyfile)

### TODO: test this -> when building documentation, set in environment using svnversion before calling doxygen
#set(CF_PACKAGE_REVISION_FILE)
#if (OPENMS_HAS_SVNVERSION)
#
#	set(svnversion_file_h ${PROJECT_SOURCE_DIR}/source/CONCEPT/openms_svnversion) # file with version number
#
#	## ensure that the file is always rebuild
#  add_custom_target(svnversion_file_target ALL ${SVNVERSION_EXECUTABLE} -n ${PROJECT_SOURCE_DIR} > ${svnversion_file_h} COMMENT "Call svnversion in project directory to generate openms_svnversion.txt")
#	set(CF_PACKAGE_REVISION_FILE ${svnversion_file_h})
#	#file(READ ${svnversion_file_h} ENV_OPENMS_REVISION) # read revision from file
#	#set($ENV{OPENMS_REVISION} ${ENV_OPENMS_REVISION}) # set environment variable e.g. for doxygen
#endif()

# are we building a shared or static lib?! (BOTH within the same BUILD-tree is NOT possible with OpenMS!!)
if (MINGW OR MSVC)
	set(BUILD_SHARED_LIBS false)
else()
	set(BUILD_SHARED_LIBS true)
endif()


## replace any variables in config.h.in with current values
set (CONFIGURED_CONFIG_H ${PROJECT_BINARY_DIR}/include/OpenMS/config.h)
configure_file(${PROJECT_SOURCE_DIR}/source/config/config.h.in ${CONFIGURED_CONFIG_H})


########################################################
###    BUILD the lib																 ###
########################################################

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

MACRO (QT4_WRAP_UI_OWN outfiles )
  QT4_EXTRACT_OPTIONS(ui_files ui_options ${ARGN})

  # create output directory (will not exist for out-of-source builds)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${directory})

  FOREACH (it ${ui_files})
    GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
    GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
    SET(outfile ${PROJECT_BINARY_DIR}/${directory}/ui_${outfile}.h)
    ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
      COMMAND ${QT_UIC_EXECUTABLE}
      ARGS ${ui_options} -o ${outfile} ${infile}
      MAIN_DEPENDENCY ${infile})
    SET(${outfiles} ${${outfiles}} ${outfile})
  ENDFOREACH (it)

ENDMACRO (QT4_WRAP_UI_OWN)

# big include file for headers and C files, which fills the OpenMS_sources variable
include (includes.cmake)

#needed for MSVC to locate Boost-Libs (put this above ADD_LIBRARY() - otherwise it will not work)
link_directories(${CONTRIB_LIB_DIR})

# add library target
# warning: set BUILD_SHARED_LIBS to decide if library is shared or static (see above)! We need the BUILD_SHARED_LIBS flag to set declspec flags for MSVC!
ADD_LIBRARY(OpenMS ${OpenMS_sources})


# always create the svn version file before compiling OpenMS
if (OPENMS_HAS_SVNVERSION)
## TODO replace this with PreProcessorMacro when compiling lib
## e.g. SET_SOURCE_FILE_PROPERTIES(libSource.c main.c COMPILE_FLAGS -DREVISON=234)
	message(STATUS "Adding dependency on SVN-Version target file...")
##	add_dependencies(OpenMS svnversion_file_target)
endif()


set(OPENMS_DEP_LIBRARIES ${QT_LIBRARIES} ${CONTRIB_CGAL} ${CONTRIB_GSL} ${CONTRIB_GSLCBLAS} ${SVM_LIBRARY} ${CONTRIB_XERCESC} ${CONTRIB_ANDIMS} ${CONTRIB_NETCDF})
if (MSVC)
	list(APPEND OPENMS_DEP_LIBRARIES opengl32.lib)
	target_link_libraries(OpenMS ${OPENMS_DEP_LIBRARIES})
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  target_link_libraries(OpenMS ${OPENMS_DEP_LIBRARIES})
endif()

### OPENMS_LIBRARIES defines the libraries used by OpenMS; this should be used to link against executables
set(OPENMS_LIBRARIES OpenMS ${OPENMS_DEP_LIBRARIES})

####### TOPP #########
set(TOPP_executables)
include(source/APPLICATIONS/TOPP/executables.cmake)



set(TOPP_DIR source/APPLICATIONS/TOPP/)
foreach(i ${TOPP_executables})
	set(resource_file ${CMAKE_SOURCE_DIR}/${TOPP_DIR}${i}.rc)
        set(resource_dir ${CMAKE_SOURCE_DIR}/${TOPP_DIR}/${i}-resources/)
	## add icons to TOPPView and INIFileEditor
	if (MSVC AND EXISTS ${resource_file})
		message(STATUS "Setting resource file ${resource_file} for ${i}")
		add_executable(${i} ${TOPP_DIR}${i}.C ${resource_file})
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND EXISTS ${resource_dir}) 
    # transform TOPPView and INIFileEditor to real Mac OS X applications
    add_executable(${i} MACOSX_BUNDLE ${TOPP_DIR}${i}.C)          

    # Need to copy the icon file and the custom Info.plist
    # if version info changes, please update also the 
    # info.plist file located in APPLICATIONS/TOPP/{i}-resources/
    if(${CMAKE_GENERATOR} MATCHES "Xcode") 
      # since Xcode builds several configurations (Debug;Release)
      # we need a special handling for the postbuild actions
      ADD_CUSTOM_COMMAND(
        TARGET ${i}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E echo "Post build action for ${i}.app -- Creating OS X bundle"
        COMMAND ${CMAKE_COMMAND} -E copy
        ${resource_dir}Info.plist
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/\$\(CONFIGURATION\)/${i}.app/Contents
        COMMAND ${CMAKE_COMMAND} -E make_directory
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/\$\(CONFIGURATION\)/${i}.app/Contents/Resources
        COMMAND ${CMAKE_COMMAND} -E copy
        ${resource_dir}${i}.icns
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/\$\(CONFIGURATION\)/${i}.app/Contents/Resources
        )
    else()
      ADD_CUSTOM_COMMAND(
        TARGET ${i}
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E echo "Post build action for ${i}.app -- Creating OS X bundle"
        COMMAND ${CMAKE_COMMAND} -E copy
        ${resource_dir}Info.plist
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i}.app/Contents
        COMMAND ${CMAKE_COMMAND} -E make_directory
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i}.app/Contents/Resources
        COMMAND ${CMAKE_COMMAND} -E copy
        ${resource_dir}${i}.icns
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i}.app/Contents/Resources
        )
    endif()
  else()
		add_executable(${i} ${TOPP_DIR}${i}.C)
	endif()
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(TOPP)
add_dependencies(TOPP ${TOPP_executables})

##### UTILS ########
set(UTILS_executables)
include(source/APPLICATIONS/UTILS/executables.cmake)

foreach(i ${UTILS_executables})
	add_executable(${i} source/APPLICATIONS/UTILS/${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(UTILS)
add_dependencies(UTILS ${UTILS_executables})


##### Doc progs ######
set(OpenMS_doc_executables)
include(doc/doxygen/parameters/executables.cmake)
foreach(i ${OpenMS_doc_executables})
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "")
	add_executable(${i} EXCLUDE_FROM_ALL ${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
endforeach(i)
add_custom_target(doc_progs)
add_dependencies(doc_progs ${OpenMS_doc_executables})

###### Examples ######
set(EXAMPLES_executables)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY source/EXAMPLES)
include(source/EXAMPLES/executables.cmake)
foreach(i ${EXAMPLES_executables})
  add_executable(${i} source/EXAMPLES/${i}.C)
  target_link_libraries(${i} ${OPENMS_LIBRARIES})
endforeach(i)
add_custom_target(EXAMPLES)
add_dependencies(EXAMPLES ${EXAMPLES_executables})


#### Tests #####
enable_testing()
INCLUDE(Dart)
add_subdirectory(source/TEST)

##### TOPPtests #####
add_subdirectory(source/TEST/TOPP)

#### configs for external include files, e.g. to develop code with OpenMS
#get_target_property(OPENMS_LIBRARIES_EXTERNAL OpenMS LOCATION)
#set(OPENMS_INCLUDE_EXTERNAL ${PROJECT_BINARY_DIR}/include)

configure_file(
"${PROJECT_SOURCE_DIR}/OpenMSConfig.cmake.in"
"${PROJECT_SOURCE_DIR}/OpenMSConfig.cmake"
@ONLY
)

include(CMakeExportBuildSettings)
cmake_export_build_settings(${OpenMS_BINARY_DIR}/OpenMSBuildSettings.cmake)
export_library_dependencies(${OpenMS_BINARY_DIR}/OpenMSLibraryDepends.cmake)
# TODO install these files to a specific location to have a standard way to be used with external code

##### targets list #####
add_custom_target(targets
	                COMMAND ${CMAKE_COMMAND} -E echo ""
	                COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
	                COMMAND ${CMAKE_COMMAND} -E echo ""
                  COMMAND ${CMAKE_COMMAND} -E echo "The following make targets are available:"
                  COMMAND ${CMAKE_COMMAND} -E echo "    [no target]    builds the OpenMS library, TOPP tools and UTILS tools"
                  COMMAND ${CMAKE_COMMAND} -E echo "    OpenMS         builds the OpenMS library"
                  COMMAND ${CMAKE_COMMAND} -E echo "    TOPP           builds the TOPP tools"
                  COMMAND ${CMAKE_COMMAND} -E echo "    UTILS          builds the UTILS tools"
                  COMMAND ${CMAKE_COMMAND} -E echo "    TEST           builds the OpenMS tests"
                  COMMAND ${CMAKE_COMMAND} -E echo "    test           executes OpenMS and TOPP tests"
                  COMMAND ${CMAKE_COMMAND} -E echo "                   make sure they are built using the TEST target"
                  COMMAND ${CMAKE_COMMAND} -E echo "    EXAMPLES       builds the example executables in source/EXAMPLES"
                  COMMAND ${CMAKE_COMMAND} -E echo "    help           list all available targets (very long)"
                  COMMAND ${CMAKE_COMMAND} -E echo ""
                  COMMAND ${CMAKE_COMMAND} -E echo "Single TOPP tools have their own target, e.g. TOPPView"
                  COMMAND ${CMAKE_COMMAND} -E echo "The same applies for class tests and UTILS tools."
	                COMMAND ${CMAKE_COMMAND} -E echo ""
	                COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
	                COMMAND ${CMAKE_COMMAND} -E echo ""
                  COMMENT "The most important targets for OpenMS"
                  VERBATIM
  )

##### Message after OpenMS has been built #####
if (MSVC)
  add_custom_command(TARGET OpenMS
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "The OpenMS library has been built."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "You should now build the TOPP tools and tests."
                    COMMAND ${CMAKE_COMMAND} -E echo "Then you should test your installation by executing the tests."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Execute the 'targets' project to see prominent targets!"
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMENT "message after library is built"
                    VERBATIM)
else()
  add_custom_command(TARGET OpenMS
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "The OpenMS library has been built."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "You should now build the TOPP tools and tests."
                    COMMAND ${CMAKE_COMMAND} -E echo "Then you should test your installation by executing the tests."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Make sure to add the OpenMS and contrib lib/ path"
                    COMMAND ${CMAKE_COMMAND} -E echo "to your LD_LIBRARY_PATH environment variable."
                    COMMAND ${CMAKE_COMMAND} -E echo "Otherwise the tests and TOPP tools will not work."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "For a full list of targets execute:"
                    COMMAND ${CMAKE_COMMAND} -E echo "make targets"
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMENT "message after library is built"
                    VERBATIM)
endif()

##### Message after TOPP has been built #####
if (MSVC)
  add_custom_command(TARGET TOPP
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "The TOPP tools have been built and installed to the bin/ folder."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Then you should test your installation by executing the tests."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Execute the 'targets' project to see prominent targets!"
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMENT "message after TOPP is built"
                    VERBATIM)
else()
  add_custom_command(TARGET TOPP
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "The TOPP tools have been built and installed to the bin/ folder."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Then you should test your installation by executing the tests."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "Make sure to add the OpenMS and contrib lib/ path"
                    COMMAND ${CMAKE_COMMAND} -E echo "to your LD_LIBRARY_PATH environment variable."
                    COMMAND ${CMAKE_COMMAND} -E echo "Otherwise the tests and TOPP tools will not work."
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "For a full list of make targets execute:"
                    COMMAND ${CMAKE_COMMAND} -E echo "make targets"
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
                    COMMAND ${CMAKE_COMMAND} -E echo ""
                    COMMENT "message after TOPP is built"
                    VERBATIM)
endif()


##### Messages at the end of cmake #####
MESSAGE(STATUS "")
MESSAGE(STATUS "-----------------------------------------------------------------")
MESSAGE(STATUS "")
MESSAGE(STATUS "You have successfully configured OpenMS and TOPP.")
MESSAGE(STATUS "")
if (MSVC)
  MESSAGE(STATUS "Execute the 'targets' project to see prominent targets!")
else()
  MESSAGE(STATUS "For a full list of make targets execute:")
  MESSAGE(STATUS "'make targets'")
endif()
MESSAGE(STATUS "")
MESSAGE(STATUS "-----------------------------------------------------------------")
MESSAGE(STATUS "")

##### install/copy?? #####
if(NOT MSVC)
  # TODO:
  #  - this code breaks the windows cmake call
  #  - install compiled documentation ??
  #  - install share
  set(DO_INSTALL false CACHE BOOL "Set to 'ON' if you want the CMAKE_INSTALL_PREFIX to be used and the library to be installed. On 'OFF' the library and Tools will be build but not installed!")

  ## CPack installation and packaging procedures
  install(TARGETS OpenMS LIBRARY DESTINATION lib)

  # install utils
  foreach(util ${UTILS_executables})
    install(TARGETS ${util}
      RUNTIME DESTINATION bin
      BUNDLE DESTINATION bin
      COMPONENT applications)
  endforeach()

  # install TOPP Tools
  foreach(TOPP_exe ${TOPP_executables})
    INSTALL(TARGETS ${TOPP_exe} 
      RUNTIME DESTINATION bin
      BUNDLE DESTINATION bin
      COMPONENT applications)
  endforeach()

  # install headers
  foreach(OpenMS_header ${OpenMS_sources_h})
    #extract path
    string(REGEX MATCH
      "^include/(.*)/"
      path_name ${OpenMS_header})

    # get path below include/..
    set(path_name ${CMAKE_MATCH_1})

    INSTALL(FILES
      ${OpenMS_header}
      DESTINATION include/${path_name}
      COMPONENT headers)
  endforeach()

  # since config h is not contained in OpenMS_sources_h we need to address it separtly 
  install(FILES ${CONFIGURED_CONFIG_H} DESTINATION include/OpenMS)

  # install mocced sources into include
  foreach(mocced ${mocced_sources})
    #extract path
    string(REGEX MATCH
      "^/.*/include/(.*)/"
      path_name ${mocced})

    # get path below include/..
    set(path_name ${CMAKE_MATCH_1})

    INSTALL(FILES
      ${mocced} 
      DESTINATION include/${path_name}
      COMPONENT headers)
  endforeach()

  # install ui files
  foreach(uiced ${uiced_sources})
    #extract path
    string(REGEX MATCH
      "^/.*/include/(.*)/"
      path_name ${uiced})
    set(path_name ${CMAKE_MATCH_1})
    INSTALL(FILES
      ${uiced} 
      DESTINATION include/${path_name}
      COMPONENT headers)
  endforeach()

  # install share
  INSTALL(DIRECTORY share/
    DESTINATION share
    COMPONENT share)

  # most of the interesting CPack stuff is done below
  include(OpenMSCPack.cmake)
endif()
