# CMake installation file
# Created 2008/10 during OpenMS retreat by CB && AB

PROJECT("OpenMS")

# TODO LIST:
# - test for contrib headers (CONTRIB_HEADER_DIR ccmake variable - test one header for each required library)

########################################################
### entries meant to be configured using CMake cache ###
### - do NOT hardcode them here!                     ###
### - edit them within CMakeCache.txt using ccmake	 ###
########################################################
# CONTRIB_CUSTOM_DIR
# OPENMS_BUILD_TYPE
# CONTRIB_CUSTOM_DIR
# STL_DEBUG
# DB_TEST
# QT_DB_PLUGIN
# MT_CUDA_BUILD_TYPE

########################################################
###    manual entries (edit this for new release)    ###
########################################################

set(OPENMS_PACKAGE_VERSION_MAJOR "1")
set(OPENMS_PACKAGE_VERSION_MINOR "9")
set(OPENMS_PACKAGE_VERSION_PATCH "0")

########################################################
###            end manual entries                    ###
########################################################

## Heart of the BUILD system : only edit when you know what you´re doing (we don´t)
## quick manual for most commands: http://www.cmake.org/cmake/help
## useful predefined variables: http://www.paraview.org/Wiki/CMake_Useful_Variables
##
## naming conventions:
##
## prefix a variable with 'CF_' if it is used to configure a file!
## e.g., CF_LibOpenMSExport



cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
set(CF_OPENMS_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION_MAJOR}.${OPENMS_PACKAGE_VERSION_MINOR}.${OPENMS_PACKAGE_VERSION_PATCH}" CACHE INTERNAL "OpenMS VERSION" FORCE)

## sanity check: sometimes CMAKE_SIZEOF_VOID_P just vanishes (when updating CMake).
if (NOT CMAKE_SIZEOF_VOID_P)
	message(FATAL_ERROR "'CMAKE_SIZEOF_VOID_P' is undefined. Thus you should delete CMakeFiles (the directory) and the CMakeCache.txt and rerun CMake again! This is some weird CMake bug that seems to appear when updating the CMake version.")
endif()

if (CMAKE_SIZEOF_VOID_P MATCHES "8")
	set(OPENMS_64BIT_ARCHITECTURE 1 CACHE INTERNAL "Architecture-bits")
else()
	set(OPENMS_64BIT_ARCHITECTURE 0 CACHE INTERNAL "Architecture-bits")
endif()

# some of our own macros (OPENMS_CHECKLIB, QT4_WRAP_UI_OWN, ...)
include (cmake/macros.cmake)

########################################################
###    compiler flags																 ###
########################################################

## Fill this with compile flags that external projects should use as well
## for OpenMS internal flags (not promoted to external compiler flags) append to CMAKE_CXX_FLAGS
set(CF_OPENMS_ADDCXX_FLAGS) ## see OpenMS/cmake/OpenMSConfig.cmake.in to see how its configured and used (i.e. as OPENMS_ADDCXX_FLAGS)

if (CMAKE_COMPILER_IS_GNUCXX)

	# determine gcc version
	execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GNUCXX_VERSION)
	string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GNUCXX_VERSION})
	list(GET GCC_VERSION_COMPONENTS 0 GNUCXX_MAJOR_VERSION)
	list(GET GCC_VERSION_COMPONENTS 1 GNUCXX_MINOR_VERSION)

  add_definitions(-Wall -Wextra -Wno-non-virtual-dtor -Wno-long-long -Wno-variadic-macros)
  if (NOT MT_ENABLE_CUDA)  # necessary since CUDA contains non-pedantic code
		add_definitions(--pedantic)
	endif()

	## Recommended setting for eclipse, see http://www.cmake.org/Wiki/CMake:Eclipse
	if (CMAKE_GENERATOR STREQUAL "Eclipse CDT4 - Unix Makefiles")
		add_definitions(-fmessage-length=0)
	endif()

	if (NOT OPENMS_64BIT_ARCHITECTURE AND ${GNUCXX_MAJOR_VERSION} MATCHES "4" AND ${GNUCXX_MINOR_VERSION} MATCHES "3")
		add_definitions(-march=i486)
	endif()

elseif (MSVC)
	## do not use add_definitions
	## add definitions also lands in stuff like RC_DEFINITION which tend to fail if you use
	## Eclipse CDT 4 - NMAKE generator
	## use set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} ...") instead

	## disable dll-interface warning
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4251 /wd4275")

	## disable deprecated functions warning (e.g. for POSIX functions)
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4996")

	## disable explicit template instantiation request for partially defined classes
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4661")

	## disable warning: "decorated name length exceeded, name was truncated"
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4503")
	
	## don't warn about unchecked std::copy()
	add_definitions(/D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_DEPRECATE)

	## xerces bug workaround
	add_definitions(/DOPENMS_XERCESDLL)

	## removed: (GSL is static now) compile such that GSL is in DLL mode
	#add_definitions(/DGSL_DLL)

	## FeatureFinder.obj is huge and won't compile in VS2008 debug otherwise:
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")

	## use multiple Cores (if available)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")

  if (NOT OPENMS_64BIT_ARCHITECTURE)
    ## enable SSE1 on 32bit, on 64bit the compiler flag does not exist
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE")
  endif()
elseif ("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_COMPILER_IS_CLANG true CACHE INTERNAL "Is CLang compiler (clang++)")
else()
	set(CMAKE_COMPILER_IS_INTELCXX true CACHE INTERNAL "Is Intel C++ compiler (icpc)")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CF_OPENMS_ADDCXX_FLAGS}")


#########################################################
#                SOME MAC OS X STUFF                    #
#########################################################

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	message(STATUS "Setting CMAKE_OSX_DEPLOYMENT_TARGET to 10.5")
	set(CMAKE_OSX_DEPLOYMENT_TARGET 10.5 CACHE INTERNAL "Set OS X Deployment Target to 10.5 to build packages compatible to leopard" FORCE )
  # force 10.5 SDK
  set (CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.5.sdk")
	# due to some awkward updates of cmake in the 2.8.x series we need to add the following
  find_library(CARBON NAMES Carbon)
  find_package(OpenGL)
  set(APPLE_EXTRA_LIBS ${CARBON} ${OPENGL_LIBRARY})

  # force cmake > 2.8.1 to avoid bugs with OSX_DEPLOYMENT_TARGET
  cmake_minimum_required(VERSION 2.8.1 FATAL_ERROR)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

##### Parallelization things
#################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

## cuda checking
set(MT_ENABLE_CUDA OFF CACHE BOOL "Should CUDA support be enabled (version <= 2.1 currently supported).")
set(MT_CUDA_BUILD_TYPE Device CACHE STRING "Switch between Device and Emulation mode.")
if (MT_ENABLE_CUDA)
	MESSAGE(STATUS "NVIDIA CUDA: ${MT_CUDA_BUILD_TYPE} mode")
	find_package(Cuda)
else()
  MESSAGE(STATUS "NVIDIA CUDA: ${MT_ENABLE_CUDA}")
endif()


if (FOUND_CUDART AND MT_ENABLE_CUDA) #CUDART and not CUDA!!!
	add_definitions(/DOPENMS_HAS_CUDA)
	MESSAGE (STATUS "Found CUDA header files in: ${FOUND_CUDA_NVCC_INCLUDE}") #Not CUDA_INCLUDE_DIR
	MESSAGE (STATUS "Found CUDA library at: ${FOUND_CUDART}")
	INCLUDE_DIRECTORIES(${FOUND_CUDA_NVCC_INCLUDE})
  set (CUDA_NVCC_INCLUDE_ARGS ${CUDA_NVCC_INCLUDE_ARGS} -I ${PROJECT_BINARY_DIR}/include/)
endif()

## TBB
set(MT_ENABLE_TBB OFF CACHE BOOL "Should Intel Threading Building Blocks support be enabled.")
set(MT_TBB_INCLUDE_DIR CACHE PATH "Intel Threading Building Blocks 'include' directory.")
set(MT_TBB_LIBRARY_DIR CACHE PATH "Intel Threading Building Blocks libraries directory.")


MESSAGE(STATUS "Intel TBB: ${MT_ENABLE_TBB}")
if (MT_ENABLE_TBB)
	find_package(TBB)
endif()

if (TBB_FOUND)
	INCLUDE_DIRECTORIES(${MT_TBB_INCLUDE_DIR})
	add_definitions(/DOPENMS_HAS_TBB)
endif()

## OpenMP
set(MT_ENABLE_OPENMP ON CACHE BOOL "Should OpenMP support be enabled")
if (MT_ENABLE_OPENMP)
	find_package(OpenMP)
endif()

if (OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}") ## do NOT use add_definitions() here, because RC.exe on windows will fail
	if (NOT MSVC)
		set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
endif()

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)
if (NOT WIN32) # we only want fPIC on non-windows systems (fPIC is implicitly true there)
	CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC)
	if (WITH_FPIC)
		add_definitions(-fPIC)
	endif()
endif()

## -Wconversion flag for GCC
set(CXX_WARN_CONVERSION OFF CACHE BOOL "Enables warnings for type conversion problems (GCC only)")
if (CXX_WARN_CONVERSION)
	if (CMAKE_COMPILER_IS_GNUCXX)
		add_definitions(-Wconversion)
	endif()
endif()
Message(STATUS "Compiler checks for conversion: ${CXX_WARN_CONVERSION}")


## default build type (for Makefile generators): Release
set(OPENMS_BUILD_TYPE "Release" CACHE STRING "Build type: 'Debug' or 'Release'")
if ("${OPENMS_BUILD_TYPE}" STREQUAL "Debug")
	Message(STATUS "OpenMS build type: ${OPENMS_BUILD_TYPE} (compiler debug info and OpenMS debug mode enabled)")
	set(CMAKE_BUILD_TYPE "Debug" CACHE INTERNAL "This must be equal to OPENMS_BUILD_TYPE!" FORCE)
	set(OPENMS_ASSERTIONS 1 CACHE INTERNAL "Enables debug messages (precondition and postconditions are enabled, a bit slower) - this is NOT changing any compiler flags!" FORCE)
elseif ("${OPENMS_BUILD_TYPE}" STREQUAL "Release")
	Message(STATUS "OpenMS build type: ${OPENMS_BUILD_TYPE} (optimization enabled)")
	set(CMAKE_BUILD_TYPE "Release" CACHE INTERNAL "This must be equal to OPENMS_BUILD_TYPE!" FORCE)
	set(OPENMS_ASSERTIONS 0 CACHE INTERNAL "Enables debug messages (precondition and postconditions are enabled, a bit slower) - this is NOT changing any compiler flags!" FORCE)
else()
	Message(FATAL_ERROR "OPENMS_BUILD_TYPE must be 'Debug' or 'Release'!")
endif()

## STL-DEBUG (only for GCC and in debug mode)
set(STL_DEBUG OFF CACHE BOOL "[GCC only] Enable STL-DEBUG mode (very slow).")
if (STL_DEBUG)
  if (CMAKE_COMPILER_IS_GNUCXX)
		if ("${OPENMS_BUILD_TYPE}" STREQUAL "Debug")
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
	Message(STATUS "[GCC only] STL debug mode: ${STL_DEBUG}")
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

set(CONTRIB_CUSTOM_DIR CACHE DOC "User defined location of contrib dir. If left empty we assume the contrib to be in OpenMS/contrib!")
set(CONTRIB_DIR ${PROJECT_SOURCE_DIR}/contrib/ CACHE INTERNAL "Final contrib path after looking at custom_contrib_path. defaults to OpenMS/contrib")

#message("CUSTOM contrib is: " ${CONTRIB_CUSTOM_DIR})

if ("${CONTRIB_CUSTOM_DIR}" STREQUAL "")
	message(STATUS "Using DEFAULT setting for contrib directory: ${CONTRIB_DIR}")
  message(STATUS "If you want to set your own path to contrib use:")
	message(STATUS "   cmake <..more options..> -D CONTRIB_CUSTOM_DIR:PATH=\"/my/path/to/contrib\" \"${PROJECT_SOURCE_DIR}\"")
else()
  get_filename_component(CONTRIB_CUSTOM_DIR_PATH ${CONTRIB_CUSTOM_DIR} ABSOLUTE)
	set(CONTRIB_DIR ${CONTRIB_CUSTOM_DIR_PATH})
	message(STATUS "Using user defined contrib directory: ${CONTRIB_CUSTOM_DIR_PATH}")
endif()

# maybe we should change this method to a more robust one
FIND_PATH(CONTRIB_CONFIRM_DIR README_contrib1.7.txt PATHS ${CONTRIB_DIR} NO_DEFAULT_PATH)

if ("${CONTRIB_CONFIRM_DIR}" STREQUAL "CONTRIB_CONFIRM_DIR-NOTFOUND")
	message(FATAL_ERROR "Error finding the contrib! Either\n  - The path to the Contrib binary(!) directory is not correct (we looked at '${CONTRIB_DIR}')\n  - You need to compile the contrib first!\n  - You need to update your contrib!\nNote: If you are using a relative path to your Contrib binary directory, it must be relative to the OpenMS source directory!")
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

#set which library extensions are preferred (we want static libraries)
if(NOT MSVC)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
endif()
if (APPLE)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib;.a")
endif()

OPENMS_CHECKLIB(CONTRIB_XERCESC "xerces-c_3;xerces-c_static_3;libxerces-c;xerces-c" "xerces-c_3D;xerces-c_static_3D;libxerces-c;xerces-c" "xerces_c")

OPENMS_CHECKLIB(CONTRIB_GSL "gsl" "gsl_d;gsl" "GSL")

OPENMS_CHECKLIB(CONTRIB_GSLCBLAS "cblas;gslcblas" "cblas_d;gslcblas" "GSL-CBLAS")

OPENMS_CHECKLIB(CONTRIB_CBC1 "libCbc;Cbc" "libCbcd;Cbc" "COIN-OR Cbc")
OPENMS_CHECKLIB(CONTRIB_CBC2 "libCgl;Cgl" "libCgld;Cgl" "COIN-OR Cgl")
OPENMS_CHECKLIB(CONTRIB_CBC3 "libClp;Clp" "libClpd;Clp" "COIN-OR Clp")
OPENMS_CHECKLIB(CONTRIB_CBC4 "libCoinUtils;CoinUtils" "libCoinUtilsd;CoinUtils" "COIN-OR Utils")
OPENMS_CHECKLIB(CONTRIB_CBC5 "libOsi;Osi" "libOsid;Osi" "COIN-OR Osi")
OPENMS_CHECKLIB(CONTRIB_CBC6 "libOsiClp;OsiClp" "libOsiClpd;OsiClp" "COIN-OR OsiClp")
set(CONTRIB_CBC ${CONTRIB_CBC1} ${CONTRIB_CBC2} ${CONTRIB_CBC3} ${CONTRIB_CBC4} ${CONTRIB_CBC5} ${CONTRIB_CBC6} )

OPENMS_CHECKLIB(CONTRIB_ZLIB "zlib;z" "zlib_d;z" "ZLIB")
OPENMS_CHECKLIB(CONTRIB_BZIP2 "libbz2;bz2" "libbz2d;bz2" "BZIP2")

## IMS
if (MSVC)
	OPENMS_CHECKLIB(CONTRIB_IMSLIB "ims" "imsd" "imslib")
else()
	# if enabled stl-debug we need a imslib which is also built with stl-debug
	# "normal" debug mode is only supported manually
	if (STL_DEBUG)
		OPENMS_CHECKLIB(CONTRIB_IMSLIB "imsd" "imsd" "imslib")
	else()
		OPENMS_CHECKLIB(CONTRIB_IMSLIB "ims" "ims" "imslib")
	endif()
endif()

## ANDIMS
if (USE_ANDIMS)
	OPENMS_CHECKLIB(CONTRIB_ANDIMS "ms10lib.a" "ms10lib.a" "ms10lib")

	OPENMS_CHECKLIB(CONTRIB_NETCDF_CPP "libnetcdf_c++.a" "libnetcdf_c++.a" "netCDF_c++")
	OPENMS_CHECKLIB(CONTRIB_NETCDF_C "libnetcdf.a" "libnetcdf.a" "netCDF")
	set(CONTRIB_NETCDF ${CONTRIB_NETCDF_CPP} ${CONTRIB_NETCDF_C})
endif()

## libSVM
if (MSVC)
	## we created a lib here, because we cannot just add an object ...
	OPENMS_CHECKLIB(SVM_LIBRARY "svm" "svmd" "LibSVM")
else()
	## there is no FIND_OBJECT ... so just hard-code it (FIND_LIBRARY won't work)
	set(SVM_LIBRARY general ${CONTRIB_LIB_DIR}svm.o CACHE INTERNAL "object file of svm library") # do not distinguish between optimized and debug
endif()

if (MSVC)
	## needed to locate libs (put this above ADD_LIBRARY() - otherwise it will not work)
	link_directories(${CONTRIB_LIB_DIR})
endif()

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include/)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/include/) ## for configured files, e.g. config.h
INCLUDE_DIRECTORIES(${CONTRIB_INCLUDE_DIR})

## QT
SET(QT_USE_QTXML 1)
SET(QT_USE_QTNETWORK 1)
SET(QT_USE_QTSQL 1)
SET(QT_USE_QTOPENGL 1)
SET(QT_USE_QTSVG 1)
SET(QT_USE_QTTEST 1)
SET(QT_MIN_VERSION "4.5.0")
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})
INCLUDE(UseQt4)
## FIX
## - QT4 library list (VS2010 requires semicolon separated lists, previous versions use space, thus a lib named "opengl32.lib glu32.lib gdi32.lib user32.lib delayimp.lib" now poses a problem)
## - this fix will be obsolete with CMake 2.8.1
set(QT_LIBRARIES_TMP "")
foreach (qtlib ${QT_LIBRARIES})
	if (${qtlib} STREQUAL "opengl32.lib glu32.lib gdi32.lib user32.lib delayimp.lib")
		list(APPEND QT_LIBRARIES_TMP "opengl32.lib;glu32.lib;gdi32.lib;user32.lib;delayimp.lib")
		message(STATUS "Fixing QT library list... done")
	else()
		list(APPEND QT_LIBRARIES_TMP ${qtlib})
	endif()
endforeach()
set(QT_LIBRARIES ${QT_LIBRARIES_TMP})
## ENDFIX

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

if ("${SVNVERSION_EXECUTABLE}" STREQUAL "SVNVERSION_EXECUTABLE-NOTFOUND")
	message(STATUS "The programm svnversion could not be found. SVN-revision information will not be available! Add the location of svnversion(.exe) to your PATH environment variable if you require SVN-revision.")
	set(OPENMS_HAS_SVNVERSION OFF CACHE INTERNAL "SVNVersion(.exe) present?")
else()
	set(OPENMS_HAS_SVNVERSION ON CACHE INTERNAL "SVNVersion(.exe) present?")
endif()

## generate include header for version - it needs to be present!
## Warning: SVNRevisionUpdater(.exe) expects the line in openms_svn_revision.h to contain exactly 3(!!!) spaces up until the comment section "// ..."!
set(SVN_REVISION_FILE "${PROJECT_BINARY_DIR}/include/OpenMS/openms_svn_revision.h")
file(WRITE ${SVN_REVISION_FILE}
"#define OPENMS_SVN_REVISION \"exported\" // generated by the build system! do not edit this manually!
")

########################################################
###    configure config.h														 ###
########################################################

set(CMAKE_INSTALL_PREFIX "" CACHE INTERNAL "This must be equal to INSTALL_PREFIX!" FORCE)
set(INSTALL_PREFIX "." CACHE PATH "Installation prefix for OpenMS and TOPP")

## define some directories
if ("${INSTALL_PREFIX}" STREQUAL ".")
	set(CF_OPENMS_DATA_PATH ${PROJECT_SOURCE_DIR}/share/OpenMS CACHE INTERNAL "Path to the shared documents of OpenMS.")
	set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
else()
	set(CF_OPENMS_DATA_PATH ${INSTALL_PREFIX}/share/OpenMS CACHE INTERNAL "Path to the shared documents of OpenMS.")
	set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX})
endif()



set(CF_OPENMS_TEST_DATA_PATH ${PROJECT_SOURCE_DIR}/source/TEST/data/ CACHE INTERNAL "Path to the test data")

## check for Microsoft Visual Studio compiler
if (MSVC)
	set(OPENMS_COMPILER_MSVC "1" CACHE INTERNAL "Do we use Microsoft Compiler?")
endif()
## check for G++
if (CMAKE_COMPILER_IS_GNUCXX)
	set(OPENMS_COMPILER_GXX "1" CACHE INTERNAL "Do we use G++ Compiler?")
endif()

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(OPENMS_BIG_ENDIAN)

## check 32/64 bit architecture (defined above!)
if (NOT DEFINED OPENMS_64BIT_ARCHITECTURE)
	message(FATAL_ERROR "Cmake script was re-ordered and is now invalid! Please make sure that OPENMS_64BIT_ARCHITECTURE is defined when config.h.in is configured!")
endif()

include(CheckTypeSize) ## Check sizeof a type
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
	## search for another Int32 type
	if (SIZE_INT MATCHES "4")
		set(CF_OPENMS_INT32_TYPE int)
	elseif (SIZE_SHORT MATCHES "4")
		set(CF_OPENMS_INT32_TYPE short)
	elseif (SIZE_LONG MATCHES "4")
		set(CF_OPENMS_INT32_TYPE long)
	else()
		Message(FATAL_ERROR "Cannot find signed 32bit integer type. Please contact the developers!")
	endif()
endif()

CHECK_TYPE_SIZE("int64_t" SIZE_INT64)
if (HAVE_SIZE_INT64)
	set(CF_OPENMS_INT64_TYPE int64_t)
else()
	## search for another Int64 type
	if (SIZE_INT MATCHES "8")
		set(CF_OPENMS_INT64_TYPE int)
	elseif (SIZE_LONG MATCHES "8")
		set(CF_OPENMS_INT64_TYPE long)
	elseif (SIZE_LONGLONG MATCHES "8")
		set(CF_OPENMS_INT64_TYPE "long long")
	else()
		Message(FATAL_ERROR "Cannot find signed 64bit integer type. Please contact the developers!")
	endif()
endif()

CHECK_TYPE_SIZE("uint8_t" SIZE_UINT8)
if (HAVE_SIZE_UINT8)
	set(CF_OPENMS_BYTE_TYPE uint8_t)
else()
	## search for another uint8 type
	if (SIZE_UCHAR MATCHES "1")
		set(CF_OPENMS_BYTE_TYPE "unsigned char")
	elseif (SIZE_USHORT MATCHES "1")
		set(CF_OPENMS_BYTE_TYPE "unsigned short")
	else()
		Message(FATAL_ERROR "Cannot find unsigned 8bit integer (byte) type. Please contact the developers!")
	endif()
endif()


CHECK_TYPE_SIZE("uint64_t" SIZE_UINT64)
if (HAVE_SIZE_UINT64)
	set(CF_OPENMS_UINT64_TYPE uint64_t)
else()
	## search for another uint64 type
	if (SIZE_ULONG MATCHES "8")
		set(CF_OPENMS_UINT64_TYPE "unsigned long")
	elseif (SIZE_ULONGLONG MATCHES "8")
		set(CF_OPENMS_UINT64_TYPE "unsigned long long")
	else()
		Message(FATAL_ERROR "Cannot find uint64 type. Please contact the developers!")
	endif()
endif()

## system headers:
include(CheckIncludeFileCXX) ## Check if the include file exists.

CHECK_INCLUDE_FILE_CXX("unistd.h" OPENMS_HAS_UNISTD_H)
CHECK_INCLUDE_FILE_CXX("process.h" OPENMS_HAS_PROCESS_H)

CHECK_INCLUDE_FILE_CXX("time.h" OPENMS_HAS_TIME_H)
CHECK_INCLUDE_FILE_CXX("sys/types.h" OPENMS_HAS_SYS_TYPES_H)
CHECK_INCLUDE_FILE_CXX("sys/times.h" OPENMS_HAS_SYS_TIMES_H)
CHECK_INCLUDE_FILE_CXX("sys/time.h"  OPENMS_HAS_SYS_TIME_H)
CHECK_INCLUDE_FILE_CXX("stdint.h"  OPENMS_HAS_STDINT_H)

include(CheckFunctionExists)
## in MinGW we have the signal.h header, but no kill() as in Linux, so we need to check for the kill() function
CHECK_FUNCTION_EXISTS("kill" OPENMS_HAS_KILL)
CHECK_FUNCTION_EXISTS("sysconf" OPENMS_HAS_SYSCONF)

## ANDIMS flag USE_ANDIMS defined above!
if (NOT DEFINED USE_ANDIMS)
	message(FATAL_ERROR "Cmake script was re-ordered and is now invalid! Please make sure that USE_ANDIMS is defined when config.h.in is configured!")
endif()

## user flag with default "QMYSQL"
set(QT_DB_PLUGIN "QMYSQL" CACHE STRING "User switch to change the Qt database plugin.")

## Configures what is tested (DB-test, MS2-search engines)
include(cmake/TestConfig.cmake)

## are we building a shared or static lib?! (BOTH within the same BUILD-tree is NOT possible with OpenMS!!)
set(BUILD_SHARED_LIBS true)

## replace any variables in config.h.in with current values
set (CONFIGURED_CONFIG_H ${PROJECT_BINARY_DIR}/include/OpenMS/config.h)
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/config.h.in ${CONFIGURED_CONFIG_H})

## replace any variables in openms_package_version.h.in with current values
set (CONFIGURED_OPENMS_PACKAGE_VERSION_H ${PROJECT_BINARY_DIR}/include/OpenMS/openms_package_version.h)
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/openms_package_version.h.in ${CONFIGURED_OPENMS_PACKAGE_VERSION_H})

## add configured files (used in 'cmake/includes.cmake')
set (OPENMS_CONFIGURED_FILES ${SVN_REVISION_FILE} ${CONFIGURED_CONFIG_H} ${CONFIGURED_OPENMS_PACKAGE_VERSION_H})

########################################################
###    Documentation                                 ###
########################################################

include(cmake/doc.cmake)

########################################################
###    BUILD the lib                                 ###
########################################################

#package type
set(VALID_PACKAGE_TYPES "none" "rpm" "deb" "dmg")
set(PACKAGE_TYPE "none" CACHE STRING "Package type (internal): ${VALID_PACKAGE_TYPES}")
# check build type
list(FIND VALID_PACKAGE_TYPES ${PACKAGE_TYPE} list_pos)
if( ${list_pos} EQUAL -1 )
	message(STATUS "The PACKAGE_TYPE ${PACKAGE_TYPE} is invalid")
	message(STATUS "Valid PACKAGE_TYPEs are:")
	foreach( BT ${VALID_PACKAGE_TYPES} )
		message(STATUS " * ${BT}")
	endforeach()
	message(FATAL_ERROR "Aborting ...")
endif()


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(OPENMS_BINARY_DIR "${PROJECT_BINARY_DIR}/bin")
set(OPENMS_WIN32_DLL_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

## big include file for headers and C files, which fills the OpenMS_sources variable
include (cmake/includes.cmake)

if (MSVC)
	## use OpenMSd.dll in debug mode
	SET(CMAKE_DEBUG_POSTFIX d)
endif()

## add library target
## warning: set BUILD_SHARED_LIBS to decide if library is shared or static (see above)! We need the BUILD_SHARED_LIBS flag to set declspec flags for MSVC!
if (FOUND_CUDART AND MT_ENABLE_CUDA) # we need both conditions due to possible present cached entries
	CUDA_ADD_LIBRARY(OpenMS ${OpenMS_sources} ${Cuda_sources})
else()
	ADD_LIBRARY(OpenMS ${OpenMS_sources})
endif()

### OpenMS_GUI lib
ADD_LIBRARY(OpenMS_GUI ${OpenMSVisual_sources})
if (MSVC)
	set (GUI_lnk_flags "/FORCE:MULTIPLE /INCREMENTAL:NO /ignore:4006 /ignore:4088")
	set_target_properties(OpenMS_GUI PROPERTIES LINK_FLAGS ${GUI_lnk_flags}) ## allow multiple definitions of symbols (e.g. from template instanciations or STL-derived classes)
endif()	


## update the openms_svn_revision.h file before compiling OpenMS
if (OPENMS_HAS_SVNVERSION)
	# build the updater executable
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin/svn")
	add_executable(SVNRevisionUpdater ${PROJECT_SOURCE_DIR}/tools/svn/SVNRevisionUpdater.C)
	# need QtCore for RegEx
	target_link_libraries(SVNRevisionUpdater ${QT_LIBRARIES})
	get_target_property(SVNRevisionUpdater_LOCATION SVNRevisionUpdater LOCATION)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin") ## reset output dir

	# run SvnRevisionUpdater to update svn_revision.h
	add_custom_target(svn_revision_update COMMENT "update the openms_svn_revision.h with current svn revision number")
	add_dependencies(svn_revision_update SVNRevisionUpdater)
	add_custom_command(TARGET svn_revision_update
										 PRE_BUILD
										 COMMAND ${SVNRevisionUpdater_LOCATION} ${PROJECT_SOURCE_DIR} ${SVN_REVISION_FILE})


	add_dependencies(OpenMS svn_revision_update)
endif()


## all the dependency libraries are linked into libOpenMS.so, except Qt and CUDA which are still dynamic
set(OPENMS_DEP_LIBRARIES ${CONTRIB_CBC} ${CONTRIB_GSL} ${CONTRIB_GSLCBLAS} ${SVM_LIBRARY} ${CONTRIB_XERCESC} ${CONTRIB_ANDIMS} ${CONTRIB_NETCDF} ${CONTRIB_BOOST_IO} ${CONTRIB_IMSLIB} ${CONTRIB_BZIP2} ${CONTRIB_ZLIB})
if (TBB_FOUND)
	list(APPEND OPENMS_DEP_LIBRARIES ${TBB_LIBRARIES})
endif()
if (MSVC)
	list(APPEND OPENMS_DEP_LIBRARIES opengl32.lib)
endif()

target_link_libraries(OpenMS ${APPLE_EXTRA_LIBS} ${QT_LIBRARIES} ${OPENMS_DEP_LIBRARIES})
target_link_libraries(OpenMS_GUI OpenMS ${QT_LIBRARIES})

## OPENMS_LIBRARIES defines the libraries used by OpenMS; this should be used to link against executables
## somehow the link dependencies of cmake are broken such that when using POSTFIX names for libs the dependencies are not forwarded...
## we fix this by adding the dependencies of OpenMS directly to the executables as well.
set(OPENMS_LIBRARIES optimized OpenMS${CMAKE_RELEASE_POSTFIX} debug OpenMS${CMAKE_DEBUG_POSTFIX} ${QT_LIBRARIES})
if (MSVC) ## Windows
	list(APPEND OPENMS_LIBRARIES ${OPENMS_DEP_LIBRARIES})
endif()

## link libraries for GUI Tools
set(OPENMS_GUI_LIBRARIES optimized OpenMS_GUI${CMAKE_RELEASE_POSTFIX} debug OpenMS_GUI${CMAKE_DEBUG_POSTFIX} OpenMS ${QT_LIBRARIES})

## directory for OpenMS(d).lib (required for linking executables)
if (MSVC)
	link_directories(${PROJECT_BINARY_DIR})
endif()

####### TOPP #########
set(TOPP_DIR source/APPLICATIONS/TOPP)
include(${TOPP_DIR}/executables.cmake)
foreach(i ${TOPP_executables})
	add_executable(${i} ${TOPP_DIR}/${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
	if (OPENMP_FOUND AND NOT MSVC)
		set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif()
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

##### GUI tools ########
set(GUI_DIR source/VISUAL/APPLICATIONS/GUITOOLS)
include(${GUI_DIR}/executables.cmake)
foreach(i ${GUI_executables})
	set(resource_file ${CMAKE_SOURCE_DIR}/${GUI_DIR}/${i}.rc)
	set(resource_dir ${CMAKE_SOURCE_DIR}/${GUI_DIR}/${i}-resources/)
	## add icons to TOPPView and INIFileEditor
	if (MSVC AND EXISTS ${resource_file})
		Message(STATUS "Setting resource file ${resource_file} for ${i}")
		add_executable(${i} ${GUI_DIR}/${i}.C ${resource_file})
	elseif(APPLE AND EXISTS ${resource_dir})
		add_executable(${i} ${GUI_DIR}/${i}.C)
		# extra work for pseudo-mac-bundle
		# move prepared package to RUNTIME_OUTPUT
		file(COPY "${PROJECT_SOURCE_DIR}/cmake/MacOSX/${i}-resources/${i}.app"
		DESTINATION "${PROJECT_BINARY_DIR}/bin")
	else()
		add_executable(${i} ${GUI_DIR}/${i}.C)
	endif()

	## append visual lib as dependency for GUI tools
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
	if (OPENMP_FOUND AND NOT MSVC)
		set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif()
endforeach(i)
add_custom_target(GUI)
add_dependencies(GUI ${GUI_executables})

## some regular TOPP/UTILS need the GUI lib
set(executables_withGUIlibDep ${TOPP_executables_with_GUIlib} ${UTILS_executables_with_GUIlib})
foreach(i ${executables_withGUIlibDep}) 
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
endforeach(i)

	
##### Doc progs ####
set(OpenMS_doc_executables)
include(doc/doxygen/parameters/executables.cmake)
foreach(i ${OpenMS_doc_executables})
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "doc/doxygen/parameters")
	add_executable(${i} EXCLUDE_FROM_ALL doc/doxygen/parameters/${i}.C)
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")
endforeach(i)
add_custom_target(doc_progs)
add_dependencies(doc_progs ${OpenMS_doc_executables} TOPP UTILS)

##### Examples #####
if ("${PACKAGE_TYPE}" STREQUAL "none")
	set(EXAMPLES_executables)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
	include(source/EXAMPLES/executables.cmake)
	foreach(i ${EXAMPLES_executables})
	  add_executable(${i} source/EXAMPLES/${i}.C)
	  target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
	endforeach(i)
	add_custom_target(Tutorials_build)
	add_dependencies(Tutorials_build ${EXAMPLES_executables})

	if (MSVC_IDE)
		set(EXAMPLES_EXE_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$(OutDir))
	else()
		set(EXAMPLES_EXE_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
	endif()


	## Create single target for each example, add it as a dependency to the main target
	foreach(i ${EXAMPLES_executables})
		list(APPEND EX_LIST COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/source/EXAMPLES ${EXAMPLES_EXE_DIR}/${i})
	endforeach(i)
	list(APPEND EX_LIST "echo all done")
	add_custom_target(Tutorials_exec
		COMMENT "Executes the examples"
	)
	separate_arguments(EX_LIST)
	add_custom_command(TARGET Tutorials_exec
                    POST_BUILD
										${EX_LIST}
										)

endif()

##### Tests #####
if ("${PACKAGE_TYPE}" STREQUAL "none")
	INCLUDE(Dart) ## for Nighlty Build log
	add_subdirectory(source/TEST EXCLUDE_FROM_ALL) ## configure the Test project (but do not add it to OpenMS itself)
endif()

#### CPP CHECK TESTS #### 
find_package(cppcheck)
if( CPPCHECK_FOUND )
	include( cmake/cppcheck.cmake )
endif()

##### targets list #####
if (MSVC)
	add_custom_target(targets
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "The following make targets are available:"
										COMMAND ${CMAKE_COMMAND} -E echo "    ALL_BUILD       [Visual Studio only] builds the OpenMS library, TOPP tools and UTILS tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    [no target]     [NMake only]         builds the OpenMS library, TOPP tools and UTILS tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    OpenMS          builds the OpenMS library"
										COMMAND ${CMAKE_COMMAND} -E echo "    TOPP            builds the TOPP tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    UTILS           builds the UTILS tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    GUI             builds the GUI tools (TOPPView,...)"
										COMMAND ${CMAKE_COMMAND} -E echo "    RUN_TESTS       [Visual Studio only] executes OpenMS and TOPP tests (*)"
										COMMAND ${CMAKE_COMMAND} -E echo "    test            [NMake only]         executes OpenMS and TOPP tests (*)"
										COMMAND ${CMAKE_COMMAND} -E echo "                    *) make sure they are built using the 'test_build' target (see below)"
										COMMAND ${CMAKE_COMMAND} -E echo "    Tutorials_build builds the tutorials in source/EXAMPLES"
										COMMAND ${CMAKE_COMMAND} -E echo "    Tutorials_exec  executes the tutorials in source/EXAMPLES"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc             builds the doxygen documentation and tutorials"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc_tutorials   builds the pdf tutorials"
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "Single TOPP tools and UTILS have their own target, e.g. TOPPView"
										COMMAND ${CMAKE_COMMAND} -E echo "The class tests have their own project in ./source/TEST (project test_build)."
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "The most important targets for OpenMS"
										VERBATIM
		)
else()
	add_custom_target(targets
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "The following make targets are available:"
										COMMAND ${CMAKE_COMMAND} -E echo "    [no target]     builds the OpenMS library, TOPP tools and UTILS tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    OpenMS          builds the OpenMS library"
										COMMAND ${CMAKE_COMMAND} -E echo "    TOPP            builds the TOPP tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    UTILS           builds the UTILS tools"
										COMMAND ${CMAKE_COMMAND} -E echo "    GUI             builds the GUI tools (TOPPView,...)"
										COMMAND ${CMAKE_COMMAND} -E echo "    test_build      builds the OpenMS tests"
										COMMAND ${CMAKE_COMMAND} -E echo "    test            executes OpenMS and TOPP tests"
										COMMAND ${CMAKE_COMMAND} -E echo "                    make sure they are built using the 'test_build' target"
										COMMAND ${CMAKE_COMMAND} -E echo "    Tutorials_build builds the tutorials in source/EXAMPLES"
										COMMAND ${CMAKE_COMMAND} -E echo "    Tutorials_exec  executes the tutorials in source/EXAMPLES"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc             builds the doxygen documentation and tutorials"
										COMMAND ${CMAKE_COMMAND} -E echo "    doc_tutorials   builds the pdf tutorials"
										COMMAND ${CMAKE_COMMAND} -E echo "    help            list all available targets (very long)"
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "Single TOPP tools and UTILS have their own target, e.g. TOPPView"
										COMMAND ${CMAKE_COMMAND} -E echo "The class tests have their own project in ./source/TEST."
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMAND ${CMAKE_COMMAND} -E echo "=========================================================================="
										COMMAND ${CMAKE_COMMAND} -E echo ""
										COMMENT "The most important targets for OpenMS"
										VERBATIM
		)
endif()

##### Message after OpenMS has been built #####
if (MSVC)
	## copy OpenMS.dll to test executables dir
	get_target_property(WIN32_DLLLOCATION OpenMS LOCATION)
	get_filename_component(WIN32_DLLPATH ${WIN32_DLLLOCATION} PATH)

	## copy OpenMS.dll to test executables dir "$(TargetFileName)" is a placeholder filled by VS at runtime
	file(TO_NATIVE_PATH "${WIN32_DLLPATH}/$(TargetFileName)" DLL_SOURCE)
	file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}/source/TEST/bin/$(ConfigurationName)/$(TargetFileName)" DLL_TARGET)
	# create target path if not exists
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
                    COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_BINARY_DIR}/source/TEST/bin/$(ConfigurationName)"
					COMMAND copy ${DLL_SOURCE} ${DLL_TARGET} /Y
                    COMMENT "Library build. Copying OpenMS(d).dll to test binary dir"
                    VERBATIM)

	## copy OpenMS_GUI.dll
	add_custom_command(TARGET OpenMS_GUI
                    POST_BUILD
					COMMAND copy ${DLL_SOURCE} ${DLL_TARGET} /Y
                    COMMENT "GUI Library build. Copying OpenMS_GUI(d).dll to test binary dir"
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

######################################################
#### build settings and configs for external code ####
######################################################
if ("${PACKAGE_TYPE}" STREQUAL "none")
	## grep definitions (-D<XXX>) that were set using 'add_definitions()'
	get_property(OPENMS_CMP_DEFS DIRECTORY . PROPERTY COMPILE_DEFINITIONS)
	## convert to CXX flags
	foreach(i ${OPENMS_CMP_DEFS})
		set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} -D${i}")
	endforeach(i)
	## for add_definition commands (currently empty, but who knows...)
	set(CF_OPENMS_ADD_DEFINITIONS)
	## the filename for library settings (this file will be included in OpenMSConfig.cmake)
	set(CF_LibOpenMSExport ${PROJECT_BINARY_DIR}/cmake/OpenMSLibraryExport.cmake)

	## create code that allows to call find_package()
	## (must come after the above code that sets the variables!)
	get_directory_property(OPENMS_INCLUDE_DIRS  INCLUDE_DIRECTORIES)
	configure_file(
		"${PROJECT_SOURCE_DIR}/cmake/OpenMSConfig.cmake.in"
		"${PROJECT_BINARY_DIR}/cmake/OpenMSConfig.cmake"
		@ONLY
	)

	## Create OpenMSLibraryExport.cmake:
	## Note: export_build_settings and
	##       export_library_dependencies are deprecated and should not be used
	if ("${INSTALL_PREFIX}" STREQUAL ".")
		set(CVERSION "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}")
		if (CVERSION VERSION_LESS "2.8.0")
			message(STATUS "Cannot register OpenMS with CMake! For external code, set the path to OpenMS during find_package() manually.")
		else()
			## register OpenMS with CMake so that it can be found easily
			export(PACKAGE OpenMS)
		endif()
		## export our lib
		export(TARGETS OpenMS FILE ${CF_LibOpenMSExport})
	else()
		include(cmake/install.cmake)
	endif()

else()
######################################################
####          install/copy (cpack mainly)         ####
######################################################
	## general definitions
	set(CPACK_PACKAGE_NAME "OpenMS")
	set(CPACK_PACKAGE_VENDOR "OpenMS.de")
	set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMS - A framework for mass spectrometry")
	set(CPACK_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION_MAJOR}.${OPENMS_PACKAGE_VERSION_MINOR}.${OPENMS_PACKAGE_VERSION_PATCH}")
	set(CPACK_PACKAGE_VERSION_MAJOR "${OPENMS_PACKAGE_VERSION_MAJOR}")
	set(CPACK_PACKAGE_VERSION_MINOR "${OPENMS_PACKAGE_VERSION_MINOR}")
	set(CPACK_PACKAGE_VERSION_PATCH "${OPENMS_PACKAGE_VERSION_PATCH}")
	set(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMS-${CPACK_PACKAGE_VERSION}")
	set(CPACK_PACKAGE_DESCRIPTION_FILE ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageDescriptionFile.cmake)
	set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/License.txt)
	set(CPACK_RESOURCE_FILE_WELCOME ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceWelcomeFile.txt)
	set(CPACK_RESOURCE_FILE_README ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceReadme.txt)

	# install routines for MacOSX
	if("${PACKAGE_TYPE}" STREQUAL "dmg")
		include(cmake/package_dragndrop_dmg.cmake)
	endif()

	if("${PACKAGE_TYPE}" STREQUAL "rpm")
		include(cmake/package_rpm.cmake)
	endif()

	if("${PACKAGE_TYPE}" STREQUAL "deb")
		include(cmake/package_deb.cmake)
	endif()

endif()

