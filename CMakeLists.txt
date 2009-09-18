# CMake installation file
# Created 2008/10 during OpenMS retreat by CB && AB

PROJECT("OpenMS")

# TODO LIST:
# - test for contrib headers (CONTRIB_HEADER_DIR ccmake variable - test one header for each required library)
# - parallel build of test executables
# - target to execute all EXAMPLES

########################################################
### entries meant to be configured using CMake cache ###
### - do NOT hardcode them here!										 ###
### - edit them within CMakeCache.txt using ccmake	 ###
########################################################
# OPENMS_BUILD_TYPE
# CONTRIB_CUSTOM_DIR
# STL_DEBUG
# QT_DB_PLUGIN
# MT_CUDA_BUILD_TYPE

########################################################
###    manual entries (edit this for new release)    ###
########################################################

set(OPENMS_PACKAGE_VERSION_MAJOR "1")
set(OPENMS_PACKAGE_VERSION_MINOR "5")
set(OPENMS_PACKAGE_VERSION_PATCH "0")

########################################################
###    end manual	entries													   ###
########################################################

# Heart of the BUILD system : only edit when you know what you´re doing (we don´t)
# quick manual for most commands: http://www.cmake.org/cmake/help/cmake2.6docs.html
# useful predefined variables: http://www.paraview.org/Wiki/CMake_Useful_Variables

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
set(CF_OPENMS_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION_MAJOR}.${OPENMS_PACKAGE_VERSION_MINOR}.${OPENMS_PACKAGE_VERSION_PATCH}" CACHE INTERNAL "OpenMS VERSION" FORCE)

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

if (CMAKE_COMPILER_IS_GNUCXX)
	if (MT_ENABLE_CUDA)
		add_definitions(-Wall -Wextra -Wno-non-virtual-dtor -Wno-long-long) # necessary since CUDA contains non-pedantic code
	else()
		add_definitions(-Wall -Wextra -Wno-non-virtual-dtor -Wno-long-long --pedantic)
	endif()

	## Recommended setting for eclipse, see http://www.cmake.org/Wiki/CMake:Eclipse
	if (CMAKE_GENERATOR STREQUAL "Eclipse CDT4 - Unix Makefiles")
		add_definitions(-fmessage-length=0)
	endif()

	## CGAL needs extra definitions to avoid pre and postconditions 
	## (see mails on open-ms-general of Hendrik Weisser. Thx very much!)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS -DCGAL_NO_ASSERTIONS")


elseif (MSVC)
	#add_definitions(/Wall) ## disable for now.. its just too much!

	## disable dll-interface warning
	add_definitions(/wd4251 /wd4275)
	
	## disable deprecated functions warning (e.g. for POSIX functions)
	add_definitions(/wd4996)

	## disable explicit template instantiation request for partially defined classes
	add_definitions(/wd4661)

	
	## don't warn about unchecked std::copy()
	add_definitions(/D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_DEPRECATE)

	## disable min and max macros by default
	## see CGAL and BOOST configs for more elaborate explanations
	add_definitions(/DNOMINMAX)

	## xerces bug workaround
	add_definitions(/DOPENMS_XERCESDLL)
	
	## compile such that GSL is in DLL mode
	#add_definitions(/DGSL_DLL)

	## FeatureFinder.obj is huge and won't compile in VS2008 debug otherwise:
	add_definitions(/bigobj)
	
	## minimal code rebuild
	#add_definitions(/Gm /Zi)

else()
	set(CMAKE_COMPILER_IS_INTELCXX true CACHE INTERNAL "Is Intel C++ compiler (icpc)")
endif()


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
		set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_CXX_FLAGS})
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
  message(STATUS "If you want to set your own path to contrib use -D CONTRIB_CUSTOM_DIR:PATH=\"/my/path/to/contrib\"")
else()
	set(CONTRIB_DIR ${CONTRIB_CUSTOM_DIR})
	message(STATUS "Using user defined contrib directory: ${CONTRIB_DIR}")
endif()

FIND_PATH(CONTRIB_CONFIRM_DIR README_contrib1.5.txt PATHS ${CONTRIB_DIR} NO_DEFAULT_PATH)

if ("${CONTRIB_CONFIRM_DIR}" STREQUAL "CONTRIB_CONFIRM_DIR-NOTFOUND")
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

OPENMS_CHECKLIB(CONTRIB_CGAL "cgal-vc80-mt;cgal-vc90-mt;CGAL" "cgal-vc80-mt-gd;cgal-vc90-mt-gd;CGAL" "CGAL")

OPENMS_CHECKLIB(CONTRIB_CBC1 "libCbc;Cbc" "libCbcd;Cbc" "COIN-OR Cbc")
OPENMS_CHECKLIB(CONTRIB_CBC2 "libCgl;Cgl" "libCgld;Cgl" "COIN-OR Cgl")
OPENMS_CHECKLIB(CONTRIB_CBC3 "libClp;Clp" "libClpd;Clp" "COIN-OR Clp")
OPENMS_CHECKLIB(CONTRIB_CBC4 "libCoinUtils;CoinUtils" "libCoinUtilsd;CoinUtils" "COIN-OR Utils")
OPENMS_CHECKLIB(CONTRIB_CBC5 "libOsi;Osi" "libOsid;Osi" "COIN-OR Osi")
OPENMS_CHECKLIB(CONTRIB_CBC6 "libOsiClp;OsiClp" "libOsiClpd;OsiClp" "COIN-OR OsiClp")

# if enabled stl-debug we need a imslib which is also built with stl-debug
# "normal" debug mode is only supported manually
if (STL_DEBUG)
	OPENMS_CHECKLIB(CONTRIB_IMSLIB "imsd" "imsd" "imslib")
else()
	if (MSVC)
		OPENMS_CHECKLIB(CONTRIB_IMSLIB "ims" "imsd" "imslib")
	else()
		OPENMS_CHECKLIB(CONTRIB_IMSLIB "ims" "ims" "imslib")
	endif()
endif()

set(CONTRIB_CBC ${CONTRIB_CBC1} ${CONTRIB_CBC2} ${CONTRIB_CBC3} ${CONTRIB_CBC4} ${CONTRIB_CBC5} ${CONTRIB_CBC6} )

if (NOT MSVC AND USE_ANDIMS)
	OPENMS_CHECKLIB(CONTRIB_ANDIMS "ms10lib.a" "ms10lib.a" "ms10lib")

	OPENMS_CHECKLIB(CONTRIB_NETCDF_CPP "libnetcdf_c++.a" "libnetcdf_c++.a" "netCDF_c++")
	OPENMS_CHECKLIB(CONTRIB_NETCDF_C "libnetcdf.a" "libnetcdf.a" "netCDF")
	set(CONTRIB_NETCDF ${CONTRIB_NETCDF_CPP} ${CONTRIB_NETCDF_C})
endif()

## libSVM support
if (MSVC) 
	## we created a lib here, because we cannot just add an object..
	OPENMS_CHECKLIB(SVM_LIBRARY "svm" "svmd" "LibSVM")
else()
	## there is no FIND_OBJECT.. so just hard-code it (FIND_LIBRARY won't work)
	set(SVM_LIBRARY general ${CONTRIB_LIB_DIR}svm.o CACHE INTERNAL "object file of svm library") # do not distinguish between optimized and debug
endif()

if (MSVC)
	## needed to locate libs (put this above ADD_LIBRARY() - otherwise it will not work)
	link_directories(${CONTRIB_LIB_DIR})
endif()


INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/include/) ## for configured files, i.e. config.h
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include/)
INCLUDE_DIRECTORIES(${CONTRIB_INCLUDE_DIR})

## QT
SET(QT_USE_QTXML 1)
SET(QT_USE_QTNETWORK 1)
SET(QT_USE_QTSQL 1)
SET(QT_USE_QTOPENGL 1)
SET(QT_MIN_VERSION "4.3.0")
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

if ("${SVNVERSION_EXECUTABLE}" STREQUAL "SVNVERSION_EXECUTABLE-NOTFOUND")
	message(STATUS "The programm svnversion could not be found. SVN-revision information will not be available! Add the location of svnversion(.exe) to your PATH environment variable if you require SVN-revision.")
	set(OPENMS_HAS_SVNVERSION OFF CACHE INTERNAL "SVNVersion(.exe) present?")
else()
	set(OPENMS_HAS_SVNVERSION ON CACHE INTERNAL "SVNVersion(.exe) present?")
endif()

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

## DB test settings
set(DB_TEST OFF CACHE BOOL "If true, the DB tests are enabled.")
set(DB_TEST_HOST "localhost" CACHE STRING "Test database server name (only used if DB_TEST is true).")
set(DB_TEST_PORT "3307" CACHE STRING "Test database server port (only used if DB_TEST is true).")
set(DB_TEST_DB "OPENMS_TEST_DB" CACHE STRING "Test database name (only used for DB_TEST).")
set(DB_TEST_USER "openms_test_user" CACHE STRING "Test database user name (only used for DB_TEST).")
set(DB_TEST_PW "openms_test_password" CACHE STRING "Test database user password (only used for DB_TEST).")
if (DB_TEST)
	# output
	message(STATUS "DB testing enabled")
	message(STATUS "DB testing - creating credentials files ...")
	# create assorted credentials files
	configure_file(${PROJECT_SOURCE_DIR}/source/TEST/DB_credentials.txt.in ${PROJECT_BINARY_DIR}/source/TEST/DB_credentials.txt)
	configure_file(${PROJECT_SOURCE_DIR}/source/TEST/TOPP/DBImporter.ini.in ${PROJECT_BINARY_DIR}/source/TEST/TOPP/DBImporter.ini)
	configure_file(${PROJECT_SOURCE_DIR}/source/TEST/TOPP/DBExporter.ini.in ${PROJECT_BINARY_DIR}/source/TEST/TOPP/DBExporter.ini)
endif()

########################################################
###    Documentation    														 ###
########################################################

include(cmake/doc.cmake)

########################################################
###    BUILD the lib																 ###
########################################################

#package type
set(PACKAGE_TYPE "none" CACHE STRING "Package type (internal): 'none', 'rpm', 'deb', 'dmg' or 'nsis'")

## are we building a shared or static lib?! (BOTH within the same BUILD-tree is NOT possible with OpenMS!!)
set(BUILD_SHARED_LIBS true)

## replace any variables in config.h.in with current values
set (CONFIGURED_CONFIG_H ${PROJECT_BINARY_DIR}/include/OpenMS/config.h)
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/config.h.in ${CONFIGURED_CONFIG_H})

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

## always create the svn version file before compiling OpenMS
if (OPENMS_HAS_SVNVERSION)
## TODO replace this with PreProcessorMacro when compiling lib
## e.g. SET_SOURCE_FILE_PROPERTIES(libSource.c main.c COMPILE_FLAGS -DREVISON=234)
#	Message(STATUS "Adding dependency on SVN-Version target file...")
#	add_dependencies(OpenMS svnversion_file_target)
endif()


## all the dependency libraries are linked into libOpenMS.so, except Qt and CUDA which are still dynamic
set(OPENMS_DEP_LIBRARIES ${CONTRIB_CBC} ${CONTRIB_CGAL} ${CONTRIB_GSL} ${CONTRIB_GSLCBLAS} ${SVM_LIBRARY} ${CONTRIB_XERCESC} ${CONTRIB_ANDIMS} ${CONTRIB_NETCDF} ${CONTRIB_BOOST_IO} ${CONTRIB_IMSLIB})
if (TBB_FOUND)
	set (OPENMS_DEP_LIBRARIES ${OPENMS_DEP_LIBRARIES} ${TBB_LIBRARIES})
endif()

if (MSVC)
	list(APPEND OPENMS_DEP_LIBRARIES opengl32.lib)
endif()

if (NOT APPLE AND NOT MSVC) ##=Linux
	target_link_libraries(OpenMS ${OPENMS_DEP_LIBRARIES}) ## do not link OpenMS against QT on Linux
else()
	target_link_libraries(OpenMS ${QT_LIBRARIES} ${OPENMS_DEP_LIBRARIES})
endif()

## OPENMS_LIBRARIES defines the libraries used by OpenMS; this should be used to link against executables
set(OPENMS_LIBRARIES OpenMS)
if (MSVC) ## Windows
	## somehow the link dependencies of cmake are broken such that when using POSTFIX names for libs the dependencies are not forwarded...
	## we fix this by adding the dependencies of OpenMS directly to the executables as well.
	set(OPENMS_LIBRARIES optimized OpenMS${CMAKE_RELEASE_POSTFIX} debug OpenMS${CMAKE_DEBUG_POSTFIX} ${QT_LIBRARIES} ${OPENMS_DEP_LIBRARIES})
elseif (NOT APPLE) ## Linux
	set(OPENMS_LIBRARIES optimized OpenMS${CMAKE_RELEASE_POSTFIX} debug OpenMS${CMAKE_DEBUG_POSTFIX} ${QT_LIBRARIES})
endif()

## directory for OpenMS(d).lib (required for linking executables)
if (MSVC)
	link_directories(${PROJECT_BINARY_DIR})
endif()

####### TOPP #########
include(source/APPLICATIONS/TOPP/executables.cmake)

set(TOPP_DIR source/APPLICATIONS/TOPP/)
foreach(i ${TOPP_executables})
	set(resource_file ${CMAKE_SOURCE_DIR}/${TOPP_DIR}${i}.rc)
	set(resource_dir ${CMAKE_SOURCE_DIR}/${TOPP_DIR}/${i}-resources/)
	## add icons to TOPPView and INIFileEditor
	if (MSVC AND EXISTS ${resource_file})
		Message(STATUS "Setting resource file ${resource_file} for ${i}")
		add_executable(${i} ${TOPP_DIR}${i}.C ${resource_file})
  elseif(APPLE AND EXISTS ${resource_dir}) 
    set_source_files_properties(${resource_dir}/${i}.icns 
      PROPERTIES
      MACOSX_PACKAGE_LOCATION Resources)

    add_executable(${i} MACOSX_BUNDLE ${TOPP_DIR}${i}.C ${resource_dir}/${i}.icns)         
    ## transform TOPPView and INIFileEditor to real Mac OS X applications
    set_target_properties(${i} PROPERTIES MACOSX_BUNDLE_INFO_STRING "${i} - Version 1.5")
    set_target_properties(${i} PROPERTIES MACOSX_BUNDLE_SHORT_VERSION_STRING "1.5")
    set_target_properties(${i} PROPERTIES MACOSX_BUNDLE_VERSION "1.5")
    set_target_properties(${i} PROPERTIES MACOSX_BUNDLE_LONG_VERSION_STRING "1.5")
    set_target_properties(${i} PROPERTIES MACOSX_BUNDLE_ICON_FILE "${i}.icns")
    set_target_properties(${i} PROPERTIES MACOSX_BUNDLE_GUI_IDENTIFIER "de.openms.${i}")
    set_target_properties(${i} PROPERTIES MACOSX_BUNDLE_BUNDLE_NAME "${i}") 
  else()
		add_executable(${i} ${TOPP_DIR}${i}.C)
	endif()
	
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

##### Doc progs ####
set(OpenMS_doc_executables)
include(doc/doxygen/parameters/executables.cmake)
foreach(i ${OpenMS_doc_executables})
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "doc/doxygen/parameters")
	add_executable(${i} EXCLUDE_FROM_ALL doc/doxygen/parameters/${i}.C)
	target_link_libraries(${i} ${OPENMS_LIBRARIES})
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
	  target_link_libraries(${i} ${OPENMS_LIBRARIES})
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


#### configs for external include files, e.g. to develop code with OpenMS
configure_file(
	"${PROJECT_SOURCE_DIR}/cmake/OpenMSConfig.cmake.in"
	"${PROJECT_BINARY_DIR}/cmake/OpenMSConfig.cmake"
	@ONLY
)


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
	file(TO_NATIVE_PATH "${WIN32_DLLPATH}/$(TargetFileName)" DLL_SOURCE)
	file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}/source/TEST/bin/$(OutDir)/$(TargetFileName)" DLL_TARGET)
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
                    COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_BINARY_DIR}/source/TEST/bin/$(OutDir)/"
										COMMAND copy ${DLL_SOURCE} ${DLL_TARGET} /Y                    
                    COMMENT "message after library is built & copy OpenMS(d).dll to test binary dir"
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

##### install/copy (cpack mainly) #####
## TODO:
##  - this code breaks the windows cmake call (already code for win support added .. needs testing)

if (NOT "${PACKAGE_TYPE}" STREQUAL "none")

	# general definitions
	set(CPACK_PACKAGE_NAME "OpenMS")
	set(CPACK_PACKAGE_VENDOR "OpenMS.de")
	set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMS - A framework for mass spectrometry")
	set(CPACK_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION_MAJOR}.${OPENMS_PACKAGE_VERSION_MINOR}.${OPENMS_PACKAGE_VERSION_PATCH}")
	set(CPACK_PACKAGE_VERSION_MAJOR "${OPENMS_PACKAGE_VERSION_MAJOR}")
	set(CPACK_PACKAGE_VERSION_MINOR "${OPENMS_PACKAGE_VERSION_MINOR}")
	set(CPACK_PACKAGE_VERSION_PATCH "${OPENMS_PACKAGE_VERSION_PATCH}")
	set(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMS-${CPACK_PACKAGE_VERSION}")
	set(CPACK_PACKAGE_DESCRIPTION_FILE ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageDescriptionFile.cmake)
	if (APPLE)
		set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/cmake/MacOSX/License.txt)
	else()
		set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/License.txt)
	endif()
	set(CPACK_RESOURCE_FILE_WELCOME ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceWelomeFile.txt)
	set(CPACK_RESOURCE_FILE_README ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceReadme.txt)
	
	# install routines for MacOSX
	if("${PACKAGE_TYPE}" STREQUAL "dmg")
		#include(cmake/package_dmg.cmake)
		include(cmake/package_dragndrop_dmg.cmake)
	endif()

	if("${PACKAGE_TYPE}" STREQUAL "rpm")
		include(cmake/package_rpm.cmake)
	endif()

	if("${PACKAGE_TYPE}" STREQUAL "deb")
		include(cmake/package_deb.cmake)
	endif()

else()
	## grep definitions that were set using 'add_definitions()'
	get_property(OPENMS_CMP_DEFS DIRECTORY . PROPERTY COMPILE_DEFINITIONS)
	
	## convert to CXX flags
	set(OPENMS_ADDCXX_FLAGS)
	#
	foreach(i ${OPENMS_CMP_DEFS})
		set(OPENMS_ADDCXX_FLAGS "${OPENMS_ADDCXX_FLAGS} -D${i}")
	endforeach(i)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMS_ADDCXX_FLAGS}")
	#Message(STATUS "Definitions: ${OPENMS_ADDCXX_FLAGS}")
	
	# build settings and  configs for external code
	include(CMakeExportBuildSettings)
	cmake_export_build_settings(${PROJECT_BINARY_DIR}/cmake/OpenMSBuildSettings.cmake)
	export_library_dependencies(${PROJECT_BINARY_DIR}/cmake/OpenMSLibraryDepends.cmake)
endif()
