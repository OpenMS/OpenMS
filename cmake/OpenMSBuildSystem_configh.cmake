
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

## user flag with default "QMYSQL" (put in config.h)
set(CF_QT_DB_PLUGIN "QMYSQL" CACHE STRING "User switch to change the Qt database plugin.")

## replace any variables in config.h.in with current values
set (CONFIGURED_CONFIG_H ${PROJECT_BINARY_DIR}/include/OpenMS/config.h)
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/config.h.in ${CONFIGURED_CONFIG_H})

## replace any variables in openms_package_version.h.in with current values
set (CONFIGURED_OPENMS_PACKAGE_VERSION_H ${PROJECT_BINARY_DIR}/include/OpenMS/openms_package_version.h)
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/openms_package_version.h.in ${CONFIGURED_OPENMS_PACKAGE_VERSION_H})
