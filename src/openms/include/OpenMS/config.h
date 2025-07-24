// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// ==========================================================================
//
// IMPORTANT:
// This is config.h / config.h.in
// Please do ONLY change config.h.in, as config.h is automagically
// created by CMAKE from config.h.in
//
// Use appropriate options to configure instead of changing config.h.
// Changes made in config.h will be lost as soon as you call CMAKE again.
//
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONFIG_H
#define OPENMS_CONFIG_H

// include OPENMS_DLLAPI macros
#include <OpenMS/OpenMSConfig.h>

#include <boost/current_function.hpp>

// Here are some global configuration flags for OpenMS

// Define compiler specifics (used in VERY few places only)
// Microsoft Visual Studio .NET, 2005, 2008
/* #undef OPENMS_COMPILER_MSVC */
// GNU g++
#define OPENMS_COMPILER_GXX

// __PRETTY_FUNCTION__ is a GNU G++ extension so we use the alternative in boost
#define OPENMS_PRETTY_FUNCTION BOOST_CURRENT_FUNCTION

// __FILENAME__ prints the path that is used during compilation of the object. Per default in CMake it is the full absolute path.
// This is confusing for debug logs from executables on user machines and unnecessary.
// This shifts the pointer to the char array of the filepath to the filename only.
// WARNING: Does not work for files outside of the src directory.
// Needs the pragmas for CLangs superstrict compiler warnings. But the pragmas lead to failures in Win and Lnx.
// Therefore it is not used for now.

/*
#define OPENMS_SHORT_FILE _Pragma("GCC diagnostic push") \
_Pragma("GCC diagnostic push") \
_Pragma("GCC diagnostic ignored \"-Wstring-plus-int\"") \
( __FILE__ + 19 ) \
_Pragma("GCC diagnostic pop")
*/

using cstr = const char *;

static constexpr cstr past_last_slash(cstr str, cstr last_slash)
{
  return
    *str == '\0' ? last_slash :
    *str == '/'  ? past_last_slash(str + 1, str + 1) :
    *str == '\\'  ? past_last_slash(str + 1, str + 1) :
    past_last_slash(str + 1, last_slash);
}

static constexpr cstr past_last_slash(cstr str)
{
  return past_last_slash(str, str);
}

// OPENMS_ASSERTIONS enables some debugging methods within some OpenMS classes
#ifdef OPENMS_COMPILER_MSVC
// we define this using NDEBUG on MSVC as there are multiple build types simultaneously in the Solution file,
// thus setting one value will not fit them all
#	ifndef NDEBUG // hopefully defined automatically by MS-compiler in Debug Mode
#		define OPENMS_ASSERTIONS
#	endif
#else // linux & Co (only one build type at a time)
#  if (1)
#    define OPENMS_ASSERTIONS
#  endif
#endif


// let Cmake decide if we are using Windows (i.e. if windows.h is available).
//   GCC and MSVC have pre-defined macros for this as well but using -ansi in GCC will disable those, thus asking the compiler is fragile
#ifndef WIN32  //avoid warning of redefinition
/* #undef WIN32 */
#endif

#ifdef WIN32   //should be true on: MinGW (32bit + 64bit) & MSVS compiler
#define OPENMS_WINDOWSPLATFORM 1
#endif

// are we building a shared lib?
#define BUILD_SHARED_LIBS

/* #undef OPENMS_BIG_ENDIAN */

#define OPENMS_HASDOXYGENDOT 0

// Define on 64 bit architectures
#define OPENMS_64BIT_ARCHITECTURE

#define PointerSizeInt int64_t
#define PointerSizeUInt uint64_t

#define OPENMS_HAS_UNISTD_H
/* #undef OPENMS_HAS_PROCESS_H */

#define OPENMS_HAS_TIME_H
#define OPENMS_HAS_SYS_TIMES_H
#define OPENMS_HAS_SYS_TIME_H
#define OPENMS_HAS_SYS_RESOURCE_H

#define OPENMS_HAS_KILL
#define OPENMS_HAS_SYSCONF

#define ENABLE_UPDATE_CHECK

#define ENABLE_DOCS

// library versions
#define OPENMS_LIBSVM_VERSION 3.1.2
#define OPENMS_LIBSVM_VERSION_MAJOR 3
#define OPENMS_LIBSVM_VERSION_MINOR 1

#define OPENMS_BOOST_VERSION_MAJOR 1
#define OPENMS_BOOST_VERSION_MINOR 78
#define OPENMS_BOOST_VERSION_SUBMINOR 0
#define OPENMS_BOOST_VERSION 1.78.0

#define OPENMS_HAS_COINOR
#define OPENMS_HAS_COIN_INCLUDE_SUBDIR_IS_COIN

#define OPENMS_GLPK_VERSION 
#define OPENMS_GLPK_VERSION_MAJOR 
#define OPENMS_GLPK_VERSION_MINOR 

// custom compile flags
#define OPENMS_MORECXX_FLAGS "<none>"

// class & TOPP tests

#ifdef _OPENMP
  #define IF_MASTERTHREAD if (omp_get_thread_num() ==0)
#else
  #define IF_MASTERTHREAD
#endif

#define MT_ENABLE_NESTED_OPENMP 1

/* #undef WITH_CRAWDAD */

// NOTE: This is a temporary hack. The aim is that OpenMS should not care
//       but currently we need this information for the ToolHandler.
#define WITH_GUI 1

#endif // OPENMS_CONFIG_H
