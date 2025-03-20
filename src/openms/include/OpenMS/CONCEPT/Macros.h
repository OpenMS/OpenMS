// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once


#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <string>
#include <cstring>

/**
    @brief General helper macros

    @{
*/

#define STRINGIFY(a) #a

#ifdef _OPENMP

// Pragma string literals are compiler-specific: 
// gcc and clang use _Pragma while MSVS uses __pragma
// the MSVS pragma does not need a string token somehow.
#ifdef OPENMS_COMPILER_MSVC
#define OPENMS_THREAD_CRITICAL(name) \
    __pragma(omp critical (name))
#else
#define OPENMS_THREAD_CRITICAL(name) \
    _Pragma( STRINGIFY( omp critical (name) ) )
#endif

#else

#define OPENMS_THREAD_CRITICAL(name) 

#endif

/** @} */ // end of helpers


/**
    @defgroup Conditions Condition macros

    @brief Macros used for to enforce preconditions and postconditions.

    These macros are enabled if debug info is enabled and optimization is disabled in configure.
    Otherwise they are replaced by an empty string, so they won't cost any performance.

    The macros throw Exception::Precondition or Exception::Postcondition respectively if the condition fails.

    @ingroup Concept

    @{
*/

#ifdef OPENMS_ASSERTIONS

/**
    @brief Precondition macro.

    @hideinitializer
*/
#define OPENMS_PRECONDITION(condition, message) \
  if (!(condition)) \
  { \
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, # condition " " # message); \
  } \

/**
    @brief Postcondition macro.

    @hideinitializer
*/
#define OPENMS_POSTCONDITION(condition, message) \
  if (!(condition)) \
  { \
    throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, # condition " " # message); \
  } \

#else

/**
    @brief Precondition macro.

    @hideinitializer
*/
#define OPENMS_PRECONDITION(condition, message)

/**
    @brief Postcondition macro.

    @hideinitializer
*/
#define OPENMS_POSTCONDITION(condition, message)

#endif

/** @} */ // end of Conditions

