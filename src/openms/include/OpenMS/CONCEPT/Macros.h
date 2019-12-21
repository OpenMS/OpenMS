// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    Exception::Precondition e(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, # condition); \
    if (std::strcmp(message, "") != 0) \
    { \
      ::std::string tmp(e.getMessage()); \
      tmp += " "; \
      tmp += ::std::string(message); \
      e.setMessage(tmp); \
    } \
    throw e; \
  } \

/**
    @brief Postcondition macro.

    @hideinitializer
*/
#define OPENMS_POSTCONDITION(condition, message) \
  if (!(condition)) \
  { \
    Exception::Postcondition e(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, # condition); \
    if (std::strcmp(message, "") != 0) \
    { \
      std::string tmp(e.getMessage()); \
      tmp += " "; \
      tmp += std::string(message); \
      e.setMessage(tmp); \
    } \
    throw e; \
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

