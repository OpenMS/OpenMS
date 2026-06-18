// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors:  Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

// Small, dependency-light definitions shared by TOPPBase: the Citation struct and
// the TOPPBase-specific parameter exceptions. Split out of TOPPBase.h so they can
// be used without pulling in the full TOPPBase interface.

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>

#include <string>

namespace OpenMS
{
  /**
    @brief Stores Citations for individual TOPP tools.

    An example would be
    \code{.cpp}
      Citation c = {"Pfeuffer J, Bielow C, Wein S, Jeong K, Netz E, Walter A, Alka O et al.",
                    "OpenMS 3 enables reproducible analysis of large-scale mass spectrometry data",
                    "Nat Methods 21, 365–367 (2024)",
                    "10.1038/s41592-024-02197-7"};
    \endcode
    Suggested format is AMA, e.g. https://www.lib.jmu.edu/citation/amaguide.pdf
  */
  struct Citation
  {
    std::string authors;    ///< list of authors in AMA style, i.e. "surname initials", ...
    std::string title;      ///< title of article
    std::string when_where; ///< suggested format: journal. year; volume, issue: pages
    std::string doi;        ///< plain DOI (no urls), e.g. 10.1021/pr100177k

    /// mangle members to string
    std::string toString() const
    {
      return authors + ". " + title + ". " + when_where + ". doi:" + doi + ".";
    }
  };

  namespace Exception
  {
    /// An unregistered parameter was accessed
    class OPENMS_DLLAPI UnregisteredParameter :
      public Exception::BaseException
    {
public:
      UnregisteredParameter(const char* file, int line, const char* function, const std::string& parameter) :
        BaseException(file, line, function, "UnregisteredParameter", parameter)
      {
        GlobalExceptionHandler::getInstance().setMessage(what());
      }

    };
    /// A parameter was accessed with the wrong type
    class OPENMS_DLLAPI WrongParameterType :
      public Exception::BaseException
    {
public:
      WrongParameterType(const char* file, int line, const char* function, const std::string& parameter) :
        BaseException(file, line, function, "WrongParameterType", parameter)
      {
        GlobalExceptionHandler::getInstance().setMessage(what());
      }

    };
    /// A required parameter was not given
    class OPENMS_DLLAPI RequiredParameterNotGiven :
      public Exception::BaseException
    {
public:
      RequiredParameterNotGiven(const char* file, int line, const char* function, const std::string& parameter) :
        BaseException(file, line, function, "RequiredParameterNotGiven", parameter)
      {
        GlobalExceptionHandler::getInstance().setMessage(what());
      }

    };
  }

} // namespace OpenMS
