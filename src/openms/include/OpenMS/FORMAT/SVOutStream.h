// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <ostream>
#include <fstream>      // std::ofstream
#include <sstream>
#include <boost/math/special_functions/fpclassify.hpp> // because isfinite not supported on Mac

namespace OpenMS
{
  /// custom newline indicator
  enum Newline {nl};

  /**
      @brief Stream class for writing to comma/tab/...-separated values files.

      Automatically inserts separators between items and handles quoting of strings. Requires @p nl (preferred) or @p std::endl as the line delimiter - @p "\n" won't be accepted.

      @ingroup Format
  */
  class OPENMS_DLLAPI SVOutStream :
    public std::ostream
  {
public:

    /**
      @brief Constructor

      @param file_out Output filename; will be overwritten if exists
      @param sep Separator string (typically comma, semicolon, or tab)
      @param replacement If @p quoting is @p NONE, used to replace occurrences of @p sep within strings before writing them
      @param quoting Quoting method for strings (see @p String::quote)
    */
    SVOutStream(const String& file_out,
                const String& sep = "\t",
                const String& replacement = "_",
                String::QuotingMethod quoting = String::DOUBLE);

    /**
      @brief Constructor

      @param out Output stream to write to (open file or @p cout)
      @param sep Separator string (typically comma, semicolon, or tab)
      @param replacement If @p quoting is @p NONE, used to replace occurrences of @p sep within strings before writing them
      @param quoting Quoting method for strings (see @p String::quote)
    */
    SVOutStream(std::ostream& out,
                const String& sep = "\t",
                const String& replacement = "_",
                String::QuotingMethod quoting = String::DOUBLE);

    /** 
      @brief Destructor

      Frees ofstream_* if filename c'tor was used.

    */
    ~SVOutStream() override;

    /**
      @brief Stream output operator for @p String

      The argument is quoted before writing; it must not contain the newline character
    */
    SVOutStream& operator<<(String str);    // use call-by-value here


    /**
      @brief Stream output operator for @p std::string

      The argument is quoted before writing; it must not contain the newline character
    */
    SVOutStream& operator<<(const std::string& str);


    /**
      @brief Stream output operator for @p char*

      The argument is quoted before writing; it must not contain the newline character
    */
    SVOutStream& operator<<(const char* c_str);


    /**
      @brief Stream output operator for @p char

      The argument is quoted before writing; it must not contain the newline character
    */
    SVOutStream& operator<<(const char c);

    /// Stream output operator for manipulators (used to catch @p std::endl)
    SVOutStream& operator<<(std::ostream& (*fp)(std::ostream&));

    /**
      @brief Stream output operator for custom newline (@p nl) without flushing

      Use "nl" instead of "endl" for improved performance
    */
    SVOutStream& operator<<(enum Newline);

    /// numeric types should be converted to String first to make use
    /// of StringConversion
    template<typename T>
    typename std::enable_if<std::is_arithmetic<typename std::remove_reference<T>::type>::value, SVOutStream&>::type operator<<(const T& value)
    {
      if (!newline_) static_cast<std::ostream&>(*this) << sep_;
      else newline_ = false;
      static_cast<std::ostream&>(*this) << String(value);
      return *this;
    };

    /// Generic stream output operator (for non-character-based types)
    template<typename T>
    typename std::enable_if<!std::is_arithmetic<typename std::remove_reference<T>::type>::value, SVOutStream&>::type operator<<(const T& value)
    {
      if (!newline_)
        static_cast<std::ostream &>(*this) << sep_;
      else
        newline_ = false;
      static_cast<std::ostream &>(*this) << value;
      return *this;
    };

    /// Unformatted output (no quoting: useful for comments, but use only on a line of its own!)
    SVOutStream& write(const String& str);   // write unmodified string


    /**
      @brief Switch modification of strings (quoting/replacing of separators) on/off

      @return previous modification state
    */
    bool modifyStrings(bool modify);


    /// Write a numeric value or "nan"/"inf"/"-inf", if applicable (would not be needed for Linux)
    template <typename NumericT>
    SVOutStream& writeValueOrNan(NumericT thing)
    {
      if ((boost::math::isfinite)(thing)) return operator<<(thing);

      bool old = modifyStrings(false);
      if ((boost::math::isnan)(thing)) 
      {
        operator<<(nan_);
      }
      else if (thing < 0) 
      {
        operator<<("-" + inf_);
      }
      else 
      {
        operator<<(inf_);
      }
      modifyStrings(old);
      return *this;
    }

protected:
    /// internal file stream when C'tor is called with a filename
    std::ofstream* ofs_;

    /// Separator string
    String sep_;

    /// Replacement for separator
    String replacement_;

    /// String to use for NaN values
    String nan_;

    /// String to use for Inf values
    String inf_;

    /// String quoting method
    String::QuotingMethod quoting_;

    /// On/off switch for modification of strings
    bool modify_strings_;

    /// Are we at the beginning of a line? (Otherwise, insert separator before next item.)
    bool newline_;

    /// Stream for testing if a manipulator is "std::endl"
    std::stringstream ss_;
  };

}

