// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors:  Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/ListUtils.h> // for StringList definition

#include <limits>

namespace OpenMS
{

  /**
  * 
  * Determines the width of the console automatically.
  * 
  * To manually force a certain width set the environment variable 'COLUMNS' to a desired value.
  * 
  */
  class OPENMS_DLLAPI ConsoleUtils
  {
  private:
    /// C'tor (private) -- use ConsoleUtils::getInstance()
    ConsoleUtils();

  public:
    /// Copy C'tor (deleted)
    ConsoleUtils(const ConsoleUtils&) = delete;

    /// Assignment operator (deleted)
    void operator=(ConsoleUtils const&) = delete;

    /// returns the singleton -- the only instanciation of this class
    static const ConsoleUtils& getInstance();

    /// Make a string console-friendly
    /// by breaking it into multiple lines according to the console width.
    /// The 'indentation' gives the number of spaces which is prepended beginning at the second (!)
    /// line, so one gets a left aligned block which has some space to the left.
    /// An indentation of 0 results in the native console's default behaviour: just break at the end of
    /// its width and start a new line.
    /// @p max_lines gives the upper limit of lines returned after breaking is finished.
    /// Excess lines are removed and replaced by '...', BUT the last line will be preserved.
    /// 
    /// @param input String to be split
    /// @param indentation Number of spaces to use for lines 2 until last line (should not exceed the console width)
    /// @param max_lines Limit of output lines (all others are removed)
    /// @param first_line_prefill Assume this many chars were already written in the current line of the console (should not exceed the console width)
    static StringList breakStringList(const String& input, const Size indentation, const Size max_lines, const Size first_line_prefill = 0);

    /// same as breakStringList(), but concatenates the result using '\n' for convenience
    static String breakString(const String& input, const Size indentation, const Size max_lines, const Size first_line_prefill = 0);

    /// width of the console (or INTMAX on internal error)
    int getConsoleWidth() const
    {
      return console_width_;
    }

    friend struct ConsoleWidthTest; ///< allows us to set console_width to a fixed value for testing

  private:
    /// width of console we are currently in (if not determinable, set to INTMAX, i.e. not breaks)
    int console_width_ = std::numeric_limits<int>::max();

    /// read console settings for output shaping
    int readConsoleSize_();

    /// returns a console friendly version of input
    StringList breakString_(const String& input, const Size indentation, const Size max_lines, Size first_line_prefill) const;
  };

} // namespace OpenMS

