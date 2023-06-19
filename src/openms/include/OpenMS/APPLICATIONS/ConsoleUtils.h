// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

