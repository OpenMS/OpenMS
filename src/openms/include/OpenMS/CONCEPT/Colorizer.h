// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Moritz Berger $
// $Authors: Chris Bielow, Moritz Berger $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <array>
#include <iosfwd>
#include <sstream>

namespace OpenMS
{
  /// Text colors/styles supported by Colorizer
  enum class ConsoleColor
  { // note: the order is important here! See Colorizer::colors_
    RED,
    GREEN,
    YELLOW,
    BLUE,
    MAGENTA,
    CYAN,
    UNDERLINE,
    BRIGHT, ///< keeps the foreground color, but makes it brighter
    INVERT, ///< invert foreground and background color (inverting twice stays inverted)
  };

  /**
    @brief Color and style the fonts shown on cout/cerr (or other streams)
   
    Allows to color the console fonts' foreground color (making the current color brighter, or setting a new color)
    and/or add an underline.
    There are predefined Colorizer objects for your convenience, e.g. 'red' or 'underline'.
    They are named identically to the possible values of enum ConsoleColor.
    Multiple styles can be combined (with limitations for nesting colors, see below).
    To undo a color/style, simply call its '.undo()' function.
    To undo all modifications and return to a normal console color/style, call 'undoAll()'.

    e.g.
    \code
      std::cout << "normal text" << red() << "this is red" << underline() << "red and underlined" << red.undo() << "just underlined" << underline.undo();
    \endcode

    You can also modify a single item and immediately return to the previous setting by passing the item to the
    bracket operator, e.g.
    \code
      std::cout << "normal text << red("this part is red") << " back to normal here";
    \endcode
    
    Undo() has limitations and does <b>not support nesting of colors</b>, i.e. does not remember the previous color, but resets to default
    \code
      std::cout << red() << "this is red" << "blue() << "this is blue" << blue.undo() << " not red! but normal";
    \endcode

    <b>Redirecting cout/cerr streams</b><br>
    If std::cout or std::cerr are redirected to a file, then Colorizer will detect this and <b>not emit</b>
    any color/style information. This is to avoid the ANSI color codes showing up in the file instead of being
    filtered out by the console's color process.

    Note: all OS's we know (Windows, Linux, MacOS) only have a single color configuration for the whole
          console/terminal, independent of streams. I.e. if you apply a permanent color to std::cout
          then all subsequent output to std::cerr will also be colored! (unless std::cout is redirected to a file,
          then coloring std::cout will have no effect, not even when printing to the console using std::cerr).
   */
  class OPENMS_DLLAPI Colorizer
  {
  public:
    /// stream insertion, e.g. std::cout << red() << "this is red" << red.undo();
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& o_stream, Colorizer& col);
    friend class IndentedStream; // to manipulate the internal string data before writing to a stream

    /// Constructor
    Colorizer(const ConsoleColor color);

    /// Constructor (deleted)
    Colorizer() = delete;

    /// Copy constructor (deleted)
    //. Explicitly deleted here, since stringstream member has no copy c'tor either
    /// If you really wanted to, you can implement a copy c'tor, but there seems little use
    Colorizer(const Colorizer& rhs) = delete;

    /// Destructor
    ~Colorizer() = default;

    /// bracket operator to prepare this colorizer object to
    /// color the next stream it is applied to with this object's color 
    /// until it is reset somewhere downstream.
    /// e.g. 'cerr << red() << "this is red" << " this too ..." << red.undo();'
    Colorizer& operator()()
    {
      input_.str(""); // clear the internal buffer
      undo_ = false;
      undo_all_ = false;
      undos_only = false;
      return *this;
    }

    /// Consume @p s, convert it to a string and insert this string with color into the next
    /// stream this Colorizer is applied to. Reset the stream right away.
    /// e.g. 'cerr << red("make this red!") << " this is not red ";
    template<typename T>
    Colorizer& operator()(T s)
    {
      input_.str(""); // clear the internal buffer
      input_ << s;    // add new data
      undo_ = true;
      undo_all_ = false;
      undos_only = false;
      return *this;
    }

    /// prepare this colorizer to undo the effect of this Colorizer object to the next stream it is applied to
    /// e.g. 'cerr << red.undo() << "not red anymore"'
    Colorizer& undo();

    /// prepare this colorizer to reset the console to its default colors/style.
    Colorizer& undoAll();

    /// A wrapper for POSIX 'isatty()' or the equivalent on Windows.
    /// Checks if the stream is written to the terminal/console
    /// or being redirected to a file/NULL/nul (i.e. not visible).
    /// This only works for std::cout and std::cerr. Passing any other stream
    /// will always return 'false'.
    static bool isTTY(const std::ostream& stream);

  protected:
    /// write the content of input_ to the stream
    void outputToStream_(std::ostream& o_stream);

    /// color the @p stream using the given color
    static void colorStream_(std::ostream& stream, const char* ANSI_command);

    /// color of the stream (const; set in C'tor)
    const ConsoleColor color_;

    /// clear the color/style of this Colorizer from a stream (usually cout/cerr) upon the next call of
    /// 'std::ostream& operator<<(std::ostream& o_stream, OpenMS::Colorizer& col)'
    bool undo_ = true;

    /// clear all color/styles from a stream (usually cout/cerr) upon the next call of
    /// 'std::ostream& operator<<(std::ostream& o_stream, OpenMS::Colorizer& col)'
    bool undo_all_ = true;

    /// optimization to prevent coloring a stream, print nothing and then immediately
    /// uncoloring is, e.g. when calling c.undo() or c.undoAll()
    bool undos_only = false;

    /// internal methods used by friends to manipulate the internal data (if present)
    /// This avoids exposing the stream coloring itself (which is easy to get wrong)
    const std::stringstream& getInternalChars_() const
    {
      return input_;
    }


    /// internal methods used by friends to manipulate the internal data (if present)
    /// This avoids exposing the stream coloring itself (which is easy to get wrong)
    void setInternalChars_(const std::string& data)
    {
      input_.str(data);
    }


    /// holds ANSI codes to switch to a certain color or undo this effect
    /// The undo may not be a perfect fit when nesting,
    /// e.g. 'cout << red() << blue() << blue.undo()' will not yield red, but the default color
    struct ColorWithUndo_
    {
      const char* enable; ///< ANSI code to activate the color/style
      const char* disable; ///< ANSO code to undo the color/style
    };
    
    /**
     * @brief ANSI colors/styles, corresponding to values of enum ConsoleColor
     */
    inline static const constexpr std::array<ColorWithUndo_, 9> colors_ {
      ColorWithUndo_ {"\033[91m", "\033[39m"}, // red
      ColorWithUndo_ {"\033[92m", "\033[39m"}, // green
      ColorWithUndo_ {"\033[93m", "\033[39m"}, // yellow
      ColorWithUndo_ {"\033[94m", "\033[39m"}, // blue
      ColorWithUndo_ {"\033[95m", "\033[39m"}, // magenta
      ColorWithUndo_ {"\033[96m", "\033[39m"}, // cyan
      ColorWithUndo_ {"\033[4m", "\033[24m"},  // underline
      ColorWithUndo_ {"\033[1m", "\033[22m"},  // bright/intensified
      ColorWithUndo_ {"\033[7m",  "\033[27m"}, // invert FG/BG
    };
    const char* color_undo_all_ = "\033[0m"; ///< resets all attributes to default

    /// internal string buffer when using operator()(T data)
    /// This data will be colored
    std::stringstream input_;

  }; // class Colorizer

  // declaration of all global colorizer objects
  extern OPENMS_DLLAPI Colorizer red;
  extern OPENMS_DLLAPI Colorizer green;
  extern OPENMS_DLLAPI Colorizer yellow;
  extern OPENMS_DLLAPI Colorizer blue;
  extern OPENMS_DLLAPI Colorizer magenta;
  extern OPENMS_DLLAPI Colorizer cyan;
  extern OPENMS_DLLAPI Colorizer bright;
  extern OPENMS_DLLAPI Colorizer underline;
  extern OPENMS_DLLAPI Colorizer invert;

  /// stream operator for Colorizers, which will (depending on the state of @p col)
  /// add color to the @p o_stream, and print some internal string, 
  /// e.g. 'cout << red("this is red")' or 'cout << red() << "stream stays red until dooms day";'
  /// or reset the color
  /// e.g. 'cout << red.undo() << "not red anymore"'
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& o_stream, OpenMS::Colorizer& col);
} // namespace OpenMS
