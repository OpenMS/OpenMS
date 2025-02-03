// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Moritz Berger $
// $Authors: Chris Bielow, Moritz Berger $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Colorizer.h>
#include <iostream>

#ifdef OPENMS_WINDOWSPLATFORM
  #include <OpenMS/APPLICATIONS/ConsoleUtils.h>
  #include <windows.h>
#endif

#ifdef OPENMS_HAS_UNISTD_H
  #include <unistd.h>  // for isatty(), STDOUT_FILENO, STDERR_FILENO
#endif


namespace OpenMS
{
  /// Initialize the Console (Windows OS) and upon closing the program make sure
  /// that the TTY colors are restored to default again
  struct InitConsole
  {
  
    /// D'tor: reset console colors when exiting the program
    ~InitConsole()
    {
      // use a local Colorizer! The global ones may not exist anymore when this Dtor is called
      Colorizer undo(ConsoleColor::BLUE); // any color will do...
      // reset color to default using both cerr and cout
      // We call both, to ensure the reset is invoked in case either of them is
      // redirected (not a TTY), hence our ANSI code will not reach the stream.
      // See Colorizer::colorStream_()
      std::cout << undo.undoAll();
      std::cerr << undo.undoAll();

      //std::cout << "\nundone coloring\n";
      //std::cerr << "\nundone coloring\n";
    }

#ifdef OPENMS_WINDOWSPLATFORM
    // Windows 10 (since it Anniversary Edition from 2016) understands ANSI codes, but only
    // if we tell it to ...
    InitConsole()
    {
      initStream(STD_OUTPUT_HANDLE);
      initStream(STD_ERROR_HANDLE);
    }

    // Set output mode to handle virtual terminal sequences
    DWORD initStream(DWORD handle)
    {
      
      HANDLE hOut = GetStdHandle(handle);
      if (hOut == INVALID_HANDLE_VALUE)
      {
        //std::cerr << "no " << handle << "\n";
        return GetLastError();
      }

      DWORD dwMode = 0;
      if (!GetConsoleMode(hOut, &dwMode))
      {
        //std::cerr << "no mode get for " << handle << "\n";
        return GetLastError();
      }

      dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
      if (!SetConsoleMode(hOut, dwMode))
      {
        //std::cerr << "no mode set for " << handle << "\n";
        return GetLastError();
      }
      return dwMode;
    }
#endif // end OPENMS_WINDOWSPLATFORM
  };

#ifdef OPENMS_WINDOWSPLATFORM
  /// implementation of isatty() for Windows
  /// Returns 'true' if the stream is shown on the console
  /// and 'false' if redirected somewhere else (file or NUL)
  bool isattyWin(const std::ostream& stream)
  {
    DWORD h_stream;
    if (&stream == &std::cout)
      h_stream = STD_OUTPUT_HANDLE;
    else if (&stream == &std::cerr)
      h_stream = STD_ERROR_HANDLE;
    else
      return false;

    HANDLE hOut = GetStdHandle(h_stream);
    if (hOut == INVALID_HANDLE_VALUE)
    {
      //std::cerr << "no handle for " << h_stream << "\n";
      return false;
    }

    DWORD dwMode = 0;
    if (!GetConsoleMode(hOut, &dwMode))
    {
      return false;
    }
    return true;
  }

#endif

  // our local object, which will be initialized and destroyed at startup/teardown
  static InitConsole windows_console_prep_and_restore;


  Colorizer::Colorizer(const ConsoleColor color) 
    : color_(color)
  {
  }

  Colorizer& Colorizer::undo()
  {
    this->input_.str("");
    undo_ = true;
    undos_only = true;
    return *this;
  }

  Colorizer& Colorizer::undoAll()
  {
    this->input_.str("");
    undo_all_ = true;
    undos_only = true;
    return *this;
  }

  void Colorizer::colorStream_(std::ostream& stream, const char* ANSI_command)
  {
    // check if the output is being fed to file or console
    // supress output of ANSI codes into a redirected cout/cerr file
    if (&stream == &std::cout || &stream == &std::cerr)
    {
      if (!isTTY(stream))
      { 
        return;
      }
    }
    // color cout/cerr if visible, or any other stream (mostly for testing purposes)
    // debug: stream << "(" << ANSI_command + 2 << ") ";
    stream << ANSI_command;
  }

  bool Colorizer::isTTY(const std::ostream& stream)
  {
  #ifdef OPENMS_WINDOWSPLATFORM
    return isattyWin(stream);
  #else
    if (&stream == &std::cout && isatty(STDOUT_FILENO))
    {
      return true;
    }
    if (&stream == &std::cerr && isatty(STDERR_FILENO))
    {
      return true;
    }
    return false;
  #endif
  }

  void Colorizer::outputToStream_(std::ostream& o_stream)
  {
    if (!undos_only) // undo() or undoAll() were called - do not color the stream or print empty data; it does not make sense
    {
      // color the stream (or console)
      colorStream_(o_stream, colors_[(int)color_].enable);

      // paste text
      o_stream << input_.str();
    }

    if (undo_all_)
    {
      colorStream_(o_stream, color_undo_all_);
    }
    else if (undo_)
    {
      colorStream_(o_stream, colors_[(int)color_].disable);
    }
  }

  // extern colizers with predefined colors
  OpenMS::Colorizer red(ConsoleColor::RED);
  OpenMS::Colorizer green(ConsoleColor::GREEN);
  OpenMS::Colorizer yellow(ConsoleColor::YELLOW);
  OpenMS::Colorizer blue(ConsoleColor::BLUE);
  OpenMS::Colorizer magenta(ConsoleColor::MAGENTA);
  OpenMS::Colorizer cyan(ConsoleColor::CYAN);
  OpenMS::Colorizer invert(ConsoleColor::INVERT);
  OpenMS::Colorizer bright(ConsoleColor::BRIGHT);
  OpenMS::Colorizer underline(ConsoleColor::UNDERLINE);

  std::ostream& operator<<(std::ostream& o_stream, OpenMS::Colorizer& col)
  {
    // colorize stream; dump internal string (if any); and reset the color (if col.resetColor() was used before).
    col.outputToStream_(o_stream);
    return o_stream;
  }

} // namespace OpenMS
