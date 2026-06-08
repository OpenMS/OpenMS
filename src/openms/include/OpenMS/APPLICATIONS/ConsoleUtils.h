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
    @brief Singleton helper for terminal-width-aware string wrapping in TOPP-tool console output.

    Detects the current console width on construction (so it reflects the terminal the program
    was launched in) and exposes helpers that wrap long strings to that width with hanging
    indentation — the typical pattern used by TOPP tools to format `--help` output, parameter
    descriptions, and CLI error messages.

    The detected width can be overridden at run time by setting the @c COLUMNS environment
    variable before launching the program. If the width cannot be determined (e.g. when stdout
    is not a terminal), it falls back to @c std::numeric_limits<int>::max(), which effectively
    disables wrapping.

    Construction is private; access via getInstance(). The class is non-copyable and the cached
    console width is read once per process — long-running tools that resize their terminal
    mid-run will continue to use the original width.

    @ingroup System
  */
  class OPENMS_DLLAPI ConsoleUtils
  {
  private:
    /// Private ctor — caches the console width once; use @ref getInstance to access
    ConsoleUtils();

  public:
    /// Copy constructor deleted (singleton)
    ConsoleUtils(const ConsoleUtils&) = delete;

    /// Assignment operator deleted (singleton)
    void operator=(ConsoleUtils const&) = delete;

    /// Return the singleton instance (constructed on first call; thread-safe per Meyers-singleton rules)
    static const ConsoleUtils& getInstance();

    /**
      @brief Wrap @p input to the current console width and return one element per output line.

      Behaviour:
        - Lines 2..N are prefixed by @p indentation spaces, giving a left-aligned hanging
          block; lines wider than the console get broken at the console width.
        - @p indentation == 0 falls back to the native console behaviour (raw line break
          at the right edge with no prefix on the next line).
        - If the wrapped output would exceed @p max_lines, excess lines are dropped and a
          @c "..." marker is inserted; the **final line is always preserved** so the tail of
          the message stays visible.
        - @p first_line_prefill lets the caller advertise how many characters are already
          present on the current console line (e.g. when this string is being appended after
          a label), so the first wrap point is calculated correctly.

      @param[in] input              String to wrap.
      @param[in] indentation        Spaces to prepend to lines 2..N (must not exceed console width).
      @param[in] max_lines          Maximum output lines; excess lines collapse into "..." with the last line kept.
      @param[in] first_line_prefill Number of characters already written on the current line (must not exceed console width).
      @return One element per wrapped output line.
    */
    static StringList breakStringList(const String& input, const Size indentation, const Size max_lines, const Size first_line_prefill = 0);

    /// Convenience wrapper around @ref breakStringList that joins the result with @c "\\n"
    static String breakString(const String& input, const Size indentation, const Size max_lines, const Size first_line_prefill = 0);

    /// Cached console width (or @c std::numeric_limits<int>::max() when detection failed and wrapping is effectively disabled)
    int getConsoleWidth() const
    {
      return console_width_;
    }

    friend struct ConsoleWidthTest; ///< Test hook: lets unit tests pin @c console_width_ to a deterministic value

  private:
    /// Cached console width; @c numeric_limits<int>::max() means "could not detect → no wrapping"
    int console_width_ = std::numeric_limits<int>::max();

    /// Probe the host environment / terminal to determine the console width (also honours the @c COLUMNS env var)
    int readConsoleSize_();

    /// Non-static implementation used by the public static wrappers; receives the cached @c console_width_
    StringList breakString_(const String& input, const Size indentation, const Size max_lines, Size first_line_prefill) const;
  };

} // namespace OpenMS

