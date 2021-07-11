// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Terminal-related functionality.
// ==========================================================================

#ifndef SEQAN_MISC_MISC_TERMINAL_H_
#define SEQAN_MISC_MISC_TERMINAL_H_

#include <cstdio>

#include <seqan/platform.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Classes, Structs, Enums, Tags
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function isTerminal()
// ----------------------------------------------------------------------------

/**
.Function.isTerminal()
..cat:Miscellaneous
..summary:Check whether we are printing to a terminal.
..signature:isTerminal()
..returns:$true$ if we are on the terminal, false otherwise.
..include:seqan/misc/misc_terminal.h
*/

#if defined(PLATFORM_WINDOWS)

#include <io.h>

inline bool isTerminal()
{
    return false;  // Windows does not understand ANSI codes.
}

#endif  // #if defined(PLATFORM_WINDOWS)

#if defined(PLATFORM_GCC)

#include <unistd.h>

inline bool isTerminal()
{
#ifdef SEQAN_NO_TERMINAL
    return false;    // explicitly disable false positive terminal detection
#else
    return isatty(fileno(stdout));
#endif
}

#endif  // #if defined(PLATFORM_GCC)

// ----------------------------------------------------------------------------
// Function isAnsiColorTerminal()
// ----------------------------------------------------------------------------

/**
.Function.isAnsiColorTerminal
..cat:Miscellaneous
..summary:Check whether we are printing to a terminal.
..signature:isTerminal()
..returns:$true$ if we are in a terminal and the terminal knows ANSI color codes.
..remarks:Currently, we assume that Unix terminals support color while Windows terminals and no terminals do not.
..include:seqan/misc/misc_terminal.h
*/

#if defined(PLATFORM_WINDOWS) || defined(PLATFORM_GCC_MINGW)

inline bool isAnsiColorTerminal()
{
    return false;
}

#else  // #if defined(PLATFORM_WINDOWS) || defined(PLATFORM_GCC_MINGW)

inline bool isAnsiColorTerminal()
{
    return isTerminal();
}

#endif  // #if defined(PLATFORM_WINDOWS) || defined(PLATFORM_GCC_MINGW)

// ----------------------------------------------------------------------------
// Function getTerminalSize()
// ----------------------------------------------------------------------------

/**
.Function.getTerminalSize
..cat:Miscellaneous
..summary:Retrieve size of terminal.
..signature:getTerminalSize(cols, rows)
..param.cols:The column count is written here.
...type:nolink:$unsigned$
..param.rows:The line count is written here.
...type:nolink:$unsigned$
..returns:$true$ on success, $false$ otherwise.
...type:nolink:$bool$
..remarks:On Windows, $rows$ contains the number of rows in the terminal *buffer*, not the window.
..include:seqan/misc/misc_terminal.h
..example:The following demonstrates the usage.
..example.code:
#include <seqan/misc/misc_terminal.h>

unsigned rows = 0, cols = 0;
bool success = getTerminalSize(cols, rows);
std::cerr << "rows == " << rows << ", cols == " << cols
          << ", success == " << success << "\n";
*/

#if defined(PLATFORM_WINDOWS)

#include <Windows.h>

// NOTE: cols actually is the buffer size :(
inline bool getTerminalSize(unsigned & cols, unsigned & rows)
{
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    int ret = GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    if (ret == 0)
        return false;

    rows = csbi.dwSize.X;
    cols = csbi.dwSize.Y;

    return true;
}

#endif  // #if defined(PLATFORM_WINDOWS)

#if defined(PLATFORM_GCC)

#include <sys/ioctl.h>
#include <unistd.h>

inline bool getTerminalSize(unsigned & cols, unsigned & rows)
{
    struct winsize w;
    w.ws_row = 0;
    w.ws_col = 0;

    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

    rows = w.ws_row;
    cols = w.ws_col;

    return true;
}

#endif  // #if defined(PLATFORM_GCC)

}  // namespace seqan

#endif // #ifndef SEQAN_MISC_MISC_TERMINAL_H_
