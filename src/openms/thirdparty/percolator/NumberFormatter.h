/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
#ifndef NUMBER_FORMATTER_H_
#define NUMBER_FORMATTER_H_

#include <charconv>
#include <cmath>
#include <cstdio>
#include <ios>
#include <locale>
#include <string>

#include "Globals.h"

namespace OpenMS { namespace Internal { namespace Percolator {


// Set by src/CMakeLists.txt to whether std::to_chars for double compiles
// and links; libc++ provides it without announcing it.
#ifndef PERCOLATOR_TO_CHARS_DOUBLE
#if defined(__cpp_lib_to_chars)
#define PERCOLATOR_TO_CHARS_DOUBLE 1
#else
#define PERCOLATOR_TO_CHARS_DOUBLE 0
#endif
#endif

/*
 * Appends numbers to a std::string with the same characters as operator<<
 * on a given stream, several times faster. The output code builds its lines
 * in a buffer with it and writes them to the stream in large blocks.
 *
 * operator<< writes a double as printf("%.*g") with the stream's precision
 * in the "C" locale, and then applies the stream's locale. std::to_chars
 * with chars_format::general and a precision is specified to give the
 * characters of printf, so the two agree if the stream has the classic
 * locale, the default format flags and no field width. Otherwise, or if
 * std::to_chars for double is missing, sameAsStream() is false and the
 * caller writes to the stream with operator<< as before. Infinities and
 * NaNs are written with snprintf, as operator<< does, because the C
 * libraries spell NaNs differently.
 */
class NumberFormatter {
 public:
  // output lines are buffered up to about this many bytes per write
  static const size_t kWriteBlockSize = 1u << 16;

  explicit NumberFormatter(const std::ios_base& stream)
      : precision_(static_cast<int>(stream.precision())) {
    const std::ios_base::fmtflags flags =
        stream.flags() &
        (std::ios_base::floatfield | std::ios_base::basefield |
         std::ios_base::showpoint | std::ios_base::showpos |
         std::ios_base::showbase | std::ios_base::uppercase);
    sameAsStream_ = PERCOLATOR_TO_CHARS_DOUBLE &&
                    (flags == 0 || flags == std::ios_base::dec) &&
                    stream.width() == 0 && precision_ >= 1 &&
                    precision_ <= 17 &&
                    stream.getloc() == std::locale::classic();
  }

  bool sameAsStream() const { return sameAsStream_; }

  void append(std::string& out, double value) const {
    char buf[64];
    int length;
#if PERCOLATOR_TO_CHARS_DOUBLE
    if (std::isfinite(value)) {
      length = static_cast<int>(
          std::to_chars(buf, buf + sizeof(buf), value,
                        std::chars_format::general, precision_).ptr - buf);
    } else
#endif
    {
      length = std::snprintf(buf, sizeof(buf), "%.*g", precision_, value);
    }
    out.append(buf, static_cast<size_t>(length));
  }

  template <class Integer>
  static void appendInteger(std::string& out, Integer value) {
    char buf[24];
    out.append(buf, static_cast<size_t>(
                        std::to_chars(buf, buf + sizeof(buf), value).ptr - buf));
  }

 private:
  int precision_;
  bool sameAsStream_;
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /* NUMBER_FORMATTER_H_ */
