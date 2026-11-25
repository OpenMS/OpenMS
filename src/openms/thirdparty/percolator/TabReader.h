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
#ifndef TABREADER_H_
#define TABREADER_H_

#include <string>
#include <string_view>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <charconv>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>

namespace OpenMS { namespace Internal { namespace Percolator {


// Set by src/CMakeLists.txt to whether std::from_chars for double compiles
// and links; libc++ never announces it with __cpp_lib_to_chars.
#ifndef PERCOLATOR_FROM_CHARS_DOUBLE
#if defined(__cpp_lib_to_chars)
#define PERCOLATOR_FROM_CHARS_DOUBLE 1
#else
#define PERCOLATOR_FROM_CHARS_DOUBLE 0
#endif
#endif

// using char pointers is much faster than istringstream
class TabReader {
 public:
  TabReader(const std::string& line)
      : f_(line.c_str()), end_(line.c_str() + line.size()), err(0) {
    errno = 0;
  }
  // reads from begin on in a line that ends with a NUL character at end
  TabReader(const char* begin, const char* end) : f_(begin), end_(end), err(0) {
    errno = 0;
  }

  // the start of the next field
  const char* position() const { return f_; }
  const char* lineEnd() const { return end_; }
  
  void advance(const char* next) {
    if (*next != '\0') {
      f_ = next + 1; // eats up the tab
    } else {
      f_ = next; // prevents pointing over the null byte
    }
  }
  
  void skip() {
    const char* pch = strchr(f_, '\t');
    if (pch == NULL) {
      err = 1;
    } else {
      advance(pch);
    }
  }
  
  void skip(size_t numSkip) {
    for (size_t i = 0; i < numSkip; ++i) skip();
  }
  
  double readDouble() {
#if PERCOLATOR_FROM_CHARS_DOUBLE
    // std::from_chars rounds correctly like strtod but is much faster. Any
    // input it does not fully accept (leading '+' or whitespace, hex,
    // out-of-range values, trailing garbage) takes the strtod path below, so
    // results and error handling are unchanged.
    {
      double d;
      std::from_chars_result res = std::from_chars(f_, end_, d);
      if (res.ec == std::errc() && res.ptr != f_ &&
          (*res.ptr == '\0' || isspace(static_cast<unsigned char>(*res.ptr)))) {
        advance(res.ptr);
        return d;
      }
    }
#else
    // without std::from_chars, most numbers of an input still avoid strtod
    {
      double d;
      const char* next;
      if (readSimpleDouble(f_, end_, d, next)) {
        advance(next);
        return d;
      }
    }
#endif
    char* next = NULL;
    errno = 0;
    double d = strtod(f_, &next);
    if (next == f_ || (*next != '\0' && !isspace(static_cast<unsigned char>(*next)))
                   || ((d == HUGE_VAL || d == -HUGE_VAL) && errno == ERANGE)) {
      err = errno ? errno : 1;
    }
    advance(next);
    return d;
  }
  
  // Clinger's fast path (W. D. Clinger, "How to read floating point numbers
  // accurately", PLDI 1990): if the digits of a decimal number form an
  // integer w <= 2^53 and its decimal exponent q has |q| <= 22, then w and
  // 10^|q| are exact doubles, and w * 10^q or w / 10^-q is a single
  // correctly rounded operation, which gives the value strtod returns.
  // Accepts an optional '-', digits with an optional '.', and an optional
  // exponent, followed by the end of the text or white space as readDouble
  // requires, and sets value and next (the end of the number). Returns false
  // for everything else, which strtod then handles: a '+' sign, leading white
  // space, hex, inf and nan, more than 19 significant digits, larger
  // integers or exponents, and platforms without plain IEEE double
  // arithmetic (x87 excess precision).
  static bool readSimpleDouble(const char* p, const char* end, double& value,
                               const char*& next) {
#if defined(FLT_EVAL_METHOD) && (FLT_EVAL_METHOD == 0 || FLT_EVAL_METHOD == 1)
    static const double kPowersOf10[] = {
        1e0,  1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,  1e8,  1e9,  1e10, 1e11,
        1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
    const bool negative = p < end && *p == '-';
    if (negative) ++p;
    unsigned long long w = 0;  // the significant digits
    int numDigits = 0, exponent = 0;
    bool anyDigit = false;
    for (; p < end && *p >= '0' && *p <= '9'; ++p) {
      anyDigit = true;
      if (w == 0 && *p == '0') continue;  // a leading zero
      if (++numDigits > 19) return false;
      w = 10 * w + static_cast<unsigned>(*p - '0');
    }
    if (p < end && *p == '.') {
      for (++p; p < end && *p >= '0' && *p <= '9'; ++p) {
        anyDigit = true;
        --exponent;
        if (w == 0 && *p == '0') continue;
        if (++numDigits > 19) return false;
        w = 10 * w + static_cast<unsigned>(*p - '0');
      }
    }
    if (!anyDigit) return false;
    if (p < end && (*p == 'e' || *p == 'E')) {
      ++p;
      bool negativeExponent = false;
      if (p < end && (*p == '-' || *p == '+')) {
        negativeExponent = *p == '-';
        ++p;
      }
      if (!(p < end && *p >= '0' && *p <= '9')) return false;
      int e = 0;
      for (; p < end && *p >= '0' && *p <= '9'; ++p) {
        if (e > 100000) return false;
        e = 10 * e + (*p - '0');
      }
      exponent += negativeExponent ? -e : e;
    }
    if (p < end && *p != '\0' && !isspace(static_cast<unsigned char>(*p))) {
      return false;
    }
    double d = 0.0;
    if (w != 0) {
      if (w > (1ULL << 53) || exponent > 22 || exponent < -22) return false;
      d = static_cast<double>(w);
      if (exponent > 0) {
        d *= kPowersOf10[exponent];
      } else if (exponent < 0) {
        d /= kPowersOf10[-exponent];
      }
    }
    value = negative ? -d : d;
    next = p;
    return true;
#else
    (void)p; (void)end; (void)value; (void)next;
    return false;
#endif
  }

  int readInt() {
    char* next = NULL;
    errno=0;
    // long long, not long: with 32-bit long (Windows), strtol saturates at
    // INT_MAX, and the range check below never fires
    long long val = strtoll(f_, &next, 10);
    if (next == f_ || (*next != '\0' && !isspace(static_cast<unsigned char>(*next)))
                   || errno == ERANGE || val < INT_MIN || val > INT_MAX) {
      err = errno ? errno : 1;
    }
    advance(next);
    return static_cast<int>(val);
  }
  
  std::string readString() {
    return std::string(readView());
  }

  // the next field, as a view into the line
  std::string_view readView() {
    const char* pch = strchr(f_, '\t');
    if (pch == NULL) {
      err = 1;
      return std::string_view(f_);
    } else {
      std::string_view s(f_, static_cast<size_t>(pch - f_));
      advance(pch);
      return s;
    }
  }

  bool error() { return err != 0; }
 private:
  const char* f_;
  const char* end_;
  int err;
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /*TABREADER_H_*/