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
#include <climits>
#include <cstring>
#include <cmath>

namespace OpenMS { namespace Internal { namespace Percolator {

// using char pointers is much faster than istringstream
class TabReader {
 public:
  TabReader(const std::string& line) : f_(line.c_str()), err(0) {
    errno = 0;
  }
  
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
    char* next = NULL;
    errno = 0;
    double d = strtod(f_, &next);
    if (next == f_ || (*next != '\0' && !isspace(*next))
                   || ((d == HUGE_VAL || d == -HUGE_VAL) && errno == ERANGE)) {
      err = errno ? errno : 1;
    }
    advance(next);
    return d;
  }
  
  int readInt() {
    char* next = NULL;
    errno=0;
    long val = strtol(f_, &next, 10);
    if (next == f_ || (*next != '\0' && !isspace(*next))
                   || val < INT_MIN || val > INT_MAX) {
      err = errno ? errno : 1;
    }
    advance(next);
    return static_cast<int>(val);
  }
  
  std::string readString() {
    const char* pch = strchr(f_, '\t');
    if (pch == NULL) {
      err = 1;
      return std::string(f_);
    } else {
      std::string s(f_, static_cast<std::basic_string<char>::size_type>(pch - f_));
      advance(pch);
      return s;
    }
  }

  bool error() { return err != 0; }
 private:
  const char* f_;
  int err;
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /*TABREADER_H_*/