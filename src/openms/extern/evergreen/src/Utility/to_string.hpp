#ifndef _TO_STRING_HPP
#define _TO_STRING_HPP

#include <string>
#include <sstream>

template <typename T>
std::string to_string(const T & rhs) {
  std::string res;
  std::ostringstream ost(res);
  ost << rhs;
  return ost.str();
}

#endif
