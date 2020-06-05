#ifndef _FROM_STRING_HPP
#define _FROM_STRING_HPP

#include <sstream>
#include <string>

double from_string(const std::string & s) {
  std::istringstream ist(s);
  double result;
  ist >> result;
  return result;
}

#endif
