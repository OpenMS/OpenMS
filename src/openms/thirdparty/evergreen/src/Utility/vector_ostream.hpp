#ifndef _VECTOR_OSTREAM_HPP
#define _VECTOR_OSTREAM_HPP

#include <iostream>
#include <vector>

template <typename T>
std::ostream & operator <<(std::ostream & os, std::vector<T> & rhs) {
  os << "{";
  for (unsigned long i=0; i<rhs.size(); ++i) {
    os << rhs[i];
    if (i+1 != rhs.size()) {
      os << ", ";
    }
  }
  os << "}";
  return os;
}

#endif
