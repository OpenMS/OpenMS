#ifndef _ANY_AND_ALL_HPP
#define _ANY_AND_ALL_HPP

inline bool any(const Vector<bool> & rhs) {
  for (unsigned long k=0; k<rhs.size(); ++k)
    if (rhs[k])
      return true;
  return false;
}

inline bool all(const Vector<bool> & rhs) {
  for (unsigned long k=0; k<rhs.size(); ++k)
    if (! rhs[k])
      return false;
  return true;
}

#endif
