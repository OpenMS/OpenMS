#ifndef _SIGN_HPP
#define _SIGN_HPP

template <typename T>
int sign(T val) {
  return (val > T(0)) - (val < T(0));
}

#endif
