#ifndef _ARRAYSHAPE_HPP
#define _ARRAYSHAPE_HPP

// Used to expand the shape of an array into a variadic template pack
// at compile time. This particular version puts the shape into a
// Vector.

template <typename T, typename ARR, unsigned long ...SHAPE>
struct ArrayShape {
  static Vector<unsigned long> eval(ARR arg) {
    return ArrayShape<T,decltype(arg[0]), SHAPE..., sizeof(arg) / sizeof(arg[0])>::eval(arg[0]);
  }
};

template <typename T, unsigned long ...SHAPE>
struct ArrayShape<T,T, SHAPE...> {
  static Vector<unsigned long> eval(T element) {
    return Vector<unsigned long>({SHAPE...});
  }
};

#endif
