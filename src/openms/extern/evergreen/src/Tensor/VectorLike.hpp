#ifndef _VECTORLIKE_HPP
#define _VECTORLIKE_HPP

template <typename T>
class VectorView;

template <typename T>
class WritableVectorView;

// Never instantiate these; always pass by reference & or &&:

template <typename T, template <typename> class VECTOR >
class VectorLike {
public:
  unsigned long size() const {
    return static_cast<const VECTOR<T> &>(*this).size();
  }
  const T & operator [](unsigned long i) const {
    return static_cast<const VECTOR<T> &>(*this)[i];
  }
  operator const T*const() const {
    return (const T*const) static_cast<const VECTOR<T> &>(*this);
  }

  VectorView<T> start_at_const(unsigned long start) const {
    return static_cast<const VECTOR<T> &>(*this).start_at_const(start);
  }
  VectorView<T> start_at_const(unsigned long start, unsigned long length) const {
    return static_cast<const VECTOR<T> &>(*this).start_at_const(start, length);
  }
};

template <typename T, template <typename> class VECTOR >
class WritableVectorLike : public VectorLike<T, VECTOR> {
public:
  T & operator [](unsigned long i) {
    return static_cast<VECTOR<T> &>(*this)[i];
  }
  void fill(T val) {
    for (unsigned long k=0; k<this->size(); ++k)
      (*this)[k] = val;
  }
  operator T*const() const {
    return (T*const)static_cast<const VECTOR<T> &>(*this);
  }
  WritableVectorView<T> start_at(unsigned long start) const {
    return static_cast<const VECTOR<T> &>(*this).start_at(start);
  }
  WritableVectorView<T> start_at(unsigned long start, unsigned long length) const {
    return static_cast<const VECTOR<T> &>(*this).start_at(start, length);
  }
};

template <typename T, typename S, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
// rvalue reference so that it can accept temporary values; however,
// the function is non-destructive:
void copy(WritableVectorLike<T, VECTOR_A> & lhs, const VectorLike<S, VECTOR_B> & rhs) {
  #ifdef SHAPE_CHECK
  assert(lhs.size() >= rhs.size());
  #endif
  for (unsigned int k=0; k<rhs.size(); ++k)
    lhs[k] = (T)rhs[k];
}

template <typename T, typename S, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
// rvalue reference so that it can accept temporary values; however,
// the function is non-destructive:
void copy(WritableVectorLike<T, VECTOR_A> && lhs, const VectorLike<S, VECTOR_B> & rhs) {
  // Calls & version above:
  copy(lhs, rhs);
}

template <typename T, template <typename> class VECTOR>
std::ostream & operator <<(std::ostream & os, const VectorLike<T, VECTOR> & rhs) {
  os << "[";
  for (unsigned long k=0; k<rhs.size(); ++k) {
    os << rhs[k];
    if (k != rhs.size()-1)
      os << ", ";
  }
  os << "]";
  return os;
}

#endif
