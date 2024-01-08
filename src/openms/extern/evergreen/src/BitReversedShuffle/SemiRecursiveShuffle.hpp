#ifndef _SEMIRECURSIVESHUFFLE_HPP
#define _SEMIRECURSIVESHUFFLE_HPP

#include "RecursiveShuffle.hpp"

// Identical to RecursiveShuffle but limits the number of recursions
// to 1. Therefore, it will compile longer, but can be slightly faster
// (empirically tested, reason unclear).
template <typename T, unsigned char NUM_BITS, unsigned char RECURSIONS_REMAINING>
class SemiRecursiveShuffle {
public:
  inline static void apply(T* __restrict x) {
    // & 1 is the same as % 2:
    if ((NUM_BITS & 1) == 1) {
      // allocate buffer and perform single LSB --> MSB:
      lsb_to_msb<T, NUM_BITS>(x);

      SemiRecursiveShuffle<T, NUM_BITS-1, RECURSIONS_REMAINING>::apply(x);
      SemiRecursiveShuffle<T, NUM_BITS-1, RECURSIONS_REMAINING>::apply(x+(1ul<<(NUM_BITS-1)));
    }
    else {
      constexpr unsigned char SUB_NUM_BITS = NUM_BITS>>1;
      constexpr unsigned long SUB_N = 1ul<<SUB_NUM_BITS;
      
      for (unsigned long k=0; k<SUB_N; ++k)
	SemiRecursiveShuffle<T, SUB_NUM_BITS, RECURSIONS_REMAINING-1>::apply(x+(k<<SUB_NUM_BITS));

      MatrixTranspose<T>::apply_square(x, SUB_N);
      
      for (unsigned long k=0; k<SUB_N; ++k)
	SemiRecursiveShuffle<T, SUB_NUM_BITS, RECURSIONS_REMAINING-1>::apply(x+(k<<SUB_NUM_BITS));
    }
  }
};

template <typename T, unsigned char NUM_BITS>
class SemiRecursiveShuffle<T, NUM_BITS, 0> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, NUM_BITS>::apply(x);
  }
};

// NUM_BITS <= 9, simply use UnrolledShuffle:
template <typename T, unsigned char RECURSIONS_REMAINING>
class SemiRecursiveShuffle<T, 9, RECURSIONS_REMAINING> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 9>::apply(x);
  }
};

template <typename T, unsigned char RECURSIONS_REMAINING>
class SemiRecursiveShuffle<T, 8, RECURSIONS_REMAINING> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 8>::apply(x);
  }
};

template <typename T, unsigned char RECURSIONS_REMAINING>
class SemiRecursiveShuffle<T, 7, RECURSIONS_REMAINING> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 7>::apply(x);
  }
};

template <typename T, unsigned char RECURSIONS_REMAINING>
class SemiRecursiveShuffle<T, 6, RECURSIONS_REMAINING> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 6>::apply(x);
  }
};

template <typename T, unsigned char RECURSIONS_REMAINING>
class SemiRecursiveShuffle<T, 5, RECURSIONS_REMAINING> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 5>::apply(x);
  }
};

template <typename T, unsigned char RECURSIONS_REMAINING>
class SemiRecursiveShuffle<T, 4, RECURSIONS_REMAINING> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 4>::apply(x);
  }
};

template <typename T, unsigned char RECURSIONS_REMAINING>
class SemiRecursiveShuffle<T, 3, RECURSIONS_REMAINING> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 3>::apply(x);
  }
};

template <typename T, unsigned char RECURSIONS_REMAINING>
class SemiRecursiveShuffle<T, 2, RECURSIONS_REMAINING> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 2>::apply(x);
  }
};

template <typename T, unsigned char RECURSIONS_REMAINING>
class SemiRecursiveShuffle<T, 1, RECURSIONS_REMAINING> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 1>::apply(x);
  }
};

template <typename T, unsigned char RECURSIONS_REMAINING>
class SemiRecursiveShuffle<T, 0, RECURSIONS_REMAINING> {
public:
  inline static void apply(T* __restrict x) {
  }
};

// Special cases to prevent ambiguity in template specialization:
template <typename T>
class SemiRecursiveShuffle<T, 9, 0> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 9>::apply(x);
  }
};

template <typename T>
class SemiRecursiveShuffle<T, 8, 0> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 8>::apply(x);
  }
};

template <typename T>
class SemiRecursiveShuffle<T, 7, 0> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 7>::apply(x);
  }
};

template <typename T>
class SemiRecursiveShuffle<T, 6, 0> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 6>::apply(x);
  }
};

template <typename T>
class SemiRecursiveShuffle<T, 5, 0> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 5>::apply(x);
  }
};

template <typename T>
class SemiRecursiveShuffle<T, 4, 0> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 4>::apply(x);
  }
};

template <typename T>
class SemiRecursiveShuffle<T, 3, 0> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 3>::apply(x);
  }
};

template <typename T>
class SemiRecursiveShuffle<T, 2, 0> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 2>::apply(x);
  }
};

template <typename T>
class SemiRecursiveShuffle<T, 1, 0> {
public:
  inline static void apply(T* __restrict x) {
    UnrolledShuffle<T, 1>::apply(x);
  }
};

template <typename T>
class SemiRecursiveShuffle<T, 0, 0> {
public:
  inline static void apply(T* __restrict x) {
  }
};

#endif
