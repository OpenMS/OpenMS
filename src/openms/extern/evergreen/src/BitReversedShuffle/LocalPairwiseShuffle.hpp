#ifndef _LOCALPAIRWISESHUFFLE_HPP
#define _LOCALPAIRWISESHUFFLE_HPP

#include <algorithm>

// From Jośe M. Ṕerez-Jord́a1 1997
template<typename T, unsigned char LOG_N, unsigned char LOG_SUB_N>
class LocalPairwiseShuffleHelper {
public:
  inline static void apply(T* __restrict const v) {
    constexpr unsigned long SUB_N = 1ul << LOG_SUB_N;
    constexpr unsigned long RECURSION_DEPTH = LOG_N-LOG_SUB_N;

    //    RECURSION_DEPTH=0: skip 1 (start at 1), += 2 , do blocks of 1
    //    RECURSION_DEPTH=1: skip 2 (start at 2), += 4 , do blocks of 2
    // ...

    // Find indices with bitstrings ending with 1 (end defined by
    // LOG_SUB_N). Swap them with the matching reversed index.
    for (unsigned long index=1ul<<RECURSION_DEPTH; index<(SUB_N>>1); index+=(1ul<<RECURSION_DEPTH)) {
      for (unsigned long block=0; block<1ul<<RECURSION_DEPTH; ++block, ++index) {
	unsigned long pair_bit_reversed = (index & ~(1ul<<RECURSION_DEPTH)) | (1ul<<(LOG_SUB_N-1));
	std::swap(v[index], v[pair_bit_reversed]);
      }
    }
      
    // Recursively apply to first and second half of the list:
    LocalPairwiseShuffleHelper<T, LOG_N, LOG_SUB_N-1>::apply(v);
    LocalPairwiseShuffleHelper<T, LOG_N, LOG_SUB_N-1>::apply(v + (1ul<<(LOG_SUB_N-1)) );
  }
};

template<typename T, unsigned char LOG_N>
class LocalPairwiseShuffleHelper<T, LOG_N, 0> {
public:
  inline static void apply(T* __restrict const v) {
    // Do nothing.
  }
};

template<typename T, unsigned char LOG_N>
class LocalPairwiseShuffle {
public:
  inline static void apply(T* __restrict const v) {
    LocalPairwiseShuffleHelper<T, LOG_N, LOG_N>::apply(v);
  }
};

#endif
