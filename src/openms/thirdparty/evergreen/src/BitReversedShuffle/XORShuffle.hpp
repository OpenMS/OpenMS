#ifndef _XORSHUFFLE_HPP
#define _XORSHUFFLE_HPP

#include <algorithm>
#include "BitReversal.hpp"

template<typename T, unsigned char LOG_N>
class XORShuffle {
public:
  inline static void apply(T* __restrict const v) {
    constexpr unsigned long N = 1ul << LOG_N;

    unsigned long reversed = 0;
    for (unsigned long index=0; index<(N-1); ) {
      // Comparison ensures swap is performed only once per unique pair:
      if (index<reversed)
	std::swap(v[index], v[reversed]);

      BitReversal<LOG_N>::advance_index_and_reversed(index, reversed);
    }
  }
};

#endif
