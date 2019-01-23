#ifndef _TABLESHUFFLE_HPP
#define _TABLESHUFFLE_HPP

#include <algorithm>
#include "BitReversal.hpp"

template<typename T, unsigned char LOG_N>
class TableShuffle {
public:
  inline static void apply(T* __restrict const v) {
    constexpr unsigned long N = 1ul << LOG_N;

    for (unsigned long index=1; index<(N-1); ++index) {
      unsigned long reversed = BitReversal<LOG_N>::reverse_bytewise(index);

      // Comparison ensures swap is performed only once per unique pair:
      if (index<reversed)

	std::swap(v[index], v[reversed]);
    }
  }
};

#endif
