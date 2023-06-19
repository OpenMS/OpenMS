#ifndef _STOCKHAMSHUFFLE_HPP
#define _STOCKHAMSHUFFLE_HPP

#include <algorithm>
#include "BitReversal.hpp"
#include "RecursiveShuffle.hpp"
#include "../Tensor/Tensor.hpp"

template<typename T, unsigned char LOG_N>
class StockhamShuffle {
public:
  inline static void apply_with_existing_buffer(T* __restrict const v, T* __restrict const buffer) {
    lsb_to_msb_with_existing_buffer<T, LOG_N>(v, buffer);
    StockhamShuffle<T,LOG_N-1>::apply_with_existing_buffer(v, buffer);
    StockhamShuffle<T,LOG_N-1>::apply_with_existing_buffer(v+(1ul<<LOG_N>>1), buffer+(1ul<<LOG_N>>2));
  }
  inline static void apply_out_of_place(T* __restrict const v) {
    T* __restrict const buffer = (T*)aligned_malloc<T>(1ul<<LOG_N>>1);
    apply_with_existing_buffer(v, buffer);
    free(buffer);
  }
};

template<typename T>
class StockhamShuffle<T, 1> {
public:
  inline static void apply_with_existing_buffer(T* __restrict const v, T* __restrict const buffer) {
    // Do nothing
  }
  inline static void apply_out_of_place(T* __restrict const v) {
    // Do nothing
  }
};

#endif
