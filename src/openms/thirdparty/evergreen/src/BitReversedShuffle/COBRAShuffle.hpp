#ifndef _COBRASHUFFLE_HPP
#define _COBRASHUFFLE_HPP

#include <algorithm>
#include "BitReversal.hpp"
#include "../Tensor/Tensor.hpp"

// From Carter and Gatlin 1998
template<typename T, unsigned char LOG_N, unsigned char LOG_BLOCK_WIDTH>
class COBRAShuffle {
public:
  inline static void apply(T* __restrict const v) {
    constexpr unsigned char NUM_B_BITS = LOG_N - 2*LOG_BLOCK_WIDTH;
    constexpr unsigned long B_SIZE = 1ul << NUM_B_BITS;
    constexpr unsigned long BLOCK_WIDTH = 1ul << LOG_BLOCK_WIDTH;

    T* __restrict buffer = aligned_malloc<T>(BLOCK_WIDTH*BLOCK_WIDTH);

    for (unsigned long b=0; b<B_SIZE; ++b) {
      unsigned long b_rev = BitReversal<NUM_B_BITS>::reverse_bytewise(b);
      
      // Copy block to buffer:
      for (unsigned long a=0; a<BLOCK_WIDTH; ++a) {
	unsigned long a_rev = BitReversal<LOG_BLOCK_WIDTH>::reverse_bytewise(a);
	
	for (unsigned long c=0; c<BLOCK_WIDTH; ++c)
	  buffer[ (a_rev << LOG_BLOCK_WIDTH) | c ] = v[ (a << NUM_B_BITS << LOG_BLOCK_WIDTH) | (b << LOG_BLOCK_WIDTH) | c ];
      }

      // Swap v[rev_index] with buffer:
      for (unsigned long c=0; c<BLOCK_WIDTH; ++c) {
      	// Note: Typo in original pseudocode by Carter and Gatlin at
      	// the following line:
      	unsigned long c_rev = BitReversal<LOG_BLOCK_WIDTH>::reverse_bytewise(c);
	
      	for (unsigned long a_rev=0; a_rev<BLOCK_WIDTH; ++a_rev) {
      	  unsigned long a = BitReversal<LOG_BLOCK_WIDTH>::reverse_bytewise(a_rev);
      	  // To guarantee each value is swapped only one time:
      	  // index < reversed_index <-->
      	  // a b c < c' b' a' <-->
      	  // a < c' ||
      	  // a <= c' && b < b' ||
      	  // a <= c' && b <= b' && a' < c

	  bool index_less_than_reverse = a < c_rev || (a == c_rev && b < b_rev) || (a == c_rev && b == b_rev && a_rev < c);
	  if ( index_less_than_reverse )
	    std::swap( v[(c_rev << NUM_B_BITS << LOG_BLOCK_WIDTH) | (b_rev<<LOG_BLOCK_WIDTH) | a_rev], buffer[ (a_rev<<LOG_BLOCK_WIDTH) | c ]);
	}
      }

      // Copy changes that were swapped into buffer above:
      for (unsigned long a=0; a<BLOCK_WIDTH; ++a) {
	unsigned long a_rev = BitReversal<LOG_BLOCK_WIDTH>::reverse_bytewise(a);
	for (unsigned long c=0; c<BLOCK_WIDTH; ++c) {
	  unsigned long c_rev = BitReversal<LOG_BLOCK_WIDTH>::reverse_bytewise(c);
	  bool index_less_than_reverse = a < c_rev || (a == c_rev && b < b_rev) || (a == c_rev && b == b_rev && a_rev < c);
	  
	  if (index_less_than_reverse)
	    std::swap(v[ (a << NUM_B_BITS << LOG_BLOCK_WIDTH) | (b << LOG_BLOCK_WIDTH) | c ], buffer[ (a_rev << LOG_BLOCK_WIDTH) | c ]);
	}
      }
    }
    free(buffer);
  }

  inline static void apply_out_of_place(T* __restrict const v) {
    T* __restrict const result = (T*)aligned_malloc<T>(1ul<<LOG_N);
    
    constexpr unsigned char NUM_B_BITS = LOG_N - 2*LOG_BLOCK_WIDTH;
    constexpr unsigned long B_SIZE = 1ul << NUM_B_BITS;
    constexpr unsigned long BLOCK_WIDTH = 1ul << LOG_BLOCK_WIDTH;

    T* __restrict buffer = (T*)aligned_malloc<T>(BLOCK_WIDTH*BLOCK_WIDTH);

    for (unsigned long b=0; b<B_SIZE; ++b) {
      unsigned long b_rev = BitReversal<NUM_B_BITS>::reverse_bytewise(b);
      
      // Copy block to buffer:
      for (unsigned long a=0; a<BLOCK_WIDTH; ++a) {
	unsigned long a_rev = BitReversal<LOG_BLOCK_WIDTH>::reverse_bytewise(a);

	for (unsigned long c=0; c<BLOCK_WIDTH; ++c)
	  buffer[ (a_rev << LOG_BLOCK_WIDTH) | c ] = v[ (a << NUM_B_BITS << LOG_BLOCK_WIDTH) | (b << LOG_BLOCK_WIDTH) | c ];
      }

      // Swap from buffer:
      for (unsigned long c=0; c<BLOCK_WIDTH; ++c) {
      	// Note: Typo in original pseudocode by Carter and Gatlin at
      	// the following line:
      	unsigned long c_rev = BitReversal<LOG_BLOCK_WIDTH>::reverse_bytewise(c);
	
      	for (unsigned long a_rev=0; a_rev<BLOCK_WIDTH; ++a_rev)
      	  result[(c_rev << NUM_B_BITS << LOG_BLOCK_WIDTH) | (b_rev<<LOG_BLOCK_WIDTH) | a_rev] = buffer[ (a_rev<<LOG_BLOCK_WIDTH) | c ];
      }
    }
    free(buffer);

    memcpy(v, result, (1ul<<(LOG_N))*sizeof(T));
    free(result);
  }
};

#endif
