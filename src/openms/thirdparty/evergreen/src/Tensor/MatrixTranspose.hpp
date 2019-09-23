#ifndef _MATRIXTRANSPOSE_HPP
#define _MATRIXTRANSPOSE_HPP

#include <algorithm>

// Note: This code is a good candidate to perform on a GPU.

// Implements cache-oblivious strategy for transposition. In cases
// where recursion can simply be unrolled into a loop (e.g., tall,
// thin matrix or short, wide matrix) it is unrolled for greater
// performance.
template <typename T>
class MatrixTranspose {
private:
  // Base case for recursion; should fit in L1 cache. In practice, it
  // just needs to be large enough to amortize out the cost of the
  // recursions. 
  static constexpr unsigned int BLOCK_SIZE = 128 / sizeof(T);;

  // Square (in place):
  static void square_helper(T* __restrict const mat, const unsigned long N, const unsigned long r_start, const unsigned long r_end, const unsigned long c_start, const unsigned long c_end) {
    unsigned long r_span = r_end-r_start;
    unsigned long c_span = c_end-c_start;

    if ( c_span <= BLOCK_SIZE ) {
      // Tall, narrow block: proceed row-by-row:
      for (unsigned long r=r_start; r<r_end; ++r)
    	// Force c > r to not swap multiple times:
    	for (unsigned long c=std::max(c_start, r+1); c<c_end; ++c)
    	  std::swap(mat[c*N + r], mat[r*N + c]);
    }
    else if ( r_span <= BLOCK_SIZE ) {
      // Short, fat block: proceeding column-by-column will be optimal:
      for (unsigned long c=c_start; c<c_end; ++c)
    	// Force c > r to not swap multiple times:
    	for (unsigned long r=r_start; r<std::min(r_end, c); ++r)
    	  std::swap(mat[c*N + r], mat[r*N + c]);
    }
    else {
      if (r_span > c_span) {
      	// if there are any cells c>r for the first subproblem:
      	if (c_end > r_start)
      	  square_helper(mat, N, r_start,r_start+r_span/2, c_start,c_end);
	
      	// if there are any cells c>r for the second subproblem:
      	if (c_end > r_start+r_span/2)
      	  square_helper(mat, N, r_start+r_span/2,r_end, c_start,c_end);
      }
      else {
      	// if there are any cells c>r for the first subproblem:
      	if (c_start+c_span/2 > r_start)
      	  square_helper(mat, N, r_start,r_end, c_start,c_start+c_span/2);
	
      	// if there are any cells c>r for the second subproblem:
      	if (c_end > r_start)
      	  square_helper(mat, N, r_start,r_end, c_start+c_span/2,c_end);
      }
    }
  }

  // Buffered (out of place):
  static void buffered_helper(T* __restrict const dest, T* __restrict const source, const unsigned long R, const unsigned long C, const unsigned long r_start, const unsigned long r_end, const unsigned long c_start, const unsigned long c_end) {
    unsigned long r_span = r_end-r_start;
    unsigned long c_span = c_end-c_start;
    if ( c_span <= BLOCK_SIZE ) {
      // Tall, narrow block: proceed row-by-row:
      for (unsigned long r=r_start; r<r_end; ++r)
	for (unsigned long c=c_start; c<c_end; ++c)
	  // dest[c,r] = source[r,c];
	  dest[c*R + r] = source[r*C + c];
    }
    else if ( r_span <= BLOCK_SIZE ) {
      // Short, fat block: proceeding column-by-column will be optimal:
      for (unsigned long c=c_start; c<c_end; ++c)
	for (unsigned long r=r_start; r<r_end; ++r)
	  // dest[c,r] = source[r,c];
	  dest[c*R + r] = source[r*C + c];
    }
    else {
      if (r_span > c_span) {
	buffered_helper(dest, source, R, C, r_start,r_start+r_span/2, c_start,c_end);
	buffered_helper(dest, source, R, C, r_start+r_span/2,r_end, c_start,c_end);
      }
      else {
	buffered_helper(dest, source, R, C, r_start,r_end, c_start,c_start+c_span/2);
	buffered_helper(dest, source, R, C, r_start,r_end, c_start+c_span/2,c_end);
      }
    }
  }

public:

  inline static void apply_square(T* __restrict const mat, const unsigned long N) {
    square_helper(mat, N, 0, N, 0, N);
  }

  static void apply_square_naive(T* __restrict const mat, const unsigned long N) {
    for (unsigned long r=0; r<N; ++r)
      for (unsigned long c=r+1; c<N; ++c)
	std::swap(mat[r*N+c], mat[c*N+r]);
  }

  inline static void apply_buffered(T* __restrict const dest, T* __restrict const source, const unsigned long R, const unsigned long C) {
    buffered_helper(dest, source, R, C, 0, R, 0, C);
  }

  static void apply_buffered_naive(T* __restrict const dest, const T* __restrict const source, const unsigned long R, const unsigned long C) {
    for (unsigned long r=0; r<R; ++r)
      for (unsigned long c=0; c<C; ++c)
	dest[c*R+r] = source[r*C+c];
  }
};

#endif
