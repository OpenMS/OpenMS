#ifndef _ALLOC_HPP
#define _ALLOC_HPP

//#include <stdlib.h>
//#include <stdio.h>
#include <string.h>
#include <assert.h>
#ifdef _MSC_VER
	#include <malloc.h>
	#define alloca _alloca
#else
	#include <alloca.h>
#endif

// Note: benefits from being tuned for specific architecture. Could do this with #ifdef...s for AVX512, etc.
const unsigned long int ALLOCATION_ALIGNMENT = 512;

template <typename T>
T* aligned_malloc(unsigned long num_elements) {
  /* Aligned malloc makes operations safe for AVX, etc. but results in a slowdown of roughly 3x
  #ifdef _WIN32
    T*result = (T*)_aligned_malloc(ALLOCATION_ALIGNMENT, num_elements*sizeof(T));
    assert(result != NULL);
    return result;
  #else
    void* result = NULL;
    posix_memalign(&result, ALLOCATION_ALIGNMENT, num_elements*sizeof(T));
    assert(result != NULL);
    return (T*)result;
  #endif
  */
  T*result = (T*)malloc(num_elements*sizeof(T));
  assert(result != NULL);
  return result;
}

template <typename T>
T* aligned_calloc(unsigned long num_elements) {
  T*result = aligned_malloc<T>(num_elements);
  assert(result != NULL);
  memset(result, 0, sizeof(T)*num_elements);
  return result;
}

// Variable length array stack allocation:
template <typename T>
T* vla_alloc(unsigned long num_elements) {
  return (T*)alloca(num_elements*sizeof(T));
}

#endif
