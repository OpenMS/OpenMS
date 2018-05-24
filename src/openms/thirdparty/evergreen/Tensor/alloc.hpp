#ifndef _ALLOC_HPP
#define _ALLOC_HPP

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <alloca.h>

// Note: may benefit from being tuned for specific architecture:
const unsigned long int ALLOCATION_ALIGNMENT = 128;

template <typename T>
T* aligned_malloc(unsigned long num_elements) {
  // TODO: the aligned allocation code below was meant to optimize
  // cache alignment, but it can be quite slow. malloc is used as an
  // alternative; however, the aligned_calloc function should be used
  // as the interface so that all invocations can be updated if a
  // superior solution is found.
  T*result = (T*)malloc(num_elements*sizeof(T));
  assert(result != NULL);
  return result;

  /*
  void* result = NULL;
  if (num_elements > 0) {
    bool success = !posix_memalign(&result, ALLOCATION_ALIGNMENT, num_elements*sizeof(T));
    assert( success );
  }
  return (T*)result;
  */
}

template <typename T>
T* aligned_calloc(unsigned long num_elements) {
  T*result = (T*)calloc(num_elements, sizeof(T));
  assert(result != NULL);
  return result;

  // Note: see above.
  /*
  void* result = NULL;
  if (num_elements > 0) {
    bool success = !posix_memalign(&result, ALLOCATION_ALIGNMENT, num_elements*sizeof(T));
    assert( success );
    memset(result, 0, num_elements*sizeof(T));
  }
  return (T*)result;
  */
}

// Variable length array stack allocation:
template <typename T>
T* vla_alloc(unsigned long num_elements) {
  return (T*)alloca(num_elements*sizeof(T));
}

#endif
