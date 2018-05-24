#ifndef _TEMPLATESEARCH_HPP
#define _TEMPLATESEARCH_HPP

#include <assert.h>
#include <utility>

typedef unsigned char TEMPLATE_SEARCH_INT_TYPE;

// For dynamically looking up a class and calling the static
// function apply(...). This can be preferred to LogSearch when you
// are amortizing the cost of lookoup. For instance, if the log size
// of an FFT is being looked up, then proceeding in ascending order
// from 0 will guarantee that the search cost is linear in the log
// size, meaning that it's in O(log(N)), which is dwarfed by the FFT
// cost O(N log(N)). This can therefore be more efficient when
// you're processing many short FFTs.
template <TEMPLATE_SEARCH_INT_TYPE MINIMUM, TEMPLATE_SEARCH_INT_TYPE MAXIMUM, template <TEMPLATE_SEARCH_INT_TYPE> class WORKER>
class LinearTemplateSearch {
public:
  template <typename...ARG_TYPES>
  inline static void apply(TEMPLATE_SEARCH_INT_TYPE v, ARG_TYPES && ... args) {
    if (v == MINIMUM)
      WORKER<MINIMUM>::apply(std::forward<ARG_TYPES>(args)...);
    else
      LinearTemplateSearch<MINIMUM+1, MAXIMUM, WORKER>::apply(v, std::forward<ARG_TYPES>(args)...);
  }
};

template <TEMPLATE_SEARCH_INT_TYPE MAXIMUM, template <TEMPLATE_SEARCH_INT_TYPE> class WORKER>
class LinearTemplateSearch<MAXIMUM, MAXIMUM, WORKER> {
public:
  template <typename...ARG_TYPES>
  inline static void apply(TEMPLATE_SEARCH_INT_TYPE v, ARG_TYPES && ... args) {
    assert(v == MAXIMUM);
    WORKER<MAXIMUM>::apply(std::forward<ARG_TYPES>(args)...);
  }
};

#endif
