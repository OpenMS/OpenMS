#ifndef _EVERGREEN_HPP
#define _EVERGREEN_HPP


// Create a macro for ALWAYS_INLINE (different modifiers on different compilers)
#if defined(_MSC_VER)
#define EVERGREEN_ALWAYS_INLINE __forceinline
#else
#define EVERGREEN_ALWAYS_INLINE inline __attribute__((always_inline))
#endif

// added by jpfeuffer to throw exceptions
#include <stdexcept>
#include <sstream>

// added by jpfeuffer to enable header guards before the namespace evergreen
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <array>
#include <set>
#include <cmath>
#include <list>
#include <map>
#include <utility>
#include <algorithm>
#include <functional>
#include <assert.h>
#include <string.h>

// Can be commented out for greater performance; these simply verify
// shapes of things and will not be very expensive (they are not like
// bounds checking). They should be left as-is by default, because
// they will warn of problems that could be tricky for novice users.
#ifdef NDEBUG
// nondebug
#else
// debug code
	#define SHAPE_CHECK
	#define NUMERIC_CHECK
	#define ENGINE_CHECK
#endif

// For debugging only (substantially decreases performance):
//#define BOUNDS_CHECK

namespace evergreen
{
  // for convenience:
  #include "../Utility/vector_ostream.hpp"

  // Inference engines:
  #include "../Engine/BeliefPropagationInferenceEngine.hpp"
  #include "BruteForceInferenceEngine.hpp"

  // Standard schedulers:
  #include "../Engine/PriorityScheduler.hpp"
  #include "../Engine/FIFOScheduler.hpp"
  #include "../Engine/RandomSubtreeScheduler.hpp"

  // Standard dependencies:
  #include "AdditiveDependency.hpp"
  #include "PseudoAdditiveDependency.hpp"
  #include "ConstantMultiplierDependency.hpp"
  // TableDependency will already be included.

  #include "BetheInferenceGraphBuilder.hpp"
}


#endif
