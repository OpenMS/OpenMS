#ifndef _EVERGREEN_HPP
#define _EVERGREEN_HPP

// Can be commented out for greater performance; these simply verify
// shapes of things and will not be very expensive (they are not like
// bounds checking). They should be left as-is by default, because
// they will warn of problems that could be tricky for novice users.
#define SHAPE_CHECK
#define NUMERIC_CHECK
#define ENGINE_CHECK
// For debugging only (substantially decreases performance):
//#define BOUNDS_CHECK

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

#endif
