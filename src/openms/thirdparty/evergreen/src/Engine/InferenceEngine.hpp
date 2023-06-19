#ifndef _INFERENCEENGINE_HPP
#define _INFERENCEENGINE_HPP

#include "../PMF/LabeledPMF.hpp"

template <typename VARIABLE_KEY>
class InferenceEngine {
public:
  virtual std::vector<LabeledPMF<VARIABLE_KEY> > estimate_posteriors(const std::vector<std::vector<VARIABLE_KEY> > & joint_distributions_to_retrieve) = 0;

  virtual double log_normalization_constant() = 0;

  virtual ~InferenceEngine() {}
};

#endif
