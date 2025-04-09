#ifndef _L1REGULARIZATION_HPP
#define _L1REGULARIZATION_HPP

#include "to_string.hpp"

// Note: This is for regularizing 1D variables only thus far; it could
// be generalized (but regularizing with multidimensional indicator
// variables would only be interesting if the likelihood function on
// the sum of the indicators has a lot of covariation between axes).
template <typename T>
class L1Regularization {
public:
  static void apply(InferenceGraphBuilder<T> & igb, const std::vector<T> & vars_to_regularize, const std::vector<T> & indicator_vars, const LabeledPMF<T> & sum_of_indicators, double p, unsigned long prior_maximum_copies_of_element) {
    assert(vars_to_regularize.size() == indicator_vars.size() && "Variables and indicator variables should be paired and in order");

    // Insert the indicator variables:
    for (unsigned long i=0; i<vars_to_regularize.size(); ++i) {
      const T & var = vars_to_regularize[i];
      const T & indicator_var = indicator_vars[i];
      
      igb.insert_dependency( make_uniform_indicator_for_nonneg_var(var, indicator_var, prior_maximum_copies_of_element, p) );
    }

    // Insert the additive dependency:
    
    // Build vector of singletons for additive dependency:
    std::vector<std::vector<T> > indicator_singletons(indicator_vars.size());
    for (unsigned long i=0; i<indicator_vars.size(); ++i)
      indicator_singletons[i] = {indicator_vars[i]};

    igb.insert_dependency( AdditiveDependency<T>(indicator_singletons, sum_of_indicators.ordered_variables(), p) );

    // Insert the distribution on the sum of indicators:
    igb.insert_dependency(TableDependency<T>(sum_of_indicators, p));
  }

  // Make a 2D LabeledPMF that is 2 x max_val. Where indicator
  static TableDependency<T> make_uniform_indicator_for_nonneg_var(const T & var, const T & indicator_var, unsigned long max_val, double p) {
    Tensor<double> ten({max_val+1, 2});
    unsigned long ten_size = ten.flat_size();
    // if number of element != 0 then the indicator i_e = 1 is 100%
    for (unsigned long i=0; i<ten_size; ++i) {
      ten[i] = 0.0;
      ten[++i] = 1.0;
    }
    // if number of element = 0 then the indicator i_e = 0 is 100%
    ten.flat()[0] = 1.0;
    ten.flat()[1] = 0;
    
    // Create uniform pmf with indiciator.
    PMF pmf( {0L, 0L}, ten );
    return TableDependency<T>(LabeledPMF<T>({var, indicator_var}, pmf), p);
  }
};

#endif
