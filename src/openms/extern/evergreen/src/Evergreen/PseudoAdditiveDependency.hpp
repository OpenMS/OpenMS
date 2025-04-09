#ifndef _PSEUDOADDITIVEDEPENDENCY_HPP
#define _PSEUDOADDITIVEDEPENDENCY_HPP

#include <map>
#include <limits>
#include "../Engine/ConvolutionTreeMessagePasser.hpp"
#include "TableDependency.hpp"

// PseudoAdditiveDependency behaves similar to AdditiveDependency, but
// in contrast, it works with brute force, creating a joint
// PMF. Because it creates this PMF, it is constructed with tables,
// which it uses to uncover the support of each variable can
// obtain. These variable supports (which are the intersection of
// supports of the tables provided) determine the supports of each
// axis of the table. PseudoAdditiveDependency is SLOW. It is intended
// for testing and debugging.

// Note: This implementation is not lighting fast, and regardless,
// AdditiveDependency will be much more efficient (this creates a
// table for the additive dependency). It is to be used for testing.
template <typename VARIABLE_KEY>
class PseudoAdditiveDependency : public Dependency<VARIABLE_KEY>, public PNormMixin {
protected:
  std::vector<std::vector<VARIABLE_KEY> > _inputs;
  std::vector<VARIABLE_KEY> _output;

  std::map<VARIABLE_KEY, std::array<long, 2> > _var_to_min_and_max;

  // Check if outputs are sums of inputs:
  bool is_additive(std::map<VARIABLE_KEY, long> var_to_outcome) const {
    for (unsigned int i=0; i<_output.size(); ++i) {
      const VARIABLE_KEY & var = _output[i];

      long res = var_to_outcome[var];
      
      long tot = 0;
      for (unsigned long k=0; k<_inputs.size(); ++k)
	tot += var_to_outcome[ _inputs[k][i] ];

      if ( res != tot )
	return false;
    }
    return true;
  }

public:
  PseudoAdditiveDependency(const std::vector<std::vector<VARIABLE_KEY> > & inputs, const std::vector<VARIABLE_KEY> & output, const std::vector<TableDependency<VARIABLE_KEY> > & existing_tables, const double p_param):
    PNormMixin(p_param),
    _inputs(inputs),
    _output(output)
  {
    #ifdef ENGINE_CHECK
    for (const std::vector<VARIABLE_KEY> & in : _inputs)
      assert(in.size() == output.size() && "Dimension of each tuple in additive dependency must equal dimension of output tuple");
    #endif

    // For all variables in the additive dependency, initialize
    // minimum and maximum values:
    for (const std::vector<VARIABLE_KEY> & in : inputs)
      for (const VARIABLE_KEY & var : in) {
	auto iter = _var_to_min_and_max.find(var);
	if (iter == _var_to_min_and_max.end())
	  _var_to_min_and_max[var] = {{std::numeric_limits<long>::min(), std::numeric_limits<long>::max()}};
      }
    for (const VARIABLE_KEY & var : output) {
      auto iter = _var_to_min_and_max.find(var);
      if (iter == _var_to_min_and_max.end())
	_var_to_min_and_max[var] = {{std::numeric_limits<long>::min(), std::numeric_limits<long>::max()}};
    }
    
    // Find the bounding box for all variables in the additive
    // dependency:
    for (const TableDependency<VARIABLE_KEY> & tab_dep : existing_tables) {
      auto lpmf = tab_dep.labeled_pmf();
      for (unsigned int i=0; i<lpmf.dimension(); ++i) {
	const VARIABLE_KEY & var = lpmf.ordered_variables()[i];
	long min_val = lpmf.pmf().first_support()[i];
	long max_val = min_val + lpmf.pmf().table().view_shape()[i];
	_var_to_min_and_max[var][0] = std::max(_var_to_min_and_max[var][0], min_val);
	_var_to_min_and_max[var][1] = std::min(_var_to_min_and_max[var][1], max_val);
      }
    }
  }

  LabeledPMF<VARIABLE_KEY> to_labeled_pmf() const {
    std::vector<VARIABLE_KEY> ordered_variables(_var_to_min_and_max.size());
    Vector<long> first_support(_var_to_min_and_max.size());
    Vector<unsigned long> shape(_var_to_min_and_max.size());
    unsigned long i=0;
    for (const auto & var_min_max : _var_to_min_and_max) {
      ordered_variables[i] = var_min_max.first;
      first_support[i] = var_min_max.second[0];
      shape[i] = var_min_max.second[1] - var_min_max.second[0] + 1;
      ++i;
    }

    Tensor<double> ten(shape);
    std::map<VARIABLE_KEY, long> var_to_outcome;
    enumerate_apply_tensors([this, &ordered_variables, &var_to_outcome, &first_support](const_tup_t counter, const unsigned char dim, double & val){
	for (unsigned int i=0; i<dim; ++i)
	  var_to_outcome[ordered_variables[i]] = first_support[i] + counter[i];

	if (is_additive(var_to_outcome))
	  val = 1.0;
      },
      ten.data_shape(),
      ten);
    
    PMF pmf(first_support, ten);
    return LabeledPMF<VARIABLE_KEY>(ordered_variables, pmf);
  }

  TableDependency<VARIABLE_KEY> to_table_dependency() const {
    return TableDependency<VARIABLE_KEY>(to_labeled_pmf(), this->p);
  }

  virtual HUGINMessagePasser<VARIABLE_KEY>* create_message_passer(InferenceGraphBuilder<VARIABLE_KEY> & igb) const {
    return to_table_dependency().create_message_passer(igb);
  }
};
#endif
