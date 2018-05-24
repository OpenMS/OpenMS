#ifndef _BRUTEFORCEINFERENCEENGINE_HPP
#define _BRUTEFORCEINFERENCEENGINE_HPP

#include "../Engine/InferenceEngine.hpp"
#include "TableDependency.hpp"

template <typename VARIABLE_KEY>
class BruteForceInferenceEngine : public InferenceEngine<VARIABLE_KEY> {
protected:
  LabeledPMF<VARIABLE_KEY> _joint;
  const double _p;

public:
  BruteForceInferenceEngine(const std::vector<TableDependency<VARIABLE_KEY> > & all_tables, double p):
    _p(p)
  {
    // Note: This could be performed more efficiently by first
    // allocating the result table (using the union of all variables
    // and the intersection of their supports) and then performing a
    // single product over all distributions simultaneously.
    for (const TableDependency<VARIABLE_KEY> & table : all_tables)
      _joint = _joint * table.labeled_pmf();
  }

  std::vector<LabeledPMF<VARIABLE_KEY> > estimate_posteriors(const std::vector<std::vector<VARIABLE_KEY> > & joint_distributions_to_retrieve) {
    std::vector<LabeledPMF<VARIABLE_KEY> > results;
    for (const std::vector<VARIABLE_KEY> & ordered_vars : joint_distributions_to_retrieve)
      results.push_back(_joint.marginal(ordered_vars, _p));
    return results;
  }
};

#endif
