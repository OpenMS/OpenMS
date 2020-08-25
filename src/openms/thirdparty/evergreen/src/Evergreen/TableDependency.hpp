#ifndef _TABLEDEPENDENCY_HPP
#define _TABLEDEPENDENCY_HPP

#include "Dependency.hpp"
#include "../PMF/LabeledPMF.hpp"

template <typename VARIABLE_KEY>
class TableDependency : public Dependency<VARIABLE_KEY>, public PNormMixin {
protected:
  LabeledPMF<VARIABLE_KEY> _lpmf;
  
public:
  TableDependency(const LabeledPMF<VARIABLE_KEY> & lpmf, const double p_param):
    PNormMixin(p_param),
    _lpmf(lpmf)
  { }

  virtual HUGINMessagePasser<VARIABLE_KEY>* create_message_passer(InferenceGraphBuilder<VARIABLE_KEY> & igb) const {
    // Note: Does not create hyperdges or bind to hyperedges. That is
    // the responsibility of InferenceGraphBuilder.
    return new HUGINMessagePasser<VARIABLE_KEY>(_lpmf, this->p);
  }

  virtual std::vector<VARIABLE_KEY> get_all_variables_used() const {
    return _lpmf.ordered_variables();
  }

  const LabeledPMF<VARIABLE_KEY> & labeled_pmf() const {
    return _lpmf;
  }
};

#endif
