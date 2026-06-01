#ifndef _CONSTANTMULTIPLIERDEPENDENCY_HPP
#define _CONSTANTMULTIPLIERDEPENDENCY_HPP

#include "../Engine/ConstantMultiplierMessagePasser.hpp"
#include "Dependency.hpp"

// For building dependencies of the form {Y0,Y1,...} = {X0,X1,...} * {s0,s1,...}.

// dithering_sigma is used when the outcomes map to floating point
// values. Then, that mass is distributed between the two neighboring
// integer bins using a Gaussian with std. deviation dithering_sigma.

// Interpolation is used when multiplying a factor >1. For example a
// distribution with integer support [0,10] scaled by 3 should give a
// distribution with support {0, 3, 6, ..., 30}. However, when a
// discrete distribution is being used as a proxy for a continuous
// distribution, then each bin is not simply a point value, but a
// distribution in [0,1), [1,2), ..., which scaled by 3 will yield
// [0,3), [3,6), ... For this reason, ConstantMultiplierDependency
// allows the option to interpolate when scaling (whether forwards or
// backwards). If the input distribution is truly discrete (not a
// proxy for a continuous distribution), then setting
// interpolate_scaled=false is appropriate. If the output distribution
// is truly discrete, then setting interpolate_unscaled=false is
// appropriate. Note that in the future, these could be handled by the
// distributions themselves (by specifying continuous vs. discrete
// distribution types).

template <typename VARIABLE_KEY>
class ConstantMultiplierDependency : public Dependency<VARIABLE_KEY> {
protected:
  std::vector<VARIABLE_KEY> _input;
  std::vector<VARIABLE_KEY> _output;
  Vector<double> _scale;
  bool _interpolate_scaled, _interpolate_unscaled;
  double _dithering_sigma;
  
public:
  ConstantMultiplierDependency(const std::vector<VARIABLE_KEY> & input, const std::vector<VARIABLE_KEY> & output, Vector<double> scale, bool interpolate_scaled, bool interpolate_unscaled, double dithering_sigma):
    _input(input),
    _output(output),
    _scale(scale),
    _interpolate_scaled(interpolate_scaled),
    _interpolate_unscaled(interpolate_unscaled),
    _dithering_sigma(dithering_sigma)
  {
    #ifdef ENGINE_CHECK
    assert(input.size() == output.size() && input.size() == scale.size() && "Dimension of input, output, and scale in constant multiplier dependency must match");
    #endif
  }

  virtual ConstantMultiplierMessagePasser<VARIABLE_KEY>* create_message_passer(InferenceGraphBuilder<VARIABLE_KEY> & igb) const {
    HUGINMessagePasser<VARIABLE_KEY>* hyperedge_in = igb.create_hyperedge();
    HUGINMessagePasser<VARIABLE_KEY>* hyperedge_out = igb.create_hyperedge();

    std::vector<VARIABLE_KEY>* edge_label_in = new std::vector<VARIABLE_KEY>(_input);
    std::vector<VARIABLE_KEY>* edge_label_out = new std::vector<VARIABLE_KEY>(_output);
    
    return new ConstantMultiplierMessagePasser<VARIABLE_KEY>(hyperedge_in, edge_label_in, hyperedge_out, edge_label_out, _scale, _interpolate_scaled, _interpolate_unscaled, _dithering_sigma);
  }

  virtual std::vector<VARIABLE_KEY> get_all_variables_used() const {
    std::vector<VARIABLE_KEY> result = _input;
    for (const VARIABLE_KEY & var : _output)
      result.push_back(var);

    return result;
  }
};
#endif
