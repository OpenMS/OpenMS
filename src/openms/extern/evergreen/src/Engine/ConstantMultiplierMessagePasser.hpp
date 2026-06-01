#ifndef _CONSTANTMULTIPLIERMESSAGEPASSER_HPP
#define _CONSTANTMULTIPLIERMESSAGEPASSER_HPP

#include "ContextFreeMessagePasser.hpp"

template <typename VARIABLE_KEY>
class ConstantMultiplierMessagePasser: public MessagePasser<VARIABLE_KEY> {
protected:
  Vector<double> _scale;
  Vector<double> _one_over_scale;
  PMF _message_in_received;
  PMF _message_out_received;
  bool _interpolate_scaled, _interpolate_unscaled;
  double _dithering_sigma_squared;

  void receive_message_in(unsigned long index) {
    Edge<VARIABLE_KEY>*incoming_edge = MessagePasser<VARIABLE_KEY>::_edges_in[index];
    if (index == 0) {
      // Message received by input:
      _message_in_received = incoming_edge->get_message().pmf();
    }
    else {
      // Message received by output:
      _message_out_received = incoming_edge->get_message().pmf();
    }
  }

  LabeledPMF<VARIABLE_KEY> get_message_out(unsigned long index) {
    Edge<VARIABLE_KEY>*incoming_edge = MessagePasser<VARIABLE_KEY>::_edges_in[index];
    if (index == 0) {
      // Get message out through input (multiply _message_out_received
      // by _one_over_scale):
      PMF scaled;
      if (_interpolate_unscaled)
	scaled = scaled_pmf_dither_interpolate(_message_out_received, _one_over_scale, _dithering_sigma_squared);
      else
	scaled = scaled_pmf_dither(_message_out_received, _one_over_scale, _dithering_sigma_squared);
      
      if (_message_in_received.dimension() > 0)
	scaled.narrow_support(_message_in_received.first_support(), _message_in_received.last_support());

      return LabeledPMF<VARIABLE_KEY>(*incoming_edge->variables_ptr, scaled);
    }
    else {
      // Get message out through output (multiply _message_in_received
      // by _scale):
      PMF scaled;
      if (_interpolate_scaled)
	scaled = scaled_pmf_dither_interpolate(_message_in_received, _scale, _dithering_sigma_squared);
      else
	scaled = scaled_pmf_dither(_message_in_received, _scale, _dithering_sigma_squared);

      if (_message_out_received.dimension() > 0)
	scaled.narrow_support(_message_out_received.first_support(), _message_out_received.last_support());

      return LabeledPMF<VARIABLE_KEY>(*incoming_edge->variables_ptr, scaled);
    }
  }

public:
  ConstantMultiplierMessagePasser(ContextFreeMessagePasser<VARIABLE_KEY>* input, const std::vector<VARIABLE_KEY>* input_edge_label, ContextFreeMessagePasser<VARIABLE_KEY>*output, std::vector<VARIABLE_KEY>* output_edge_label, const Vector<double> & scale_param, bool interpolate_scaled, bool interpolate_unscaled, double dithering_sigma):
    _scale(scale_param),
    _one_over_scale(1.0 / scale_param),
    _interpolate_scaled(interpolate_scaled),
    _interpolate_unscaled(interpolate_unscaled),
    _dithering_sigma_squared(dithering_sigma*dithering_sigma)
  {
    #ifdef ENGINE_CHECK
    assert(input_edge_label->size() == output_edge_label->size());
    assert(input_edge_label->size() == scale_param.size());
    #endif

    // Bind input first and output last:
    this->bind_to(input, input_edge_label);
    this->bind_to(output, output_edge_label);
  }

  void print(std::ostream & os) const {
    Edge<VARIABLE_KEY>*input_edge = MessagePasser<VARIABLE_KEY>::_edges_in[0];
    Edge<VARIABLE_KEY>*output_edge = MessagePasser<VARIABLE_KEY>::_edges_in[1];

    os << "ConstantMultiplierMessagePasser ";
    for (unsigned int i=0; i<output_edge->variables_ptr->size(); ++i)
      os << (*output_edge->variables_ptr)[i] << " ";
    os << "= " << _scale << " * ";
    for (unsigned int i=0; i<input_edge->variables_ptr->size(); ++i)
      os << (*input_edge->variables_ptr)[i] << " ";
  }

  const Vector<double> & scale() const {
    return _scale;
  }
};

#endif
