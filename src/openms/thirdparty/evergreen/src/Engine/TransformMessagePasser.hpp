#ifndef _TRANSFORMMESSAGEPASSER_HPP
#define _TRANSFORMMESSAGEPASSER_HPP

#include "ContextFreeMessagePasser.hpp"

template <typename VARIABLE_KEY>
class TransformMessagePasser : public MessagePasser<VARIABLE_KEY> {
protected:
  void receive_message_in(unsigned long index) {
    // Reorder to match _ordered_variables and pass into _ct:
    _ct.receive_message_in(index, this->_edges_in[index]->get_message().get_pmf());
  }
  LabeledPMF<VARIABLE_KEY> get_message_out(unsigned long index) {
    // Get message from _ct and then label appropriately:
    PMF pmf = _ct.get_message_out(index);

    return LabeledPMF<VARIABLE_KEY>(*this->_edges_out[index]->variables_ptr, std::move(pmf));
  }
  virtual void map_input_outcome_to_output_outcome(const Vector<unsigned long> & input_outcome, Vector<unsigned long> & output_outcome) = 0;
  virtual void map_output_outcome_to_input_outcome(Vector<unsigned long> & input_outcome, const Vector<unsigned long> & output_outcome) = 0;

public:
  TransformMessagePasser(const ContextFreeMessagePasser<VARIABLE_KEY>* & input, const std::vector<std::vector<VARIABLE_KEY>*> & input_edge_label, ContextFreeMessagePasser<VARIABLE_KEY>*output, std::vector<VARIABLE_KEY>* output_edge_label, const unsigned char dimension_param, const double p):
    MessagePasser<VARIABLE_KEY>(p),
    _ct(inputs.size(), dimension_param, p)
  {
    #ifdef ENGINE_CHECK
    assert(inputs.size() == input_edge.size());
    assert(output.size() == output_edge.size());
    assert(input.size() == output.size());
    #endif

    // Bind input first and output last:
    this->bind_to(input, input_edge_label);
    this->bind_to(output, output_edge_label);
  }
};

#endif
