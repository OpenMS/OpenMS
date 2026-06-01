#ifndef _CONVOLUTIONTREEMESSAGEPASSER_HPP
#define _CONVOLUTIONTREEMESSAGEPASSER_HPP

#include "PNormMixin.hpp"
#include "ContextFreeMessagePasser.hpp"
#include "ConvolutionTree.hpp"

template <typename VARIABLE_KEY>
class ConvolutionTreeMessagePasser : public MessagePasser<VARIABLE_KEY>, public PNormMixin {
protected:
  ConvolutionTree _ct;

  void receive_message_in(unsigned long index) {
    // Note that variable order within each PMF (which assigns axes)
    // will be determined by the edge label already (when the message
    // is passed). Therefore, the order of variables in each edge
    // label matters (not only the set of variables).
    _ct.receive_message_in(index, this->_edges_in[index]->get_message().pmf());
  }

  LabeledPMF<VARIABLE_KEY> get_message_out(unsigned long index) {
    // Get message from _ct and then label appropriately:
    PMF pmf = _ct.get_message_out(index);

    return LabeledPMF<VARIABLE_KEY>(*this->_edges_out[index]->variables_ptr, std::move(pmf));
  }

public:
  // Should only be constructed by binding to ContextFreeMessagePasser
  // types because it will involve adding edges to those MessagePasser
  // types, which could violate context.
  ConvolutionTreeMessagePasser(const std::vector<ContextFreeMessagePasser<VARIABLE_KEY>*> & inputs, const std::vector<std::vector<VARIABLE_KEY>*> & input_edge_labels, ContextFreeMessagePasser<VARIABLE_KEY>*output, std::vector<VARIABLE_KEY>* output_edge_label, const unsigned char dimension_param, const double p_param):
    PNormMixin(p_param),
    _ct(inputs.size(), dimension_param, p_param)
  {
    #ifdef ENGINE_CHECK
    assert(inputs.size() == input_edge_labels.size());

    // Note: dimensions are not verified, but will be when messages
    // are passed through edges.
    #endif

    // Bind inputs first and output last; this is what ConvolutionTree
    // is expecting:
    for (unsigned long k=0; k<inputs.size(); ++k)
      this->bind_to(inputs[k], input_edge_labels[k]);
    this->bind_to(output, output_edge_label);
  }

  virtual void print(std::ostream & os) const {
    // TODO: implement more informative version:
    os << "ConvolutionTreeMessagePasser " << int(_ct.dimension()) << " ";
    for (unsigned long i=0; i<this->_edges_in.size()-1; ++i) {
      os << "{ ";
      for (unsigned char j=0; j<_ct.dimension(); ++j) {
	os << (*this->_edges_in[i]->variables_ptr)[j];
	os << " ";
      }
      os << "}";

      if (i != this->_edges_in.size()-2)
	os << " + ";
    }

    os << " = { ";
    for (unsigned char j=0; j<_ct.dimension(); ++j) {
      os << (*this->_edges_in[this->_edges_in.size()-1]->variables_ptr)[j];
      os << " ";
    }
    os << "}";
  }
};

#endif
