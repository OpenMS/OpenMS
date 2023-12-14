#ifndef _HUGINMESSAGEPASSER_HPP
#define _HUGINMESSAGEPASSER_HPP

#include "PNormMixin.hpp"
#include "ContextFreeMessagePasser.hpp"

// Note: There is an unexploited speedup when there are only two
// edges. In that case, it is better to solve it as a Shafer-Shenoy,
// where nothing is cached, since the outgoing messages will simply be
// raw inputs (with no multiplication). This could be done in
// HUGINMessagePasser or could alternatively be the responsibility of
// the graph construction (which would choose an alternate
// Shafer-Shenoy message passer in such situations).

template <typename VARIABLE_KEY>
class HUGINMessagePasser : public ContextFreeMessagePasser<VARIABLE_KEY>, public PNormMixin {
protected:
  LabeledPMF<VARIABLE_KEY> _prior;
  LabeledPMF<VARIABLE_KEY> _product;
  std::vector<LabeledPMF<VARIABLE_KEY> > _last_messages_received;
  std::vector<bool> _ready_to_send_ab_initio;

  void add_input_and_output_edges(Edge<VARIABLE_KEY>*edge_in, Edge<VARIABLE_KEY>*edge_out) {
    MessagePasser<VARIABLE_KEY>::add_input_and_output_edges(edge_in, edge_out);

    // Push an empty last message received:
    _last_messages_received.push_back(LabeledPMF<VARIABLE_KEY>());

    // When edge labels are subset of variables in _product,
    // start as ready to send ab initio:
    bool can_send_on_construction = true;
    for (const VARIABLE_KEY & var: *edge_in->variables_ptr)
      can_send_on_construction &= _product.contains_variable(var);

    _ready_to_send_ab_initio.push_back(can_send_on_construction);
  }

  bool ready_to_send_message_ab_initio(unsigned long edge_index) const {
    return _ready_to_send_ab_initio[edge_index];
  }

  void receive_message_in(unsigned long edge_index) {
    #ifndef SHAFERSHENOY
    Edge<VARIABLE_KEY>*incoming_edge = MessagePasser<VARIABLE_KEY>::_edges_in[edge_index];

    if (_product.dimension() > 0) {
      if (_last_messages_received[edge_index].dimension()>0)
	// A message had previously been received along edge; divide
	// out the old one and multiply in the new one:
	_product = incoming_edge->get_message() * _product / _last_messages_received[edge_index];
      else
	// No message had previously been received along edge;
	// multiply in the new one:
	_product = _product * incoming_edge->get_message();
    }
    else
      // No previous messages were received until now, and no prior
      // was provided; initialize the product distribution.
      _product = incoming_edge->get_message();

    _last_messages_received[edge_index] = incoming_edge->get_message();
    #endif
  }

  LabeledPMF<VARIABLE_KEY> get_message_out(unsigned long edge_index) {
    #ifdef SHAFERSHENOY
    LabeledPMF<VARIABLE_KEY> result = _product;
    for (unsigned int i=0; i<this->number_edges(); ++i) {
      if (i != edge_index) {
	Edge<VARIABLE_KEY>*e = this->_edges_in[i];
	if (e->has_message()) {
	  auto msg = e->get_possibly_outdated_message();
	  result = result * msg;
	}
      }
    }
    Edge<VARIABLE_KEY>*outward_edge = MessagePasser<VARIABLE_KEY>::_edges_out[edge_index];
    return result.marginal(*outward_edge->variables_ptr, this->p);
    #else
    Edge<VARIABLE_KEY>*outward_edge = MessagePasser<VARIABLE_KEY>::_edges_out[edge_index];

    LabeledPMF<VARIABLE_KEY> message_out = _product.marginal(*outward_edge->variables_ptr, this->p);

    // Divide out the message along the inward edge:
    if (this->_edge_received[edge_index]) {
      const LabeledPMF<VARIABLE_KEY> & message_in = _last_messages_received[edge_index];
      message_out = message_out / message_in;
    }

    return message_out;
    #endif
  }

public:
  HUGINMessagePasser(const LabeledPMF<VARIABLE_KEY> & prior, const double p):
    PNormMixin(p),
    _prior(prior),
    _product(prior)
  { }

  HUGINMessagePasser(LabeledPMF<VARIABLE_KEY> && prior, const double p):
    PNormMixin(p),
    _prior(std::move(prior)),
    _product(_prior)
  { }

  HUGINMessagePasser(const double p):
    PNormMixin(p)
  { }

  const LabeledPMF<VARIABLE_KEY> & joint_posterior() const {
    #ifdef SHAFERSHENOY
    LabeledPMF<VARIABLE_KEY> result;
    for (unsigned int i=0; i<this->number_edges(); ++i) {
      Edge<VARIABLE_KEY>*e = this->_edges_in[i];
      if (e->has_message()) {
	auto msg = e->get_possibly_outdated_message();
	result = result * msg;
      }
    }
    LabeledPMF<VARIABLE_KEY> * res = new LabeledPMF<VARIABLE_KEY>(result*_product);
    for (unsigned int i=0; i<this->number_edges(); ++i) {
      Edge<VARIABLE_KEY>*e = this->_edges_in[i];
      if (e->has_message()) {
	auto msg = e->get_possibly_outdated_message();
      }
    }
    return *res;
    #else
    return _product;
    #endif
  }

  const LabeledPMF<VARIABLE_KEY> & prior() const {
    return _prior;
  }

  void print(std::ostream & os) const {
    os << "HUGINMessagePasser prior:" << _prior << " _joint:" << _product;
  }
};

#endif
