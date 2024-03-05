#ifndef _HYPEREDGE_HPP
#define _HYPEREDGE_HPP

#include "HUGINMessagePasser.hpp"

// Like a HUGINMessagePasser, but it is elligible to pass once every
// variable along an edge has been received.
template <typename VARIABLE_KEY>
class Hyperedge : public HUGINMessagePasser<VARIABLE_KEY> {
private:
  std::unordered_set<VARIABLE_KEY> _vars_received;
  std::vector<bool> _ready_to_send;

  bool _all_ready_to_send;

protected:
  void add_input_and_output_edges(Edge<VARIABLE_KEY>*edge_in, Edge<VARIABLE_KEY>*edge_out) {
    HUGINMessagePasser<VARIABLE_KEY>::add_input_and_output_edges(edge_in, edge_out);

    _ready_to_send.push_back(false);
  }

  // Note that this relies on the fact that MessagePasser updates by
  // performing _ready_to_send[i] = _ready_to_send[i] | other_edges_received;
  // therefore, setting _ready_to_send[edge_index] here will ensure
  // that when _vars_received is a superset of the variables along
  // the edge, then the edge will be marked as ready.
  void receive_message_in(unsigned long edge_index) {
    HUGINMessagePasser<VARIABLE_KEY>::receive_message_in(edge_index);

    // Set edges out as ready to send where appropriate.
    if (! _all_ready_to_send) {
      // For greater performance, don't bother updating if this edge has
      // already been received. 
      if (! this->_edge_received[edge_index]) {
	// Add the variables to the set _vars_received.
	Edge<VARIABLE_KEY>*incoming_edge = this->_edges_in[edge_index];
      
	for (const VARIABLE_KEY & var : *incoming_edge->variables_ptr)
	  _vars_received.insert(var);
      
	for (unsigned long i=0; i<this->number_edges(); ++i) {
	  // Don't bother waking edge opposite to the message received (it
	  // will by definition be elligible to send, but nothing will
	  // have changed since the message received will not be used to
	  // send back).
	  if (i != edge_index) {
	    bool vars_received_are_superset = true;
	  
	    Edge<VARIABLE_KEY>*other_edge = this->_edges_in[i];
	    for (const VARIABLE_KEY & var : *other_edge->variables_ptr)
	      vars_received_are_superset = vars_received_are_superset && _vars_received.find(var) != _vars_received.end();

	    _ready_to_send[i] = vars_received_are_superset;
	  }
	}

	_all_ready_to_send = true;
	for (unsigned long i=0; i<this->number_edges(); ++i)
	  _all_ready_to_send = _all_ready_to_send && _ready_to_send[i];
      }
    }
  }

  virtual bool ready_to_send_message(unsigned long edge_index) const {
    return _ready_to_send[edge_index];
  }

  bool can_potentially_pass_any_messages() const {
    return true;
  }
  
public:
  Hyperedge():
    // Hyperedges use p=1.0; they exist solely to cache products via
    // the HUGIN algorithm.
    HUGINMessagePasser<VARIABLE_KEY>(1.0),
    _all_ready_to_send(false)
  { }

  void absorb_hyperedge(Hyperedge<VARIABLE_KEY>* he_to_absorb) {
    // Add edges from he_to_absorb into this:
    for (unsigned long i=0; i<he_to_absorb->number_edges(); ++i) {
      Edge<VARIABLE_KEY>*edge = he_to_absorb->get_edge_out(i);

      MessagePasser<VARIABLE_KEY>*dest_mp = edge->dest;

      if (dest_mp != this) {
	unsigned long source_edge_index = this->number_edges();
	unsigned long dest_edge_index = edge->dest_edge_index;
	
	Edge<VARIABLE_KEY>*edge_in = new Edge<VARIABLE_KEY>(dest_mp, this, edge->variables_ptr, dest_edge_index, source_edge_index);
	Edge<VARIABLE_KEY>*edge_out = new Edge<VARIABLE_KEY>(this, dest_mp, edge->variables_ptr, source_edge_index, dest_edge_index);
	this->add_input_and_output_edges(edge_in, edge_out);
	
	// Edge into this becomes edge out of dest_mp and vice versa:
	dest_mp->rewire_edge(edge->dest_edge_index, edge_out, edge_in);
      }
    }

    delete he_to_absorb;
  }

  void print(std::ostream & os) const {
    os << "Hyperedge " << this->_product;
  }
};

#endif
