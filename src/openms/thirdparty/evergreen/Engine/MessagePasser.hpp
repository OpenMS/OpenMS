#ifndef _MESSAGEPASSER_HPP
#define _MESSAGEPASSER_HPP

#include "Edge.hpp"
#include <unordered_set>

template <typename VARIABLE_KEY>
class ContextFreeMessagePasser;

// Interface for message passers in the engine:
template <typename VARIABLE_KEY>
class MessagePasser {
protected:
  // Note: _edges_in[i] must be the reverse edge of
  // _edges_out{i]. This is ensured by only modifying via
  // add_input_and_output_edges.

  // Note: The members _edges_in and _edges_out would ideally would be
  // private for the reason above.  However, the
  // add_input_and_output_edges function, which decides whether a
  // particular edge begins ready to pass or not, modifies them. For
  // this reason, it is simpler for them to be protected. However, at
  // a later date, it could be nice to refactor so that they're
  // private and there is a virtual bool ready_to_pass function that
  // is called when adding edges. This would also require that all
  // accesses of _edges_in and _edges_out were performed through
  // accessors.
  
  std::vector<Edge<VARIABLE_KEY>* > _edges_in;
  std::vector<Edge<VARIABLE_KEY>* > _edges_out;

  virtual void add_input_and_output_edges(Edge<VARIABLE_KEY>*edge_in, Edge<VARIABLE_KEY>*edge_out) {
    _edges_in.push_back(edge_in);
    _edges_out.push_back(edge_out);
    _edge_received.push_back(false);
  }

  std::vector<bool> _edge_received;

  // For determining the edges that need to be dirtied (possibly in
  // amortized O(1)):
  unsigned long _number_edges_with_messages_received;
  bool _all_edges_out_not_up_to_date;
  bool _all_edges_out_but_one_not_up_to_date;
  long _up_to_date_edge_if_one_exists;

  // Derived classes will override the following functions:
  virtual void receive_message_in(unsigned long edge_index) = 0;
  virtual LabeledPMF<VARIABLE_KEY> get_message_out(unsigned long edge_index) = 0;

  void update_after_receiving_message_in(unsigned long edge_index) {
    // Update which messages have been received and the count:
    if ( ! _edge_received[edge_index] ) {
      _edge_received[edge_index] = true;
      ++_number_edges_with_messages_received;
    }

    // Make local vars so that these can be modified below:
    bool all_not_up_to_date = _all_edges_out_not_up_to_date;
    bool all_but_this_one_not_up_to_date = _number_edges_with_messages_received > 0 && _all_edges_out_but_one_not_up_to_date && (_up_to_date_edge_if_one_exists == (long)edge_index);

    // after receiving a message, either all edges out or all edges
    // out but one are not up to date.
    if (_edges_out[edge_index]->up_to_date()) {
      _all_edges_out_not_up_to_date = false;

      _all_edges_out_but_one_not_up_to_date = true;
      _up_to_date_edge_if_one_exists = edge_index;
    }
    else {
      _all_edges_out_not_up_to_date = true;

      _all_edges_out_but_one_not_up_to_date = false;
      // value of _up_to_date_edge_if_one_exists does not matter if
      // _all_edges_out_not_up_to_date is true.
      _up_to_date_edge_if_one_exists = -1L;
    }

    // Don't bother dirtying edges if they were all already
    // dirty. Likewise, don't bother dirtying edges if all edges that
    // would be dirtied are already dirty.
    if ( ! all_not_up_to_date && ! all_but_this_one_not_up_to_date )
      for (unsigned long i=0; i<number_edges(); ++i)
	if (i != edge_index) {
	  // Messages out along every edge (except the opposite edge of
	  // the message being received) now become invalid:
	  
	  // Mark old messages along e as no longer valid:
	  _edges_out[i]->set_not_up_to_date();
	}

  }

  // Only allows rhs to be ContextFreeMessagePasser* so that calling
  // bind_to does not violate existing context.
  void bind_to(ContextFreeMessagePasser<VARIABLE_KEY>*rhs, const std::vector<VARIABLE_KEY>*const ordered_edge_vars) {
    unsigned long num_this_edges = number_edges();
    unsigned long num_rhs_edges = rhs->number_edges();

    // Edges should only ever be created as pairs (as below):
    Edge<VARIABLE_KEY>*edge = new Edge<VARIABLE_KEY>(this, rhs, ordered_edge_vars, num_this_edges, num_rhs_edges);
    Edge<VARIABLE_KEY>*opposite_edge = new Edge<VARIABLE_KEY>(rhs, this, ordered_edge_vars, num_rhs_edges, num_this_edges);

    add_input_and_output_edges(opposite_edge, edge);
    rhs->add_input_and_output_edges(edge, opposite_edge);
  }

  MessagePasser():
    _number_edges_with_messages_received(0),
    _all_edges_out_not_up_to_date(true),
    _all_edges_out_but_one_not_up_to_date(false),
    _up_to_date_edge_if_one_exists(-1L),
    color(0)
  { }

public:
  // To permit basic graph operations by marking in O(n):
  long color;

  // Note: costs Omega(n) each call, so shouldn't be called
  // frequently; however, it's a useful shorthand to prevent duplicate
  // code in other functions that repeatedly do the same thing.
  std::unordered_set<VARIABLE_KEY> variables_used_by_incident_edges() const {
    std::unordered_set<VARIABLE_KEY> result;
    for (const Edge<VARIABLE_KEY>*edge : _edges_in)
      for (const VARIABLE_KEY & var : *edge->variables_ptr)
	result.insert( var );
    return result;
  }

  virtual ~MessagePasser() {}

  void receive_message_in_and_update(unsigned long edge_index) {
    receive_message_in(edge_index);

    Edge<VARIABLE_KEY>*incoming_edge = _edges_in[edge_index];
    update_after_receiving_message_in(incoming_edge->dest_edge_index);
  }

  LabeledPMF<VARIABLE_KEY> update_and_get_message_out(unsigned long edge_index) {
    // Assume this will be used to set _edges_out[edge_index] is up-to-date.
    _all_edges_out_but_one_not_up_to_date = _all_edges_out_not_up_to_date;
    _up_to_date_edge_if_one_exists = edge_index;
    _all_edges_out_not_up_to_date = false;
    
    return get_message_out(edge_index);
  }

  // Note: excludes ab initio messages; to check ab initio, use
  // ready_to_send_message_ab_initio.
  virtual bool ready_to_send_message(unsigned long edge_index) const {
    // To be ready to send, either all messages should be received, or
    // all but the message requested out.
    return _number_edges_with_messages_received == number_edges() || (_number_edges_with_messages_received+1 == number_edges() && !_edge_received[edge_index]);
  }

  virtual bool ready_to_send_message_ab_initio(unsigned long edge_index) const {
    return false;
  }

  // Provides access for passing on new messages:
  Edge<VARIABLE_KEY>* get_edge_out(unsigned long edge_index) const {
    return _edges_out[edge_index];
  }

  unsigned long number_edges() const {
    // Equivalent to _edges_out.size():
    return _edges_in.size();
  }

  // Note: excludes ab initio messages.
  virtual bool can_potentially_pass_any_messages() const {
    return _number_edges_with_messages_received+1 >= number_edges();
  }

  bool edge_received(unsigned long edge_index) const {
    return _edge_received[edge_index];
  }

  virtual void print(std::ostream & os) const = 0;

  // Used to replace edges; Edge type has immutable source, dest
  // pointers and integer indices, so they must be replaced to edit
  // the graph. This is primarily used to merge hyperedges within
  // InferenceGraphBuilder.
  void rewire_edge(unsigned long edge_index, Edge<VARIABLE_KEY>*new_edge_in, Edge<VARIABLE_KEY>*new_edge_out) {
    Edge<VARIABLE_KEY>*edge_in = _edges_in[edge_index];
    Edge<VARIABLE_KEY>*edge_out = _edges_out[edge_index];

    _edges_in[edge_index] = new_edge_in;
    _edges_out[edge_index] = new_edge_out;

    if (edge_in->variables_ptr != new_edge_in->variables_ptr)
      delete edge_in->variables_ptr;
    delete edge_out;
    delete edge_in;
  }
};

template <typename VARIABLE_KEY>
std::ostream & operator << (std::ostream & os, const MessagePasser<VARIABLE_KEY> & rhs) {
  rhs.print(os);
  return os;
}

#endif
