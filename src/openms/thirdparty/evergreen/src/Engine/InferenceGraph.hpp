#ifndef _INFERENCEGRAPH_HPP
#define _INFERENCEGRAPH_HPP

#include <unordered_set>
#include <list>
#include "MessagePasser.hpp"
#include "HUGINMessagePasser.hpp"
#include "../Utility/shuffled_sequence.hpp"

template <typename VARIABLE_KEY>
class InferenceGraph {
protected:
  void verify_all_connected_message_passers_included() {
    std::unordered_set<MessagePasser<VARIABLE_KEY>* > connected_mps(message_passers.begin(), message_passers.end());

    for (MessagePasser<VARIABLE_KEY>*mp : message_passers) {
      for (unsigned long edge_ind=0; edge_ind<mp->number_edges(); ++edge_ind) {
	Edge<VARIABLE_KEY>*edge = mp->get_edge_out(edge_ind);
	assert( connected_mps.find(edge->dest) != connected_mps.end() );
      }
    }
  }

  void verify_edges() {
    // Verify opposite edges:
    for (MessagePasser<VARIABLE_KEY>*mp : message_passers) {
      for (unsigned long edge_ind=0; edge_ind<mp->number_edges(); ++edge_ind) {
	Edge<VARIABLE_KEY>*edge = mp->get_edge_out(edge_ind);
	assert(edge->source == mp);
	assert(edge->source_edge_index == edge_ind);
	assert(edge->get_opposite_edge_ptr()->dest == mp);
      }
    }
  }

  void verify() {
    verify_all_connected_message_passers_included();
    verify_edges();
  }

public:
  // Permits modification of MessagePasser types via pointer, but not
  // modification of the pointers themselves:
  std::vector<MessagePasser<VARIABLE_KEY>* > message_passers;

  // Using rvalue references in constructors is efficient, but it also
  // pushes a bit for the InferenceGraph to own the underlying data
  // from here on out (InferenceGraph will delete allocated Edge and
  // MessagePasser types when it destructs).

  InferenceGraph(std::vector<MessagePasser<VARIABLE_KEY>* > && message_passers_param):
    message_passers(std::move(message_passers_param))
  {
    #ifdef ENGINE_CHECK
    // Only necessary on new construction (if constructing from
    // another InferenceGraph, that graph should have already been
    // verified).
    verify();
    #endif
  }

  InferenceGraph(InferenceGraph<VARIABLE_KEY> && ig):
    message_passers(std::move(ig.message_passers))
  { }

  // Disable copying so that destructor will not be called multiple times:
  InferenceGraph(const InferenceGraph<VARIABLE_KEY> &) = delete;

  // Disable copying so that destructor will not be called multiple times:
  const InferenceGraph & operator =(const InferenceGraph<VARIABLE_KEY> &) = delete;
  
  ~InferenceGraph() {
    // Delete _variables_ptr collections first so that edges are still
    // available. This is slightly tricky since these pointers are
    // shared between the forward and reverse edges, which could lead
    // to them being deleted multiple times. One solution would be to
    // set to NULL when deleted and only delete if not NULL, but since
    // that pointer is a *const, it cannot be assigned. Therefore,
    // implemented using a set for simplicity; can be implemented more
    // efficiently:
    std::unordered_set<const std::vector<VARIABLE_KEY>* > all_edge_labels;
    for (MessagePasser<VARIABLE_KEY>*mp : message_passers) {
      for (unsigned long k=0; k<mp->number_edges(); ++k) {
	const Edge<VARIABLE_KEY>*edge = mp->get_edge_out(k);
	all_edge_labels.insert(edge->variables_ptr);
      }
    }
    for (const std::vector<VARIABLE_KEY>*edge_label : all_edge_labels)
      delete edge_label;
  
    // Delete all edges out (ensures every edge will be deleted
    // exactly once):
    for (MessagePasser<VARIABLE_KEY>*mp : message_passers) {
      for (unsigned long k=0; k<mp->number_edges(); ++k) {
	Edge<VARIABLE_KEY>*edge = mp->get_edge_out(k);
	delete edge;
      }
    }

    // Delete message passers:
    for (MessagePasser<VARIABLE_KEY>*mp : message_passers)
      delete mp;
  }

  std::vector<Edge<VARIABLE_KEY>*> edges_ready_ab_initio() const {
    // Find edges that can pass from the first iteration (e.g.,
    // HUGINMessagePasser nodes with priors may sometimes start ready
    // to pass some of their edges):

    std::vector<Edge<VARIABLE_KEY>*> result;
  
    for (MessagePasser<VARIABLE_KEY>*mp : message_passers) {
      for (unsigned long edge_index=0; edge_index<mp->number_edges(); ++edge_index)
	if (mp->ready_to_send_message_ab_initio(edge_index)) {
	  Edge<VARIABLE_KEY>*edge = mp->get_edge_out(edge_index);

	  result.push_back(edge);
	}
    }
    return result;
  }

  // For debugging:
  void print(std::ostream & os) const {
    for (MessagePasser<VARIABLE_KEY>*mp : message_passers) {
      os << mp << " ";
      mp->print(os);
      os << std::endl;
      for (unsigned long k=0; k<mp->number_edges(); ++k) {
	Edge<VARIABLE_KEY>*edge = mp->get_edge_out(k);
	os << "\t";
	for (const VARIABLE_KEY & var : *edge->variables_ptr)
	  os << var << " ";
	os << edge->ready_to_pass() << " ";
	os << edge << ": ";
	os << edge->dest << " ";
	edge->dest->print(os);

	os << " received opposite on " << edge->get_opposite_edge_ptr() << " " << edge->source->edge_received( edge->get_opposite_edge_ptr()->dest_edge_index );
	os << std::endl;
      }
      os << std::endl;
    }
  }
};

// For applying depth and breadth first search on any lambda on
// MessagePasser<VARIABLE_KEY>*.

// Note: function is responsible for coloring edges. This is important
// to prevent infinite looping and infinite memory use (= crash).
template <typename VARIABLE_KEY, typename FUNCTION>
void node_dfs(std::list<MessagePasser<VARIABLE_KEY>* > queued_mps, FUNCTION function) {
  while (queued_mps.size() > 0) {
    MessagePasser<VARIABLE_KEY>*mp = queued_mps.front();
    queued_mps.pop_front();
    if (mp->color >= 0)
      continue;

    function(mp);

    // Visit the edges in random order:
    std::vector<unsigned long> shuffled_edge_indices = shuffled_sequence(mp->number_edges());
    for (unsigned long i : shuffled_edge_indices) {
      MessagePasser<VARIABLE_KEY>*next_mp = mp->get_edge_out(i)->dest;
      if (next_mp->color < 0)
	queued_mps.push_front(next_mp);
    }
  }
}

// To help the compiler wire an inlined list {a,b,...} to a std::list:
template <typename VARIABLE_KEY, typename FUNCTION>
void node_dfs(std::initializer_list<MessagePasser<VARIABLE_KEY>* > queued_mps_il, FUNCTION function) {
  node_dfs(std::list<MessagePasser<VARIABLE_KEY>* >(queued_mps_il), function);
}

template <typename VARIABLE_KEY, typename FUNCTION>
void node_bfs(std::list<MessagePasser<VARIABLE_KEY>* > queued_mps, FUNCTION function) {
  while (queued_mps.size() > 0) {
    MessagePasser<VARIABLE_KEY>*mp = queued_mps.front();
    queued_mps.pop_front();
    if (mp->color >= 0)
      continue;
    
    function(mp);
    
    // Visit the edges in random order:
    std::vector<unsigned long> shuffled_edge_indices = shuffled_sequence(mp->number_edges());
    for (unsigned long i : shuffled_edge_indices) {
      MessagePasser<VARIABLE_KEY>*next_mp = mp->get_edge_out(i)->dest;
      if (next_mp->color < 0)
	queued_mps.push_back(next_mp);
    }
  }
}

#include "split_connected_components.hpp"

#endif
