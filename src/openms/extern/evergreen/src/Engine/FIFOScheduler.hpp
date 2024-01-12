#ifndef _FIFOSCHEDULER_HPP
#define _FIFOSCHEDULER_HPP

#include "ListQueue.hpp"

template <typename VARIABLE_KEY>
std::ostream & operator <<(std::ostream & os, const std::vector<VARIABLE_KEY> & rhs) {
  os << "[ ";
  for (const VARIABLE_KEY & var : rhs)
    os << var << " ";
  os << "]";
  return os;
}

template <typename VARIABLE_KEY>
class FIFOScheduler : public Scheduler<VARIABLE_KEY> {
protected:
  ListQueue<VARIABLE_KEY> _queue;

public:
  FIFOScheduler(double dampening_lambda, double convergence_threshold, unsigned long maximum_iterations):
    Scheduler<VARIABLE_KEY>(dampening_lambda, convergence_threshold, maximum_iterations)
  {}

  void add_ab_initio_edges(InferenceGraph<VARIABLE_KEY> & graph){
    // todo: shuffles ab initio edges (could do them in DFS/BFS formation for greater efficiency)
    std::vector<Edge<VARIABLE_KEY>*> starters;
    for (Edge<VARIABLE_KEY>* edge : graph.edges_ready_ab_initio())
      starters.push_back(edge);

    // shuffle:
    for (unsigned int i=0; i<starters.size(); ++i) {
      int j = rand()%starters.size();
      std::swap(starters[i], starters[j]);
    }

    for (Edge<VARIABLE_KEY>* edge : starters)
      _queue.push_if_not_in_queue(edge);
  }

  unsigned long process_next_edges() {
    if ( _queue.is_empty() )
      return 0;

    Edge<VARIABLE_KEY>*edge = _queue.pop_next();

    MessagePasser<VARIABLE_KEY>*source_mp = edge->source;
    // Update the message in the edge immediately before use (in a
    // lazy manner):

    LabeledPMF<VARIABLE_KEY> new_msg = source_mp->update_and_get_message_out(edge->source_edge_index);
    if ( ! edge->has_message() || (edge->has_message() && mse_divergence(edge->get_possibly_outdated_message(), new_msg) > this->_convergence_threshold) ) {
      if (edge->has_message())
	// Dampen:
	new_msg = dampen(edge->get_possibly_outdated_message(), new_msg, this->_dampening_lambda).transposed(*edge->variables_ptr);

      edge->set_message( std::move(new_msg) );

      // Receive the message:
      MessagePasser<VARIABLE_KEY>*dest_mp = edge->dest;
      dest_mp->receive_message_in_and_update(edge->dest_edge_index);

      // Wake up other edges:

      // Do not bother trying to wake any edges if <n-1 messages have
      // been received by dest_mp:

      if (dest_mp->can_potentially_pass_any_messages()) {
	unsigned long edge_index_received = edge->dest_edge_index;
	for (unsigned long edge_index_out=0; edge_index_out<dest_mp->number_edges(); ++edge_index_out) {
	  // Do not wake edge opposite to the edge received:
	  if (edge_index_out != edge_index_received && dest_mp->ready_to_send_message(edge_index_out)) {
	    Edge<VARIABLE_KEY>*e = dest_mp->get_edge_out(edge_index_out);

	    _queue.push_if_not_in_queue(e);
	  }
	}
      }
    }
    return 1;
  }
  bool has_converged() const {
    return _queue.is_empty();
  }
};

#endif
