#ifndef _PRIORITYSCHEDULER_HPP
#define _PRIORITYSCHEDULER_HPP

#include "Scheduler.hpp"
#include "SetQueue.hpp"

// Note: for some graphs (e.g., trees and HMM-like graphs), the
// runtime may be in O(n log(n)) instead of O(n), because of
// SetQueue. These graphs will be solved in O(n) with
// FIFOScheduler.
template <typename VARIABLE_KEY>
class PriorityScheduler : public Scheduler<VARIABLE_KEY> {
protected:
  SetQueue<VARIABLE_KEY> _queue;

  void set_priority_without_updating_message_and_update_queue(Edge<VARIABLE_KEY>*e, double new_priority) {
    // If message is not up to date, it will be refreshed when
    // passing.

    if ( ! e->in_queue ) {
      // If the edge is not in the queue, only add it to the queue
      // (which will update its priority) as long as convergence has
      // not been reached:
      if (new_priority > this->_convergence_threshold) {
	// This edge has changed more than the convergence criteria
	// allows; add it to the queue.
	
	_queue.push_or_update(e, new_priority);
      }
    }
  }

  void set_message_at_edge_and_update_queue(Edge<VARIABLE_KEY>*e, LabeledPMF<VARIABLE_KEY> && msg, double priority_bias=0.0) {
    double new_priority;
    if (e->has_message()) {
      // Transpose to guarantee that the message gets correct variable
      // order (this is important for some context-dependent message
      // passers, e.g. a multidimensional ConvolutionTreeMessagePasser):

      new_priority = mse_divergence(e->get_possibly_outdated_message(), msg);
      // Note that since the edge has been awoken, it will be out of
      // date (that is what this function is addressing); therefore,
      // call get_possibly_outdated_message(), because that does not
      // enforce check of whether or not it is up to date.

      msg = dampen(e->get_possibly_outdated_message(), msg, this->_dampening_lambda).transposed(*e->variables_ptr);
    }
    else {
      // Otherwise rank with sparsest messages first:
      const Tensor<double> & tab = msg.pmf().table();

      #ifdef SHAPE_CHECK
      assert( tab.flat_size() > 0 );
      #endif

      // When priority_bias > 1.0, ensures priority > 1 >= max MSE,
      // which means this sparsity score will always trump the
      // divergence score (and therefore edges with no previous
      // message will always be prioritized earlier than messages with
      // a previous message, regardless of the respective
      // sparsity-based priority and MSE).

      // Therefore, use priority_bias=2.0 for initial edges to
      // hyperdges and priority_bias=1.0 for initial edges back from
      // hyperedges. This will ensure ab initio messages are first
      // passed to hyperedges, then edges back from hyperedges, then
      // edges woken up.
      new_priority = priority_bias + 1.0 / tab.flat_size();
    }

    if ( ! e->in_queue ) {
      // If the edge is not in the queue, only add it to the queue
      // (which will update its priority) as long as convergence has
      // not been reached:
      
      if (new_priority >= this->_convergence_threshold)
	// This edge has changed more than the convergence criteria
	// allows; add it to the queue.
	
	_queue.push_or_update(e, new_priority);
    }
    else {
      // If the edge is in the queue, it has not yet passed the old
      // message. Therefore, even if the change between the old
      // message and the new message is very small, neither have been
      // passed, and so convergence is not necessarily reached. Thus,
      // only allow the edge to move forward in the queue, but not
      // backward.
      if (new_priority > e->priority)
	_queue.push_or_update(e, new_priority);
    }

    e->set_message(std::move(msg));
  }

public:
  PriorityScheduler(double dampening_lambda, double convergence_threshold, unsigned long maximum_iterations):
    Scheduler<VARIABLE_KEY>(dampening_lambda, convergence_threshold, maximum_iterations)
  {}
  
  void add_ab_initio_edges(InferenceGraph<VARIABLE_KEY> & graph){
    for (Edge<VARIABLE_KEY>* edge : graph.edges_ready_ab_initio())
      set_priority_without_updating_message_and_update_queue(edge, 2.0);
  }

  unsigned long process_next_edges() {
    if ( _queue.is_empty() )
      return 0;

    Edge<VARIABLE_KEY>*edge = _queue.pop_max();

    // If this edge was enqueued lazily (i.e., if the message has not
    // been set) or if the edge is not up to date, set its message
    // now:
    MessagePasser<VARIABLE_KEY>*source_mp = edge->source;
    if ( ! edge->ready_to_pass() ) {
      edge->set_message( std::move(source_mp->update_and_get_message_out(edge->source_edge_index)) );
    }

    MessagePasser<VARIABLE_KEY>*dest_mp = edge->dest;

    #ifdef PRINT_MESSAGES
    std::cout << "Message Passed: " << std::endl;
    std::cout << "FROM  ";
    edge->source->print(std::cout);
    std::cout << "  TO  ";
    edge->dest->print(std::cout);
    std::cout << "  WITH  " << edge->get_message() << std::endl;
    #endif

    dest_mp->receive_message_in_and_update(edge->dest_edge_index);

    // Iterate through the outgoing edges other than the one just
    // received:

    // Relies on the fact that edges must be constructed symmetrically
    // (i.e., input and output edges must be added simultaneously, so
    // for any MessagePasser, edge->dest_edge_index ==
    // edge->get_opposite_edge_ptr()->source_edge_index).
    unsigned long edge_index_received = edge->dest_edge_index;
    for (unsigned long edge_index_out=0; edge_index_out<dest_mp->number_edges(); ++edge_index_out) {
      // Do not wake edge opposite to the edge received:
      if (edge_index_out != edge_index_received && dest_mp->ready_to_send_message(edge_index_out)) {
	Edge<VARIABLE_KEY>*e = dest_mp->get_edge_out(edge_index_out);

	set_message_at_edge_and_update_queue(e, dest_mp->update_and_get_message_out(edge_index_out));
      }
    }
    return 1;
  }
  bool has_converged() const {
    // Edges will be added to the queue in a non-lazy fashion;
    // therefore, all edges with messages that are not converged
    // should be in the queue.
    return _queue.is_empty();
  }
};

#endif
