#ifndef _RANDOMSUBTREESCHEDULER_HPP
#define _RANDOMSUBTREESCHEDULER_HPP

#include "Scheduler.hpp"
#include "random_tree_subgraph.hpp"

template <typename VARIABLE_KEY>
class RandomSubtreeScheduler : public Scheduler<VARIABLE_KEY> {
private:
  std::list<MessagePasser<VARIABLE_KEY>* > _mp_ordering_1, _mp_ordering_2;
  std::list<MessagePasser<VARIABLE_KEY>* >* _current_mp_ordering;
  bool _any_passed_this_batch;

  // Returns true if any non-convergent messages could be passed,
  // false otherwise:
  bool pass_all_messages_possible(MessagePasser<VARIABLE_KEY>*mp) {
    bool any_passed = false;

    // Note: Could save some time by ignoring cases where it's clear
    // no message can be passed:
    // if (mp->can_potentially_pass_any_messages())
    // However, this would be a little tricky since it does not
    // include the ab initio case.
    for (unsigned long i=0; i<mp->number_edges(); ++i) {
      if (mp->ready_to_send_message_ab_initio(i) || mp->ready_to_send_message(i)) {
	Edge<VARIABLE_KEY>*edge = mp->get_edge_out(i);
	
	LabeledPMF<VARIABLE_KEY> new_msg = mp->update_and_get_message_out(i);

	if ( ! edge->has_message() || (edge->has_message() && mse_divergence(edge->get_possibly_outdated_message(), new_msg) > this->_convergence_threshold) ) {
	  any_passed = true;

	  if (edge->has_message())
	    // Dampen:
	    new_msg = dampen(edge->get_possibly_outdated_message(), new_msg, this->_dampening_lambda).transposed(*edge->variables_ptr);
	
	  edge->set_message( std::move(new_msg) );
	
	  MessagePasser<VARIABLE_KEY>*dest_mp = edge->dest;
	  dest_mp->receive_message_in_and_update(edge->dest_edge_index);
	}
      }
    }

    return any_passed;
  }
public:
  RandomSubtreeScheduler(double dampening_lambda_param, double convergence_threshold_param, unsigned long maximum_iterations_param):
    Scheduler<VARIABLE_KEY>(dampening_lambda_param, convergence_threshold_param, maximum_iterations_param),
    _current_mp_ordering(NULL),
    _any_passed_this_batch(true)
  { }

  void add_ab_initio_edges(InferenceGraph<VARIABLE_KEY> & ig) {
    _mp_ordering_1 = random_tree_subgraph(ig);
    _mp_ordering_2 = random_tree_subgraph(ig);

    _current_mp_ordering = &_mp_ordering_1;
  }

  unsigned long process_next_edges() {
    unsigned long iteration = 0;
    _any_passed_this_batch = false;

    // Gather messages in:
    for (auto iter = _current_mp_ordering->rbegin(); iter != _current_mp_ordering->rend() && iteration < this->_maximum_iterations; ++iter, ++iteration) {
      bool iter_passes = pass_all_messages_possible(*iter);
      _any_passed_this_batch = _any_passed_this_batch || iter_passes;
    }
    // Scatter messages out:
    for (auto iter = _current_mp_ordering->begin(); iter != _current_mp_ordering->end() && iteration < this->_maximum_iterations; ++iter, ++iteration) {
      bool iter_passes = pass_all_messages_possible(*iter);
      _any_passed_this_batch = _any_passed_this_batch || iter_passes;
    }

    // Oscillate the current ordering between the two trees:
    if (_current_mp_ordering == &_mp_ordering_1)
      _current_mp_ordering = &_mp_ordering_2;
    else
      _current_mp_ordering = &_mp_ordering_1;

    return iteration;
  }

  bool has_converged() const {
    return ! _any_passed_this_batch;
  }
};

#endif
