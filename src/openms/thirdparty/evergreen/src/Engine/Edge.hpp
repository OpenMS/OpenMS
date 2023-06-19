#ifndef _EDGE_HPP
#define _EDGE_HPP

#include "Queueable.hpp"
#include "../PMF/LabeledPMF.hpp"

template <typename VARIABLE_KEY>
class MessagePasser;

// Note: this is currently hard coded to use MSE-based divergence; it
// could be made more general later.
template <typename VARIABLE_KEY>
class Edge : public Queueable {
public:
  MessagePasser<VARIABLE_KEY> *const source, *const dest;
  const unsigned long source_edge_index, dest_edge_index;

  const std::vector<VARIABLE_KEY> *const variables_ptr;

  long color;

protected:
  bool _up_to_date;

  // Store current and previous message for use with dampening and
  // with divergence-based priority (i.e., edges with most changed
  // messages updated first among those ready to pass).
  LabeledPMF<VARIABLE_KEY> _current_message;

public:

  Edge(MessagePasser<VARIABLE_KEY>*source_param, MessagePasser<VARIABLE_KEY>*dest_param, const std::vector<VARIABLE_KEY>*variables_ptr_param, unsigned long source_edge_index_param, unsigned long dest_edge_index_param):
    source(source_param),
    dest(dest_param),
    source_edge_index(source_edge_index_param),
    dest_edge_index(dest_edge_index_param),
    variables_ptr(variables_ptr_param),
    color(0),
    _up_to_date(false)
  { }

  void set_message(LabeledPMF<VARIABLE_KEY> && msg) {
    // To prevent exponential feedback:
    msg.reset_log_normalization_constant();
    
    _current_message = std::move(msg);
    _up_to_date = true;
  }

  Edge*get_opposite_edge_ptr() const {
    return dest->get_edge_out(dest_edge_index);
  }

  const LabeledPMF<VARIABLE_KEY> & get_message() const {
    #ifdef ENGINE_CHECK
    assert( ready_to_pass() );
    #endif

    return _current_message;
  }

  void reset_message_norm_constant() {
      _current_message.reset_log_normalization_constant();
  }

  // Does not require ready_to_pass(), only has_message():
  const LabeledPMF<VARIABLE_KEY> & get_possibly_outdated_message() const {
    #ifdef ENGINE_CHECK
    assert( has_message() );
    #endif

    return _current_message;
  }

  void set_not_up_to_date() {
    _up_to_date = false;
  }

  bool up_to_date() const {
    return _up_to_date;
  }

  bool has_message() const {
    return _current_message.dimension() > 0;
  }

  bool ready_to_pass() const {
    return has_message() && _up_to_date;
  }
};

#endif
