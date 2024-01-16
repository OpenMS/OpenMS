#ifndef _KIKUCHIGRAPH_HPP
#define _KIKUCHIGRAPH_HPP

template <typename VARIABLE_KEY>
class KikuchiGraph {
public:
  // Permits modification of MessagePasser types via pointer, but not
  // modification of the pointers themselves:
  const std::vector<MessagePasser<VARIABLE_KEY>* > _message_passers;

  KikuchiGraph(std::vector<MessagePasser<VARIABLE_KEY>* > && message_passers):
    _message_passers(std::move(message_passers))
  { }

  ~KikuchiGraph() {
    // Delete _variables_ptr collections first so that edges are still
    // available (to get opposite):
    for (MessagePasser<VARIABLE_KEY>*mp : _message_passers) {
      for (unsigned long k=0; k<mp->number_edges(); ++k) {
	Edge<VARIABLE_KEY>*edge = mp->get_edge_out(k);
	if (edge->variables_ptr != NULL) {
	  delete edge->variables_ptr;
	  edge->variables_ptr = NULL;
	  edge->opposite()->variables_ptr = NULL;
	}
      }
    }

    // Delete all edges out (ensures every edge will be deleted
    // exactly once):
    for (MessagePasser<VARIABLE_KEY>*mp : _message_passers) {
      for (unsigned long k=0; k<mp->number_edges(); ++k) {
	Edge<VARIABLE_KEY>*edge = mp->get_edge_out(k);
	delete edge;
      }
    }

    // Delete message passers:
    for (MessagePasser<VARIABLE_KEY>*mp : _message_passers) {
      delete mp;
  }
};

#endif
