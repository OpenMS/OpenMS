#ifndef _HMMSCHEDULER_HPP
#define _HMMSCHEDULER_HPP

// A small, custom-made scheduler for HMMs. HMM should be constructed
// manually (without hyperedge types that would be produced by
// BetheGraphBuilder).

#include "../../Evergreen/evergreen.hpp"

template <typename VARIABLE_KEY>
class HMMScheduler : public FIFOScheduler<VARIABLE_KEY> {
public:
  HMMScheduler():
    // HMM graphs should have no loops, and hence dampening and
    // convergence threshold are moot. Likewise, use maximum unsigned
    // long as allowed number of iterations (convergence will occur
    // when no messages are woken).
    FIFOScheduler<VARIABLE_KEY>(0.0, 1e-6, -1ul)
  {}

  void add_ab_initio_edges(InferenceGraph<VARIABLE_KEY> & graph) {
    for (Edge<VARIABLE_KEY>* edge : graph.edges_ready_ab_initio()) {
      // Only allow ab initio edges coming from leaf nodes. Note that
      // this will not guaranteee all messages will be passed on
      // general graphs.

      bool source_is_leaf = edge->source->number_edges() == 1;
      if (source_is_leaf)
	this->_queue.push_if_not_in_queue(edge);
    }
  }

  // Note: An alternative approach would be to simply hard-code
  // message passing by overriding run_until_convergence.
};

#endif
