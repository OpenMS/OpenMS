#ifndef _BELIEFPROPAGATIONINFERENCEENGINE_HPP
#define _BELIEFPROPAGATIONINFERENCEENGINE_HPP

#include "InferenceEngine.hpp"
#include "Scheduler.hpp"
#include "InferenceGraph.hpp"
#include "SetHash.hpp"
#include "Hyperedge.hpp"
#include "../Utility/to_string.hpp"
#include <unordered_map>

template <typename VARIABLE_KEY>
class BeliefPropagationInferenceEngine : public InferenceEngine<VARIABLE_KEY> {
protected:
  Scheduler<VARIABLE_KEY> & _scheduler;
  const InferenceGraph<VARIABLE_KEY> & _graph;

  bool every_nontrivial_edge_has_passed_at_least_one_message() const {
    bool res = true;
    for (MessagePasser<VARIABLE_KEY>*mp : _graph.message_passers)
      for (unsigned long k=0; k<mp->number_edges(); ++k) {
	Edge<VARIABLE_KEY>*edge = mp->get_edge_out(k);
	if (edge->source->number_edges() == 1 && static_cast<Hyperedge<VARIABLE_KEY>*>(edge->source) != NULL)
	  continue;
	if (edge->dest->number_edges() == 1 && static_cast<Hyperedge<VARIABLE_KEY>*>(edge->dest) != NULL)
	  continue;
	
	res = res && mp->edge_received(k);
      }
    return res;
  }

public:
  BeliefPropagationInferenceEngine(Scheduler<VARIABLE_KEY> & scheduler, const InferenceGraph<VARIABLE_KEY> & graph):
    _scheduler(scheduler),
    _graph(graph)
  { }

  std::vector<LabeledPMF<VARIABLE_KEY> > estimate_posteriors(const std::vector<std::vector<VARIABLE_KEY> > & joint_distributions_to_retrieve) {
    _scheduler.run_until_convergence();
    if ( ! every_nontrivial_edge_has_passed_at_least_one_message() )
      // This can happen if the graph is so large that not every edge
      // has been visited yet or if the graph contains a connected
      // component with no prior information:
      std::cerr << "Warning: Not every edge has passed a message (however posteriors may exist for the variables of interest). It may be that belief propagation hasn't yet converged (e.g., if this graph is large). If the graph is not large, check that your model doesn't add an edge using the wrong variable." << std::endl;

    std::vector<LabeledPMF<VARIABLE_KEY> > results;

    // Build a dictionary of varaibles to the message passers that
    // contain them. Set the initial number of bins with the number
    // of message passers to avoid resizing:
    std::unordered_map< std::unordered_set<VARIABLE_KEY>, const HUGINMessagePasser<VARIABLE_KEY>*, SetHash<VARIABLE_KEY> > variables_to_message_passers(_graph.message_passers.size());

    for (const MessagePasser<VARIABLE_KEY>* mp : _graph.message_passers) {
      const HUGINMessagePasser<VARIABLE_KEY>* hmp = dynamic_cast<const HUGINMessagePasser<VARIABLE_KEY>* >(mp);
      if (hmp != NULL) {
	// mp is a HUGINMessagePasser:

	const std::vector<VARIABLE_KEY> & ordered_variables = hmp->joint_posterior().ordered_variables();

	std::unordered_set<VARIABLE_KEY> unordered_variables(ordered_variables.begin(), ordered_variables.end());

	auto iter = variables_to_message_passers.find(unordered_variables);
	if ( iter == variables_to_message_passers.end() )
	  variables_to_message_passers[unordered_variables] = hmp;
      }
    }

    for (const std::vector<VARIABLE_KEY> & ordered_variables : joint_distributions_to_retrieve) {
      std::unordered_set<VARIABLE_KEY> unordered_variables(ordered_variables.begin(), ordered_variables.end());

      auto iter = variables_to_message_passers.find(unordered_variables);

      if (iter == variables_to_message_passers.end()) {
	std::string vars = "";
	for (const VARIABLE_KEY & var : unordered_variables)
	  vars += to_string(var) + " ";
	
	std::cerr << "Could not find posterior for variable set " << vars << std::endl;
	assert(false);
      }

      results.push_back(iter->second->joint_posterior().transposed(ordered_variables));
    }

    return results;
  }
};

#endif
