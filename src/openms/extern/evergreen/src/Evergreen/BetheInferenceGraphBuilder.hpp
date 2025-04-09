#ifndef _BETHEINFERENCEGRAPHBUILDER_HPP
#define _BETHEINFERENCEGRAPHBUILDER_HPP

#include "InferenceGraphBuilder.hpp"

// Note: BetheInferenceGraphBuilder is useful for medium-sized,
// densely connected graphs. For large chain graphs, simply build an
// HMM or the HMM-like tree decomposition (no automated builder for
// that yet). On those graphs, Bethe will perform poorly, and may even
// not visit all edges.

template <typename VARIABLE_KEY>
class BetheInferenceGraphBuilder : public InferenceGraphBuilder<VARIABLE_KEY> {
protected:
  void add_singleton_hyperedges_for_hugins() {
    // make a copy so that hyperedges can be inserted as you go:
    std::vector<MessagePasser<VARIABLE_KEY>* > mps = this->_message_passers;
    for (MessagePasser<VARIABLE_KEY>*mp : mps) {
      HUGINMessagePasser<VARIABLE_KEY>*hmp = dynamic_cast<HUGINMessagePasser<VARIABLE_KEY>* >(mp);
      if (hmp != NULL) {
	for (const VARIABLE_KEY & var : hmp->joint_posterior().ordered_variables()) {
	  Hyperedge<VARIABLE_KEY>*he = this->create_hyperedge();
	  hmp->bind_to(he, new std::vector<VARIABLE_KEY>{var});
	}
      }
    }
  }

  void merge_hyperedges_with_identical_incident_variable_sets() {
    // Note that the inner unordered_set could be changed to a vector
    // for speedup; you would simply check if the last element in the
    // vector was he, and if not, push_back. This exploits the fact
    // that each hyperedge is considered one at a time, and so
    // multiple insertions can only happen when it is inserted back to
    // back.
    std::unordered_map<std::unordered_set<VARIABLE_KEY>, std::unordered_set<Hyperedge<VARIABLE_KEY>* >, SetHash<VARIABLE_KEY> > var_sets_to_hyperedges;
    for (MessagePasser<VARIABLE_KEY>*mp : this->_message_passers) {
      Hyperedge<VARIABLE_KEY>*he = dynamic_cast<Hyperedge<VARIABLE_KEY>* >(mp);
      if (he != NULL) {
	std::unordered_set<VARIABLE_KEY> vars_used = he->variables_used_by_incident_edges();
	auto iter = var_sets_to_hyperedges.find(vars_used);
	if (iter == var_sets_to_hyperedges.end())
	  var_sets_to_hyperedges[vars_used] = std::unordered_set<Hyperedge<VARIABLE_KEY>* >();
	var_sets_to_hyperedges[vars_used].insert(he);
      }
    }
    std::vector<std::vector<Hyperedge<VARIABLE_KEY>* > > hes_to_merge;
    for (const std::pair<std::unordered_set<VARIABLE_KEY>, std::unordered_set<Hyperedge<VARIABLE_KEY>* > > & vars_and_hes : var_sets_to_hyperedges) {
      const std::unordered_set<Hyperedge<VARIABLE_KEY>* > & he_set = vars_and_hes.second;
      std::vector<Hyperedge<VARIABLE_KEY>* > collection_to_merge(he_set.begin(), he_set.end());
      hes_to_merge.push_back(collection_to_merge);
    }
    this->merge_hyperedges(hes_to_merge);
  }

  void bind_singletons_to_superset_hyperedges() {
    std::unordered_map<VARIABLE_KEY, std::unordered_set<Hyperedge<VARIABLE_KEY>* > > vars_to_hyperedges;
    for (MessagePasser<VARIABLE_KEY>*mp : this->_message_passers) {
      Hyperedge<VARIABLE_KEY>*he = dynamic_cast<Hyperedge<VARIABLE_KEY>* >(mp);
      if (he != NULL) {
	for (const VARIABLE_KEY & var : he->variables_used_by_incident_edges()) {
	  auto iter = vars_to_hyperedges.find(var);
	  if (iter == vars_to_hyperedges.end())
	    vars_to_hyperedges[var] = std::unordered_set<Hyperedge<VARIABLE_KEY>* >();
	  
	  vars_to_hyperedges[var].insert(he);
	}
      }
    }

    for (const std::pair<VARIABLE_KEY, std::unordered_set<Hyperedge<VARIABLE_KEY>* > > & var_and_he_set : vars_to_hyperedges) {
      const VARIABLE_KEY & var = var_and_he_set.first;
      const std::unordered_set<Hyperedge<VARIABLE_KEY>* > & he_set = var_and_he_set.second;

      // singleton_he should be the only surviving hyperedge that
      // contains the singleton set of this variable (otherwise it
      // would have been merged above):
      bool singleton_found = false;
      Hyperedge<VARIABLE_KEY>*singleton_he = NULL;
      for (Hyperedge<VARIABLE_KEY>*singleton_he_local : he_set)
	if (singleton_he_local->variables_used_by_incident_edges().size() == 1) {
	  singleton_he = singleton_he_local;
	  singleton_found = true;
	  break;
	}

      // Hyperedge with singleton is not found, but there are
      // higher-order hyperedges using that variable. E.g., one
      // Dependency could create {A,B,C} hyperedge and another could
      // create {A,B} hyperedge.
      
      // Therefore, if {A} and {B} hyperedges don't exist, then they
      // should be created so that {A,B,C} and {A,B} can be connected
      // through them.
      if (! singleton_found && he_set.size() > 1)
	singleton_he = this->create_hyperedge();

      if (singleton_he == NULL)
	continue;

      // for all he in he_set, the ones with >1 variables should be
      // bound to the canonical singleton for that variable
      for (Hyperedge<VARIABLE_KEY>*he : he_set)
	if (he->variables_used_by_incident_edges().size() > 1) {
	  singleton_he->bind_to(he, new std::vector<VARIABLE_KEY>{var});
	}
    }
  }
  
  void construct_graph_connections() {
  // 1. Add singleton hyperedges for all vars in HUGINMessagePasser types with priors:
    add_singleton_hyperedges_for_hugins();

    // 2. Merge all hyperedges that use identical variable sets (these
    // hyperedges are either singletons constructed above or were
    // constructed by the Dependency types; in either case, they are
    // atomic and cannot be broken down).
   merge_hyperedges_with_identical_incident_variable_sets();

    // 3. Bind the singleton hyperedges to higher-order hyperedges
    // that are supersets:
    bind_singletons_to_superset_hyperedges();
  }
};

#endif
