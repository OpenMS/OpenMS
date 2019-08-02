#ifndef _INFERENCEGRAPHBUILDER_HPP
#define _INFERENCEGRAPHBUILDER_HPP

#include "../Engine/InferenceGraph.hpp"
#include "../Engine/Hyperedge.hpp"
#include "../Engine/SetHash.hpp"
#include "Dependency.hpp"

template <typename VARIABLE_KEY>
class InferenceGraphBuilder {
private:
   bool _has_created_graph;

protected:
  std::vector<MessagePasser<VARIABLE_KEY>* > _message_passers;
  
  // To be defined by derived types; e.g., Bethe construction defines
  // one connection scheme, there are multiple schemes for Kikuchi
  // construction, etc.
  virtual void construct_graph_connections() = 0;

public:
  InferenceGraphBuilder():
    _has_created_graph(false)
  { }

  virtual ~InferenceGraphBuilder() {
    if ( ! _has_created_graph ) {
      // TODO: If no graph has been created, then the allocated memory
      // needs to be destroyed.
      assert(false && "InferenceGraphBuilder needs to create a graph or else it leaks memory");
    }
  }

  Hyperedge<VARIABLE_KEY>* create_hyperedge() {
    Hyperedge<VARIABLE_KEY>* hyperedge = new Hyperedge<VARIABLE_KEY>();
    _message_passers.push_back(hyperedge);

    return hyperedge;
  }

  // Arguments are a collection of vectors of hyperedges to merge. 
  void merge_hyperedges(const std::vector<std::vector<Hyperedge<VARIABLE_KEY>* > > & hes_to_merge) {
    std::vector<MessagePasser<VARIABLE_KEY>* > new_message_passers;
    // add all non-hyperedge MessagePasser types:
    for (MessagePasser<VARIABLE_KEY>*mp : _message_passers) {
      Hyperedge<VARIABLE_KEY>* he = dynamic_cast<Hyperedge<VARIABLE_KEY>* >(mp);
      if (he == NULL)
	new_message_passers.push_back(mp);
    }

    // Use the first Hyperedge in the collection of each as the
    // Hyperedge that will be kept (others will be merged into it)
    for (const std::vector<Hyperedge<VARIABLE_KEY>* > & he_vec : hes_to_merge) {
      Hyperedge<VARIABLE_KEY>* he_to_keep = he_vec[0];
      new_message_passers.push_back(he_to_keep);

      for (unsigned long i=1; i<he_vec.size(); ++i)
	he_to_keep->absorb_hyperedge(he_vec[i]);
    }

    _message_passers = new_message_passers;
  }

  const std::vector<MessagePasser<VARIABLE_KEY>* > & message_passers() const {
    return _message_passers;
  }

  void insert_dependency(const Dependency<VARIABLE_KEY> & dep) {
    MessagePasser<VARIABLE_KEY>*mp = dep.create_message_passer(*this);
    _message_passers.push_back(mp);
  }

  InferenceGraph<VARIABLE_KEY> to_graph() {
    // To prevent two InferenceGraph destructors from deleting the
    // same memory multiple times:
    assert(! _has_created_graph && "Each InferenceGraphBuilder should only be used to create a single graph; for a copy of the graph, it should be copied");
    
    construct_graph_connections();
    _has_created_graph = true;
    
    return InferenceGraph<VARIABLE_KEY>(std::move(_message_passers));
  }
};

#endif
