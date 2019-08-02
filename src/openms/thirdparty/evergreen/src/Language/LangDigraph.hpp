#ifndef _LANGGRAPH_HPP
#define _LANGGRAPH_HPP
#include <unordered_map>
#include <map>
#include "FrozenSet.hpp"

template <typename NODE>
class LangGraph {
private:
  std::unordered_map<NODE, FrozenSet<NODE> > node_to_edges;
public:
  
  void insert_edge(const NODE & u, const std::set<NODE> & v) {
    FrozenSet<std::string> set(v);
    node_to_edges.insert(std::pair<NODE, FrozenSet<NODE> >(u, set));
  }

  void insert_clique(const FrozenSet<NODE> & fs);

  std::set<NODE>dfs(const NODE & u, std::set<NODE> & connected_component, std::unordered_map<NODE, bool> & node_is_connected) const {
    
    if(node_is_connected[u]) {
      return connected_component;
    }
    connected_component.insert(u);
    node_is_connected[u] = true;
    auto adj_frozenset = node_to_edges.find(u);
    if (adj_frozenset != node_to_edges.end()) {
      std::set<NODE> adj_nodes = adj_frozenset->second.get_set();
      for (NODE v : adj_nodes) {
        dfs(v, connected_component, node_is_connected);        
      }  
    }
    return connected_component;
  } 

  std::vector<FrozenSet<NODE> > get_connected_subgraphs() const {    
    std::vector<FrozenSet<NODE> > connected_components;
    std::unordered_map<NODE, bool> node_is_connected;
    for(auto node_iter = node_to_edges.begin(); node_iter != node_to_edges.end(); ++node_iter)
      node_is_connected[node_iter->first] = false;
    for(auto node_iter = node_is_connected.begin(); node_iter != node_is_connected.end(); ++node_iter) {
      if (!node_iter->second) {
        std::set<NODE> connected_component;
        dfs(node_iter->first, connected_component, node_is_connected);
        FrozenSet<std::string> frozen_connected_component(connected_component);
        connected_components.push_back(frozen_connected_component);
      }     
    }
    return connected_components;
  }
};

#endif
