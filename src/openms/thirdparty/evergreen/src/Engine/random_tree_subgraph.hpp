#ifndef _RANDOM_TREE_SUBGRAPH_HPP
#define _RANDOM_TREE_SUBGRAPH_HPP

// Note: Assumes graph is connected (otherwise, split into connected
// compoonents and then call)
template <typename VARIABLE_KEY>
std::list<MessagePasser<VARIABLE_KEY>* > random_tree_subgraph(InferenceGraph<VARIABLE_KEY> & ig) {
  // Note: doing rand() % N can effectively discard some entropy; if
  // high-quality random numbers are necessary, then there are more
  // effective ways to do this.
  auto rand_int = [](unsigned long size) {
    return rand() % size;
  };

  // Clear node colors:
  for (unsigned long i=0; i<ig.message_passers.size(); ++i)
    ig.message_passers[i]->color = -1;

  // Choose random root:  
  MessagePasser<VARIABLE_KEY>*root = ig.message_passers[ rand_int(ig.message_passers.size()) ];

  // Build traversal order of tree via DFS:
  std::list<MessagePasser<VARIABLE_KEY>* > result;
  node_dfs({root}, [&result](MessagePasser<VARIABLE_KEY>*mp){
      result.push_back(mp);
      
      // Mark the node as visited:
      mp->color = 1;
    });

  return result;
}

#endif
