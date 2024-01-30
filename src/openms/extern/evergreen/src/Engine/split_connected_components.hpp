#ifndef _SPLIT_CONNECTED_COMPONENTS_HPP
#define _SPLIT_CONNECTED_COMPONENTS_HPP

// This can be done easily in O(n log(n)); here it is done in a
// slightly less straightforward way to get an O(n) runtime.
template <typename VARIABLE_KEY>
std::vector<InferenceGraph<VARIABLE_KEY> > split_connected_components(InferenceGraph<VARIABLE_KEY> && ig) {
  // Clear colors:
  for (unsigned long i=0; i<ig.message_passers.size(); ++i)
    ig.message_passers[i]->color = -1L;

  // Assign colors for connected components in O(n):
  unsigned long current_color = 0;
  for (unsigned long i=0; i<ig.message_passers.size(); ++i) {
    MessagePasser<VARIABLE_KEY>*mp = ig.message_passers[i];
    if (mp->color < 0) {

      // Note: could possibly simplify, queueing the nodes belonging
      // to current_color inside this closure:
      node_dfs({mp}, [current_color](MessagePasser<VARIABLE_KEY>*mp){
	  mp->color = current_color;
	});

      ++current_color;
    }
  }

  // Group nodes by color:
  std::vector<std::vector<MessagePasser<VARIABLE_KEY>* > > mps_grouped_by_color(current_color);
  for (unsigned long i=0; i<ig.message_passers.size(); ++i) {
    MessagePasser<VARIABLE_KEY>*mp = ig.message_passers[i];
    mps_grouped_by_color[mp->color].push_back(mp);
  }
      
  // Build several graphs for result:
  std::vector<InferenceGraph<VARIABLE_KEY> > result;
  for (std::vector<MessagePasser<VARIABLE_KEY>* > & mps : mps_grouped_by_color)
    result.push_back( std::move(InferenceGraph<VARIABLE_KEY>(std::move(mps))) );

  // Prevent destructor from being called multiple times:
  ig.message_passers = std::vector<MessagePasser<VARIABLE_KEY>* >();
  return result;
}

#endif
