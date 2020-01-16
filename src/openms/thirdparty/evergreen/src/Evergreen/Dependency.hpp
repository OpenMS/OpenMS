#ifndef _DEPENDENCY_HPP
#define _DEPENDENCY_HPP

template <typename VARIABLE_KEY>
class InferenceGraphBuilder;

template <typename VARIABLE_KEY>
class Dependency {
public:
  // The idea here is that dependencies create message passers and
  // bind them to context free message passers through which they
  // reach the other message passers in the graph. This allows context
  // to be respected, but simplifies the job for building a graph
  // (which can then be performed solely by linking context free
  // message passers to each other where necessary). Thus, in derived
  // classes, create_message_passer should create a message passer of
  // the appropriate type and create and bind it to the necessary
  // hyperedges.
  virtual MessagePasser<VARIABLE_KEY>* create_message_passer(InferenceGraphBuilder<VARIABLE_KEY> & igb) const = 0;
  
  //virtual std::vector<VARIABLE_KEY> get_all_variables_used() const = 0;
  
  virtual ~Dependency() {}
};

#endif
