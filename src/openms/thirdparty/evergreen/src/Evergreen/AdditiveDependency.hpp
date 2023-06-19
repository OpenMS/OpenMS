#ifndef _ADDITIVEDEPENDENCY_HPP
#define _ADDITIVEDEPENDENCY_HPP

#include "../Engine/ConvolutionTreeMessagePasser.hpp"
#include "Dependency.hpp"

template <typename VARIABLE_KEY>
class AdditiveDependency : public Dependency<VARIABLE_KEY>, public PNormMixin {
protected:
  std::vector<std::vector<VARIABLE_KEY> > _inputs;
  std::vector<VARIABLE_KEY> _output;

public:
  AdditiveDependency(const std::vector<std::vector<VARIABLE_KEY> > & inputs, const std::vector<VARIABLE_KEY> & output, const double p_param):
    PNormMixin(p_param),
    _inputs(inputs),
    _output(output)
  {
    #ifdef ENGINE_CHECK
    for (const std::vector<VARIABLE_KEY> & in : _inputs)
      assert(in.size() == output.size() && "Dimension of each tuple in additive dependency must equal dimension of output tuple");
    #endif
  }

  virtual ConvolutionTreeMessagePasser<VARIABLE_KEY>* create_message_passer(InferenceGraphBuilder<VARIABLE_KEY> & igb) const {
    std::vector<ContextFreeMessagePasser<VARIABLE_KEY>*> hyperedges_in;
    std::vector<std::vector<VARIABLE_KEY>*> edge_labels_in;

    for (const std::vector<VARIABLE_KEY> & in : _inputs) {
      HUGINMessagePasser<VARIABLE_KEY>* hyperedge = igb.create_hyperedge();
      hyperedges_in.push_back(hyperedge);

      std::vector<VARIABLE_KEY>* edge_label_in = new std::vector<VARIABLE_KEY>(in);
      edge_labels_in.push_back(edge_label_in);
    }

    ContextFreeMessagePasser<VARIABLE_KEY>* hyperedge_out = igb.create_hyperedge();
    std::vector<VARIABLE_KEY>* edge_label_out = new std::vector<VARIABLE_KEY>(_output);
    
    // Note: The above allocations are not deallocated until
    // InferenceGraph destructor.
    return new ConvolutionTreeMessagePasser<VARIABLE_KEY>(hyperedges_in, edge_labels_in, hyperedge_out, edge_label_out, (unsigned char)_output.size(), this->p);
  }

  virtual std::vector<VARIABLE_KEY> get_all_variables_used() const {
    std::vector<VARIABLE_KEY> result = flatten(get_inputs());
    for (const VARIABLE_KEY & var : get_output())
      result.push_back(var);
    return result;
  }

  const std::vector<std::vector<VARIABLE_KEY> > & get_inputs() const {
    return _inputs;
  }
  
  const std::vector<VARIABLE_KEY> & get_output() const {
    return _output;
  }
};
#endif
