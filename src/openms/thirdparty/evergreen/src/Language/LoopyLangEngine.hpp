#ifndef _LOOPYLANGENGINE_HPP
#define _LOOPYLANGENGINE_HPP

#include "LangEngine.hpp"

typedef struct LoopyLangEngineT : public LangEngineT {
  void build(std::vector<std::vector<Dependency<std::string>* > > dependencies_of_subgraphs, double dampening, double epsilon, long max_iter) {
    for(InferenceGraph<std::string>* ptr: ig_ptrs)
      delete(ptr);
    ig_ptrs.clear();

    for(Scheduler<std::string>* ptr: sched_ptrs)
      delete(ptr);

    for(InferenceEngine<std::string>* ptr : eng_ptrs)
      delete(ptr);

    eng_ptrs.resize(dependencies_of_subgraphs.size());
    sched_ptrs.resize(dependencies_of_subgraphs.size());
    ig_ptrs.resize(dependencies_of_subgraphs.size());

    for (int i = 0; i < dependencies_of_subgraphs.size(); ++i) {
      const std::vector<Dependency<std::string>* > & deps = dependencies_of_subgraphs[i];
      BetheInferenceGraphBuilder<std::string> igb;
      for (Dependency<std::string>* dep : deps)
	igb.insert_dependency(*dep);
      sched_ptrs[i] = new FIFOScheduler<std::string>(dampening, epsilon, max_iter);

      ig_ptrs[i] = new InferenceGraph<std::string>(igb.to_graph());
      sched_ptrs[i]->add_ab_initio_edges(*ig_ptrs[i]);
      eng_ptrs[i] = new BeliefPropagationInferenceEngine<std::string>(*sched_ptrs[i], *ig_ptrs[i]);
    }

    is_built = true;
  }
} LoopyLangEngine;

#endif
