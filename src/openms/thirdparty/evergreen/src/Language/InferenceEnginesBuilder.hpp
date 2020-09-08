#ifndef _INFERENCEENGINESBUILDER_HPP
#define _INFERENCEENGINESBUILDER_HPP

#include <string>
#include "../Evergreen/evergreen.hpp"

class InferenceEnginesBuilder {
public:
  virtual std::vector<InferenceEngine<std::string>* > build_engines(const std::vector<std::vector<Dependency<std::string>*> > & deps) = 0;
  virtual ~InferenceEnginesBuilder() {}
};

class BruteForceInferenceEnginesBuilder : public InferenceEnginesBuilder {
public:
  std::vector<InferenceEngine<std::string>* > build_engines(const std::vector<std::vector<Dependency<std::string>* > > & dependencies_of_subgraphs) {
    std::vector<InferenceEngine<std::string>* > result(dependencies_of_subgraphs.size());

    for (unsigned long i=0; i<dependencies_of_subgraphs.size(); ++i) {
      std::vector<TableDependency<std::string> > table_deps;
      std::vector<AdditiveDependency<std::string> > additive_deps;

      const std::vector<Dependency<std::string>* > & dependency_subgraph = dependencies_of_subgraphs[i];
      for (Dependency<std::string>* dep : dependency_subgraph) {
	if (dynamic_cast<TableDependency<std::string>*>(dep) != NULL) {
	  // TableDependency type
	  TableDependency<std::string>*table_dep = dynamic_cast<TableDependency<std::string>*>(dep);  
	  table_deps.push_back(*table_dep); 
	}
	else if (dynamic_cast<AdditiveDependency<std::string>*>(dep) != NULL) {
	  // AdditiveDependency type
	  AdditiveDependency<std::string>*additive_dep = dynamic_cast<AdditiveDependency<std::string>*>(dep);
	  additive_deps.push_back(*additive_dep);
	}
	else {
	  // error: user needs to define new dependency type in this if-else ladder
	}
      }
      result[i] = new BruteForceInferenceEngine<std::string>(table_deps, additive_deps, p);
    }

    return result;
  }
};

class BeliefPropagationInferenceEnginesBuilder : public InferenceEnginesBuilder {
protected:
  double dampening_lambda;
  double epsilon;
  long max_iter;

  std::vector<Scheduler<std::string>*> scheduler_ptrs;
  std::vector<InferenceGraph<std::string>*> graph_ptrs;
public:
  BeliefPropagationInferenceEnginesBuilder(double damp, double eps, long max_it):
    dampening_lambda(damp),
    epsilon(eps),
    max_iter(max_it)
  { }

  ~BeliefPropagationInferenceEnginesBuilder() {
    for (Scheduler<std::string>*sp : scheduler_ptrs)
      delete sp;
    scheduler_ptrs.clear();

    for (InferenceGraph<std::string>*ig : graph_ptrs)
      delete ig;
    graph_ptrs.clear();
  }

  std::vector<InferenceEngine<std::string>* > build_engines(const std::vector<std::vector<Dependency<std::string>*> > & dependencies_of_subgraphs) {
    for (Scheduler<std::string>*sp : scheduler_ptrs)
      delete sp;
    scheduler_ptrs.clear();

    for (InferenceGraph<std::string>*ig : graph_ptrs)
      delete ig;
    graph_ptrs.clear();
    
    std::vector<InferenceEngine<std::string>* > result(dependencies_of_subgraphs.size());

    for (unsigned long i=0; i<dependencies_of_subgraphs.size(); ++i) {
      const std::vector<Dependency<std::string>* > & deps = dependencies_of_subgraphs[i];
      BetheInferenceGraphBuilder<std::string> igb;
      for (Dependency<std::string>* dep : deps)
	igb.insert_dependency(*dep);
      Scheduler<std::string>* sched = new FIFOScheduler<std::string>(dampening_lambda, epsilon, max_iter);
      scheduler_ptrs.push_back(sched);

      InferenceGraph<std::string> *ig_ptr = new InferenceGraph<std::string>(igb.to_graph());
      graph_ptrs.push_back(ig_ptr);

      sched->add_ab_initio_edges(*ig_ptr);
      result[i] = new BeliefPropagationInferenceEngine<std::string>(*sched, *ig_ptr);
    }
    
    return result;
  }
};

#endif
