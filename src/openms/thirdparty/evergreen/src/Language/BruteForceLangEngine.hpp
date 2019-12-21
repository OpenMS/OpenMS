#ifndef _BRUTEFORCELANGENGINE_HPP
#define _BRUTEFORCELANGENGINE_HPP

#include "LangEngine.hpp"

typedef struct BruteForceLangEngineT : public LangEngineT {

    void build(const std::vector<std::vector<Dependency<std::string>* > > & dependencies_of_subgraphs, double dampening, double epsilon, long max_iter) {

      for(InferenceGraph<std::string>* ptr: ig_ptrs)
        delete(ptr);
      ig_ptrs.clear();

      for(Scheduler<std::string>* ptr: sched_ptrs)
        delete(ptr);
      //sched_ptrs.clear();

      for(InferenceEngine<std::string>* ptr : eng_ptrs)
        delete(ptr);
      //eng_ptrs.clear();

      eng_ptrs.resize(dependencies_of_subgraphs.size());
      sched_ptrs.resize(dependencies_of_subgraphs.size());
      ig_ptrs.resize(dependencies_of_subgraphs.size());
      //#pragma omp parallel for // FIXME: Figure out why this makes valgrind throw a fit
      for (int i = 0; i < dependencies_of_subgraphs.size(); ++i) {
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
        eng_ptrs[i] = new BruteForceInferenceEngine<std::string>(table_deps, additive_deps, p);
      }

      is_built = true;
    }
} BruteForceLangEngine;

#endif
