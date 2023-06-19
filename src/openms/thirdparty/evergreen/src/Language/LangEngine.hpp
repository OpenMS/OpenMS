#ifndef _LANGENGINE_HPP
#define _LANGENGINE_HPP

#include "InferenceEnginesBuilder.hpp"
#include "LangDigraph.hpp"
#include <unordered_map>
#include <fstream>
#include "../Utility/graph_to_dot.hpp"

typedef struct LangEngineT {
  LangGraph<std::string> graph;
  
  InferenceEnginesBuilder*ieb_ptr;

  std::unordered_map<std::string, std::vector<unsigned long> > var_to_graphs_containing;
  std::vector<Dependency<std::string>* > dependencies;

  std::vector<InferenceEngine<std::string>* > engine_ptrs;

  LangEngineT(const double & default_damp, const double & default_eps, const long & default_max_iter):
    ieb_ptr(new BeliefPropagationInferenceEnginesBuilder(default_damp, default_eps, default_max_iter))
  { }

  ~LangEngineT() {
    delete ieb_ptr;
    
    for (Dependency<std::string>*dep_ptr : dependencies)
      delete dep_ptr;
  }

  void insert_dependency(Dependency<std::string> * dep) {
    dependencies.push_back(dep);
    const std::vector<std::string> & vars_used = dep->get_all_variables_used();
    for (const std::string & var : vars_used)
      var_to_graphs_containing[var].push_back(dependencies.size()-1);
  }

  void set_engine(InferenceEnginesBuilder*ieb) {
    delete ieb_ptr;
    ieb_ptr = ieb;
  }

  void print(const std::vector<std::vector<std::string> > & result_vars) {
    const std::vector<std::string> & flat_result_vars = flatten(result_vars);
    const std::vector<std::vector<std::string> > & partitioned_subgraphs = partition_into_subgraphs(flat_result_vars);

    const std::vector<std::vector<Dependency<std::string>* > > & dependencies_of_subgraphs = get_dependencies_of_subgraphs(partitioned_subgraphs);
    engine_ptrs = ieb_ptr->build_engines(dependencies_of_subgraphs);

    std::unordered_map<std::string, int> var_to_graph_number;
    for (unsigned long i=0; i<partitioned_subgraphs.size(); ++i) {
      const std::vector<std::string> & vars_in_connected_graph = partitioned_subgraphs[i];
      for (const std::string & var : vars_in_connected_graph)
      	var_to_graph_number[var] = i;
    }

    // outer vector is for subgraph
    // middle vector is for the possibly many tuples on which we want posteriors
    // inner vector is a tuple of variables.
    std::vector<std::vector<std::vector<std::string> > > printed_partitioned_subgraphs;
    printed_partitioned_subgraphs.resize(dependencies_of_subgraphs.size());
    for (const std::vector<std::string> & result_var : result_vars) {
      int graph_num = var_to_graph_number[result_var[0]];
      for (const std::string & var : result_var)
        if (var_to_graph_number[var] != graph_num)
          std::cerr << "ERROR: Printing error, tried to print posteriors on set of vars that belong in different subgraphs." << std::endl;
      printed_partitioned_subgraphs[var_to_graph_number[result_var[0]]].push_back(result_var);
    }

    std::vector<std::vector<LabeledPMF<std::string> > > all_results_to_print(printed_partitioned_subgraphs.size());
    #pragma omp parallel for
    for (unsigned long i=0; i<printed_partitioned_subgraphs.size(); ++i)
      all_results_to_print[i] = engine_ptrs[i]->estimate_posteriors(printed_partitioned_subgraphs[i]);
    for (const std::vector<LabeledPMF<std::string> > & results_to_print : all_results_to_print)
      for (const LabeledPMF<std::string> & result_to_print : results_to_print)
        std::cout << result_to_print << std::endl; 
  }
  
  // ----------------------------------------------------
  // not for client use (i.e., private helper functions):
  // ----------------------------------------------------

  std::vector<std::vector<std::string> > partition_into_subgraphs(const std::vector<std::string> & result_vars) {
    std::vector<std::vector<std::string> > partitioned_subgraphs;
    std::set<std::string> vars_visited;
    std::set<std::string> result_vars_visited; 
    for(const std::string & result_var : result_vars) {
      if (result_vars_visited.find(result_var) == result_vars_visited.end()) {
        std::vector<std::string> subgraph;
        partition_into_single_subgraph(result_var, subgraph, vars_visited, result_vars, result_vars_visited);
        if (subgraph.size() > 0) 
          partitioned_subgraphs.push_back(subgraph);
      }
    }
    return partitioned_subgraphs;
  }

  void partition_into_single_subgraph(const std::string & result_var, std::vector<std::string> & subgraph, std::set<std::string> & vars_visited, const std::vector<std::string> & result_vars, std::set<std::string> & result_vars_visited) {
    if (vars_visited.find(result_var) == vars_visited.end()) {
      vars_visited.insert(result_var);
      subgraph.push_back(result_var);
      for(const int & dep_index : var_to_graphs_containing[result_var]) {
        std::vector<std::string> adj_vars = dependencies[dep_index]->get_all_variables_used();      
        for (const std::string & adj_var : adj_vars) {
          if (result_vars_visited.find(adj_var) == result_vars_visited.end()) {
            if (find(result_vars.begin(), result_vars.end(), adj_var) != result_vars.end())
              result_vars_visited.insert(adj_var);
            partition_into_single_subgraph(adj_var, subgraph, vars_visited, result_vars, result_vars_visited);
          }
        }
      }
    }
  }

  void get_dependencies_in_single_subgraph(const std::string & var, std::vector<bool> & deps_visited, std::vector<Dependency<std::string>* > & connected_dependencies) {
    for (const int & dep_index : var_to_graphs_containing[var]) {
      if (!deps_visited[dep_index]) {
        deps_visited[dep_index] = true;
        connected_dependencies.push_back(dependencies[dep_index]);
        const std::vector<std::string> & vars_used = dependencies[dep_index]->get_all_variables_used();
        for (const std::string & var_used : vars_used) {
          get_dependencies_in_single_subgraph(var_used, deps_visited, connected_dependencies);
        }
      }
    }
    if (connected_dependencies.size() == 0) {
      std::cerr << "ERROR: printing error, tried to print posteriors on var " << var << " that doesn't exist in any graph" << std::endl;
    }
  }

  std::vector<std::vector<Dependency<std::string>* > > get_dependencies_of_subgraphs(const std::vector<std::vector<std::string> > & partitioned_subgraphs) {
    std::vector<std::vector<Dependency<std::string>* > > dependencies_of_subgraphs;
    std::vector<bool> deps_visited;
    deps_visited.resize(dependencies.size());
    for(const std::vector<std::string> & subgraph : partitioned_subgraphs) {
      std::vector<Dependency<std::string>* > dependencies_in_single_graph;
      get_dependencies_in_single_subgraph(subgraph[0], deps_visited, dependencies_in_single_graph);
      dependencies_of_subgraphs.push_back(dependencies_in_single_graph);
    } 
    return dependencies_of_subgraphs;
  }

  void recompute_and_print_normalization_constant() {
    std::vector<std::string> all_vars;
    for (const std::pair<std::string, std::vector<unsigned long> > & p : var_to_graphs_containing)
      all_vars.push_back(p.first);
    
    std::vector<std::vector<std::string> > partitioned_subgraphs = partition_into_subgraphs(all_vars);
    std::vector<std::vector<Dependency<std::string>* > > all_partitioned_dependencies = get_dependencies_of_subgraphs(partitioned_subgraphs);

    engine_ptrs = ieb_ptr->build_engines(all_partitioned_dependencies);

    std::vector<std::vector<std::vector<std::string> > > singleton_partitions;
    for (std::vector<std::string> & subgraph_vars : partitioned_subgraphs)
      singleton_partitions.push_back(make_singletons(subgraph_vars));

    #pragma omp parallel for
    for (unsigned long i=0; i<partitioned_subgraphs.size(); ++i)
      engine_ptrs[i]->estimate_posteriors(singleton_partitions[i]);

    print_normalization_constant();
  }

  void print_normalization_constant() {
    double log_nc = 0.0;
    for (InferenceEngine<std::string>*ie : engine_ptrs)
      log_nc += ie->log_normalization_constant();

    std::cout << "Log probability of model: " << log_nc << std::endl;

    for (unsigned long i=0; i<engine_ptrs.size(); ++i)
      delete engine_ptrs[i];
  }

  void save_graph(char*str) {
    BetheInferenceGraphBuilder<std::string> igb;
    for (Dependency<std::string>* dep : dependencies)
      igb.insert_dependency(*dep);

    std::ofstream fout(str);
    graph_to_dot(igb.to_graph(), fout);
    fout.close();
  }
} LangEngine;

#endif
