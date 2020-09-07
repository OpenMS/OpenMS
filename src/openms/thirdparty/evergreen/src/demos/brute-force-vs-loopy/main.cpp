#include <string>
#include "../../Evergreen/evergreen.hpp"
#include "../../Utility/inference_utilities.hpp"

const double p = 16.0;

void brute_force(const std::vector<TableDependency<std::string> > & deps, const std::vector<std::vector<std::string> > & vars) {
  BruteForceInferenceEngine<std::string> bf(deps,p);
  estimate_and_print_posteriors(bf, vars);

  std::cout << std::endl;
}

void loopy(const std::vector<TableDependency<std::string> > & deps, const std::vector<std::vector<std::string> > & vars) {
  BetheInferenceGraphBuilder<std::string> igb;
  for (const TableDependency<std::string> & td : deps)
    igb.insert_dependency(td);
  InferenceGraph<std::string> ig = igb.to_graph();

  FIFOScheduler<std::string> sched(0.0, 1e-8, 10000);
  sched.add_ab_initio_edges(ig);
  BeliefPropagationInferenceEngine<std::string> bpie(sched, ig);
  estimate_and_print_posteriors(bpie, vars);
 
  std::cout << std::endl;
}

int main() {
  TableDependency<std::string> td1(LabeledPMF<std::string>({"a", "b"}, PMF({0L,0L}, Tensor<double>({2ul,2ul}, {.87, .13, .74, .26}))), p);
  TableDependency<std::string> td2(LabeledPMF<std::string>({"b", "c"}, PMF({0L,0L}, Tensor<double>({2ul,2ul}, {.4, .2, .1, .3}))), p);
  TableDependency<std::string> td3(LabeledPMF<std::string>({"a", "c"}, PMF({0L,0L}, Tensor<double>({2ul,2ul}, {.3, .1, .45, .15}))), p);

  std::cout << "Brute force" << std::endl;
  brute_force({td1,td2,td3}, {{"a","b"}, {"b","c"}});

  std::cout << "Loopy belief propagation" << std::endl;
  loopy({td1,td2,td3}, {{"a","b"}, {"b","c"}});
  
  return 0;
}
