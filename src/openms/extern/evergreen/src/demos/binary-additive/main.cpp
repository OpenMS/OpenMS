#include "../../Evergreen/evergreen.hpp"
#include "../../Utility/inference_utilities.hpp"

int main(int argc, char**argv) {
  if (argc != 2) {
    std::cerr << "Usage: binary_tree <LOG_N>" << std::endl;
    exit(1);
  }

  int log_n = atoi(argv[1]);

  const double p=std::numeric_limits<double>::infinity();
  BetheInferenceGraphBuilder<unsigned long> igb;

  const unsigned long n=1ul<<log_n;

  std::cout << "Creating dependencies..." << std::endl;
  for (unsigned long i=0; i<=n; ++i) {
    double prob0 = rand() % 1000 / 999.0;
    double prob[] = {prob0+0.01, 1-prob0+0.01};
    //    LabeledPMF<unsigned long> lpmf({i},PMF({0L},Tensor<double>::from_array(prob)));
    LabeledPMF<unsigned long> lpmf({i}, PMF({0L},Tensor<double>::from_array(prob)));
    igb.insert_dependency( TableDependency<unsigned long>(lpmf,p) );
  }

  std::vector<std::vector<unsigned long> > inputs;
  for (unsigned long i=0; i<n; ++i)
    inputs.push_back({i});

  igb.insert_dependency( AdditiveDependency<unsigned long>(inputs,{n},p) );

  std::cout << "Constructing graph..." << std::endl;
  InferenceGraph<unsigned long> ig = igb.to_graph();

  FIFOScheduler<unsigned long> fifo(0.001, 1e-16, 1ul<<32);
  fifo.add_ab_initio_edges(ig);
  BeliefPropagationInferenceEngine<unsigned long> bpie(fifo, ig);
  estimate_and_print_posteriors(bpie, {{0}, {n}});
}
