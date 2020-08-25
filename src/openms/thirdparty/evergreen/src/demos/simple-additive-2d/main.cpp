#include <fstream>

#include "../../Evergreen/evergreen.hpp"
#include "../../Utility/inference_utilities.hpp"
#include "../../Utility/graph_to_dot.hpp"

const double p=16;

void solve_1d_bethe(const std::vector<LabeledPMF<std::string> > & inputs, const LabeledPMF<std::string> & output, const std::vector<std::vector<std::string> > & vars_for_posteriors) {
  BetheInferenceGraphBuilder<std::string> igb;

  for (const LabeledPMF<std::string> & lpmf : inputs)
    igb.insert_dependency( TableDependency<std::string>(lpmf,p) );
  igb.insert_dependency( TableDependency<std::string>(output,p) );

  // 2x AdditiveDependency types:
  std::vector<std::vector<std::string> > input_vars_0;
  for (const LabeledPMF<std::string> & lpmf : inputs)
    input_vars_0.push_back( {lpmf.ordered_variables()[0]} );
  igb.insert_dependency( AdditiveDependency<std::string>(input_vars_0, {output.ordered_variables()[0]}, p) );

  std::vector<std::vector<std::string> > input_vars_1;
  for (const LabeledPMF<std::string> & lpmf : inputs)
    input_vars_1.push_back( {lpmf.ordered_variables()[1]} );
  igb.insert_dependency( AdditiveDependency<std::string>(input_vars_1, {output.ordered_variables()[1]}, p) );
  
  InferenceGraph<std::string> ig = igb.to_graph();

  FIFOScheduler<std::string> fifo(0.01, 1e-8, 10000);
  fifo.add_ab_initio_edges(ig);
  BeliefPropagationInferenceEngine<std::string> bpie(fifo, ig);
  estimate_and_print_posteriors(bpie, vars_for_posteriors);

  write_graph_to_dot_file(ig, "bethe_1d.dot");
}

void solve_2d_bethe(const std::vector<LabeledPMF<std::string> > & inputs, const LabeledPMF<std::string> & output, const std::vector<std::vector<std::string> > & vars_for_posteriors) {
  BetheInferenceGraphBuilder<std::string> igb;

  for (const LabeledPMF<std::string> & lpmf : inputs)
    igb.insert_dependency( TableDependency<std::string>(lpmf,p) );
  igb.insert_dependency( TableDependency<std::string>(output,p) );

  // AdditiveDependency:
  std::vector<std::vector<std::string> > input_vars;
  for (const LabeledPMF<std::string> & lpmf : inputs)
    input_vars.push_back(lpmf.ordered_variables());
  igb.insert_dependency( AdditiveDependency<std::string>(input_vars, output.ordered_variables(), p) );

  InferenceGraph<std::string> ig = igb.to_graph();

  FIFOScheduler<std::string> fifo(0.01, 1e-8, 10000);
  fifo.add_ab_initio_edges(ig);
  BeliefPropagationInferenceEngine<std::string> bpie(fifo, ig);
  estimate_and_print_posteriors(bpie, vars_for_posteriors);

  write_graph_to_dot_file(ig, "bethe_2d.dot");
}

void solve_2d_exact(const std::vector<LabeledPMF<std::string> > & inputs, const LabeledPMF<std::string> & output, const std::vector<std::vector<std::string> > & vars_for_posteriors) {
  std::vector<ContextFreeMessagePasser<std::string>* > input_mps;
  std::vector<std::vector<std::string>* > input_labels;
  for (const LabeledPMF<std::string> & lpmf : inputs) {
    input_mps.push_back( new HUGINMessagePasser<std::string>(lpmf, p) );
    input_labels.push_back( new std::vector<std::string>(lpmf.ordered_variables()) );
  }

  ContextFreeMessagePasser<std::string>*output_mp = new HUGINMessagePasser<std::string>(output, p);
  std::vector<std::string>*output_label = new std::vector<std::string>(output.ordered_variables());

  ConvolutionTreeMessagePasser<std::string>*ctmp = new ConvolutionTreeMessagePasser<std::string>(input_mps, input_labels, output_mp, output_label, 2, p);

  std::vector<MessagePasser<std::string>* > mps;
  for (ContextFreeMessagePasser<std::string>*hmp : input_mps)
    mps.push_back(hmp);
  mps.push_back(output_mp);
  mps.push_back(ctmp);

  InferenceGraph<std::string> ig(std::move(mps));
  
  FIFOScheduler<std::string> fifo(0.01, 1e-8, 10000);
  fifo.add_ab_initio_edges(ig);
  BeliefPropagationInferenceEngine<std::string> bpie(fifo, ig);
  estimate_and_print_posteriors(bpie, vars_for_posteriors);

  write_graph_to_dot_file(ig, "exact_2d.dot");
}

int main() {
  LabeledPMF<std::string> av({"A","V"},PMF({2L,1L},Tensor<double>({3,3},{1,10,9,3,7,2,1,2,6})));
  LabeledPMF<std::string> bw({"B","W"},PMF({1L,0L},Tensor<double>({3,3},{1,2,3,4,5,6,7,8,9})));
  LabeledPMF<std::string> cx({"C","X"},PMF({-1L,0L},Tensor<double>({3,2},{2,8,4,1,2,3})));
  LabeledPMF<std::string> dy({"D","Y"},PMF({0L,0L},Tensor<double>({2,3},{7,5,2,5,6,3})));
  LabeledPMF<std::string> ez({"E","Z"},PMF({0L,1L},Tensor<double>({2,3},{10,3,6,4,1,7})));
  std::cout << av << std::endl;
  std::cout << bw << std::endl;
  std::cout << cx << std::endl;
  std::cout << dy << std::endl;
  std::cout << ez << std::endl;
  // (A,V) = (B,W) + (C,X) + (D,Y) + (E,Z)

  std::cout << "2x 1D Convolution trees (Bethe construction)" << std::endl;
  solve_1d_bethe({bw, cx, dy, ez}, av, {{"A","V"}, {"E","Z"}});
  std::cout << std::endl;

  std::cout << "2D Convolution tree (Bethe construction with 1D bottlenecks)" << std::endl;
  solve_2d_bethe({bw, cx, dy, ez}, av, {{"A","V"}, {"E","Z"}});
  std::cout << std::endl;

  std::cout << "2D Convolution tree (exact)" << std::endl;
  solve_2d_exact({bw, cx, dy, ez}, av, {{"A","V"}, {"E","Z"}});
  std::cout << std::endl;
}
