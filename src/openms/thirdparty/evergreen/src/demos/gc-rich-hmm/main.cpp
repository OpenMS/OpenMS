#include <fstream>
#include <iostream>

#include "../../Evergreen/evergreen.hpp"
#include "HMM.hpp"
#include "HMMScheduler.hpp"

const double p = std::numeric_limits<double>::infinity();   // constant for p-norm approximation

std::string load_sequence(std::string file) {
  std::ifstream myfile(file);
  assert(myfile.is_open() == true && "Error: File not found");
  
  std::string line;
  std::string sequence;
  while ( std::getline(myfile,line) ) {
    std::stringstream ss_input(line);
    ss_input >> sequence;
  }
  return sequence;
}

int main() {
  
  // [Pr(H_1 = 0), Pr(H_1 = 1)]
  PMF prior({0L}, Tensor<double>({2ul}, {0.996, 0.004}));

  // [Pr(H_{i+1} = 0 | H_i = 0), Pr(H_{i+1} = 1 | H_i = 0), Pr(H_{i+1} = 0 | H_i = 1), Pr(H_{i+1} = 1 | H_i = 1)]
  PMF transition({0L,0L}, Tensor<double>({2ul,2ul},{0.99957, 0.00043,   0.00116954, 0.9988305}));

  // [Pr(D_i = G | H_i = 0), Pr(D_i = A | H_i = 0), Pr(D_i = T | H_i = 0), Pr(D_i = C | H_i = 0),
  //  Pr(D_i = G | H_i = 0), Pr(D_i = A | H_i = 0), Pr(D_i = T | H_i = 1), Pr(D_i = C | H_i = 1)]
  PMF emission({0L, 0L}, Tensor<double>({2ul,4ul},{0.209, 0.291,   0.291, 0.209,   0.331, 0.169,   0.169, 0.331}));
  
  // Data obtained from: https://www.ncbi.nlm.nih.gov/nuccore/CP000037
  std::string sequence = load_sequence("Shigella_boydii.fasta");

  std::cout << "RandomSubtreeScheduler" << std::endl;
  RandomSubtreeScheduler<std::string> rs_sched(0.0, 1e-3, -1ul);
  HMM hmm(prior, transition, emission, sequence, p, rs_sched);
  auto posteriors = hmm.solve();
  std::cout << std::endl;

  // This custom HMMScheduler is faster, but less general:
  std::cout << "HMMScheduler" << std::endl;
  HMMScheduler<std::string> hmm_sched;
  HMM hmm2(prior, transition, emission, sequence, p, hmm_sched);
  auto posteriors2 = hmm2.solve();
  std::cout << std::endl;
  
  return 0;
}
