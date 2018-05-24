#include "HydrophobicityPeptideSolver.hpp"

int main(int argc, char**argv) {
  if (argc != 5) {
    std::cout << "Usage: hydro_pep_solver <observed hydrophobicity> <hydrophobicity discretization> <maximum peptide length> <p>" << std::endl;
    exit(1);
  }
 
  double hydrophobicity = atof(argv[1]);

  double hydrophobicity_discretization = atof(argv[2]);

  unsigned long max_length = atoi(argv[3]);
  double p = atof(argv[4]);

  FIFOScheduler<std::string> sched(0.01, 1e-8, 10000);
  HydrophobicityPeptideSolver pep_solver(hydrophobicity, p, max_length, hydrophobicity_discretization, sched);
  pep_solver.solve_and_print();
  
  return 0;
}
