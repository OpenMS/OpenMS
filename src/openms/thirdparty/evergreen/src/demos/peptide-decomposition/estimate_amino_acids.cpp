#include "PeptideSolver.hpp"

int main(int argc, char**argv) {
  if (argc != 7) {
    std::cout << "Usage: pep_solver <observed mass> <observed hydrophobicity> <mass discretization> <hydrophobicity discretization> <maximum peptide length> <p>" << std::endl;
    exit(1);
  }
 
  double mass = atof(argv[1]);

  double hydrophobicity = atof(argv[2]);

  double mass_discretization = atof(argv[3]);
  double hydrophobicity_discretization = atof(argv[4]);

  unsigned long max_length = atoi(argv[5]);
  double p = atof(argv[6]);

  FIFOScheduler<std::string> sched(0.01, 1e-8, 10000);
  PeptideSolver pep_solver(mass, hydrophobicity, p, max_length, mass_discretization, hydrophobicity_discretization, sched);
  pep_solver.solve_and_print();
  
  return 0;
}
