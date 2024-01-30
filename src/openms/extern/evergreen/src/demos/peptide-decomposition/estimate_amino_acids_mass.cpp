#include "MassPeptideSolver.hpp"

int main(int argc, char**argv) {
  if (argc != 5) {
    std::cout << "Usage: mass_pep_solver <observed mass> <mass discretization> <maximum peptide length> <p>" << std::endl;
    exit(1);
  }
 
  double mass = atof(argv[1]);

  double mass_discretization = atof(argv[2]);

  unsigned long max_length = atoi(argv[3]);
  double p = atof(argv[4]);

  FIFOScheduler<std::string> sched(0.01, 1e-8, 10000);
  MassPeptideSolver pep_solver(mass, p, max_length, mass_discretization, sched);
  pep_solver.solve_and_print();
  
  return 0;
}
