#include "IsotopeQuantifier.hpp"

const Elements elements("element_isotope_list.txt");

void print_usage() {
  std::cerr << "Usage: isotope_quant <peak tsv filename> <intensity discretization> <intensity Gaussian std. dev> <maximum number copies for element> {missing, no_missing} <p> [maximum number of unique elements]" << std::endl;
  exit(1);
}

int main(int argc, char**argv) {
  if (argc != 7 && argc != 8) {
    print_usage();
  }

  std::string peak_file = argv[1];
  std::cerr << "peak_file = " << peak_file << std::endl;

  int intensity_discretization = atoi(argv[2]);
  std::cerr << "intensity_discretization = " << intensity_discretization << std::endl;
  
  double intensity_std_dev = atof(argv[3]);
  std::cerr << "intensity_std_dev = " << intensity_std_dev << std::endl;
  
  int maximum_copies_per_element = atoi(argv[4]);
  std::cerr << "maximum_copies_per_element = " << maximum_copies_per_element << std::endl;

  std::string missing_str = argv[5];
  bool include_missing = (missing_str == "missing");
  if( ! include_missing && (missing_str != "no_missing") )
    print_usage();

  double p = atof(argv[6]);

  int maximum_unique_elements = 0;
  if (argc == 8) {
    maximum_unique_elements = atoi(argv[7]);
    std::cerr << "maximum_unique_elements = " << maximum_unique_elements << std::endl;
  }
  
  FIFOScheduler<std::string> sched(0.01, 1e-16, 1000000ul);
  IsotopeQuantifier ms_solver(peak_file, elements, sched, p, intensity_discretization, intensity_std_dev, maximum_copies_per_element, include_missing, maximum_unique_elements);
  ms_solver.run_and_print_results();
  return 0;
}
