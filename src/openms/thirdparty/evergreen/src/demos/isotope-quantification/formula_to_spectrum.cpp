#include "IsotopeQuantifier.hpp"

const Elements elements("element_isotope_list.txt");

void print_usage() {
  std::cerr << "Usage:\n";
  std::cerr << "\tformula2spectrum discretize_mass=15 Ca=10 [Ar=2 ...]" << std::endl;
  exit(1);
}

int main(int argc, char**argv) {
  if (argc <= 1)
    print_usage();

  double discretization = -1;
  std::string exact_or_disc = argv[1];
  int eq = exact_or_disc.find("=");
  if (eq == -1)
    print_usage();
  else {
    exact_or_disc[eq] = ' ';
    std::istringstream ist(exact_or_disc);
    std::string garbage;
    ist >> garbage;
    if (garbage != "discretize_mass")
      print_usage();
    
    ist >> discretization;
    
    if (discretization <= 0) {
      std::cerr << "discretize_mass must be >0" << std::endl;
      return 1;
    }
  }

  std::map<std::string, unsigned int> element_to_count;

  for (int i=2; i<argc; ++i) {
    std::string element_and_count = argv[i];
    int eq = element_and_count.find("=");

    if (eq == -1)
      print_usage();

    element_and_count[eq] = ' ';

    std::string element;
    int count;
    
    std::istringstream ist(element_and_count);
    ist >> element >> count;

    if (count <= 0) {
      std::cerr << "Abundance of element must be integer > 0" << std::endl;
      return 1;
    }

    if (element_to_count.find(element) != element_to_count.end()) {
      std::cerr << "Error: " + element + "added multiple times" << std::endl;
      return 1;
    }

    element_to_count[element] = count;
  }

  std::map<double, double> peaks;

  // Discretize:
  // Use false to ignore unobserved peaks
  peaks = IsotopeQuantifier::mass_discretized_theoretical_peaks_from_chemical_formula(element_to_count, elements, discretization, false);

  // Print:
  std::cout << "mass_discretization " << discretization << std::endl;
  for (const std::pair<double, double> & x : peaks) {
    std::cout << x.first << "\t" << x.second << std::endl;
  }
  
  return 0;
}
