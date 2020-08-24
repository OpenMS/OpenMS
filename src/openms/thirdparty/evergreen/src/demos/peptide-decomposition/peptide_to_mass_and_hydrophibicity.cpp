#include "Peptide.hpp"

int main(int argc, char**argv) {
  if (argc != 2) {
    std::cout << "Usage: <peptide sequence>" << std::endl;
    exit(1);
  }

  std::string seq = argv[1];
  Peptide pep(seq);
  std::cout << pep.mass() << " " << pep.hydrophobicity() << std::endl;
  
  return 0;
}
