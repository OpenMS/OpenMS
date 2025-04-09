#ifndef _Peptide_HPP
#define _Peptide_HPP

#include <string>
#include <set>
#include <iostream>
#include <vector>
#include <map>
#include <assert.h>

class Peptide {
protected:
  std::string _amino_acids;

  double _mass;
  double _hydrophobicity;

  void verify_valid_characters() {
    std::set<char> amino_set(_amino_acids.begin(), _amino_acids.end());
    for (char c : _amino_acids) {
      if (amino_set.find(c) == amino_set.end()) {
	std::cerr << "Invalid character: " << c << std::endl;
	assert(false);
      }
    }
  }

  // Calculate the mass of the peptide.
  void init_mass() {
    std::map<char, double> amino_acid_to_mass;
    for (unsigned long i=0; i<amino_acids.size(); ++i)
      amino_acid_to_mass[amino_acids[i]] = masses[i];
    
    _mass = 0.0;
    for (char aa : _amino_acids) {
      assert(amino_acid_to_mass.find(aa) != amino_acid_to_mass.end() && "Error: Amino acid not found.");
      _mass += amino_acid_to_mass[aa];
    }
  }
  
  // Calculate the hydrophobicity of the peptide. 
  void init_hydrophobicity() {
    std::map<char, double> amino_acid_to_hydrophobicity;
    for (unsigned long i=0; i<amino_acids.size(); ++i)
      amino_acid_to_hydrophobicity[amino_acids[i]] = hydrophobicities[i];
    
    _hydrophobicity = 0.0;
    for (char aa : _amino_acids) {
      assert(amino_acid_to_hydrophobicity.find(aa) != amino_acid_to_hydrophobicity.end() && "Error: Amino acid not found.");
      _hydrophobicity += amino_acid_to_hydrophobicity[aa];
    }
  }
  
public:
  static const std::vector<char> amino_acids;
  static const std::vector<double> masses;
  static const std::vector<double> hydrophobicities;

  Peptide(const std::string & seq):
  _amino_acids(seq)
  {
    verify_valid_characters();
    init_mass();
    init_hydrophobicity();
  }

  unsigned long size() const {
    return _amino_acids.size();
  }

  char operator [] (unsigned long i) const {
    return _amino_acids[i];
  }
  
  const double & mass() const{
    return _mass;
  }
  
  const double & hydrophobicity() const{
    return _hydrophobicity;
  }

};

// {A:Ala, R:Arg, N:Asn, D:Asp, C:Cys, E:Glu, Q:Gln, G:Gly, H:His, I:Ile, L:Leu, K:Lys,
//  M:Met, F:Phe, P:Pro, S:Ser, T:Thr ,W:Trp , Y:Tyr, V:Val}
const std::vector<char> Peptide::amino_acids = {'A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V'};


// http://www.matrixscience.com/help/aa_help.html (average mass)
const std::vector<double> Peptide::masses = {71.0779, 156.1857, 114.1026,  115.0874, 103.1429, 129.114, 128.1292, 57.0513, 137.1393, 113.1576, 113.1576, 128.1723, 131.1961, 147.1739, 97.1152, 87.0773, 101.1039, 186.2099, 163.1733, 99.1311};

// wwHydrophobicity from
// https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html
const std::vector<double> Peptide::hydrophobicities = {-0.17, -0.81, -0.42, -1.23, 0.24, -2.02, -0.58, -0.01, -0.96, 0.31, 0.56, -0.99, 0.23, 1.13, -0.45, -0.13, -0.14, 1.85, 0.94, -0.07};

std::ostream & operator<<(std::ostream & os, const Peptide & rhs) {
  for (unsigned long i=0; i<rhs.size(); ++i)
    os << rhs[i];
  os << ": mass=" << rhs.mass() << " hydrophobicity=" << rhs.hydrophobicity();
  return os;
}

#endif
