#ifndef OPENMS_CHEMISTRY_NASEQUENCE_H
#define OPENMS_CHEMISTRY_NASEQUENCE_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>


#include <vector>
#include <iosfwd>

namespace OpenMS
{

class OPENMS_DLLAPI NASequence
{

public:
  NASequence(); //default constructor
  NASequence(const NASequence& seq); // copy constructor
  NASequence& operator=(const NASequence& rhs); //assignment operator
  bool operator==(const NASequence& rhs) const;
  virtual ~NASequence(); //destructor
  void setSequence(std::vector<Ribonucleotide*>& s);
  std::vector<Ribonucleotide*> getSequence() const;
  size_t size() const;
  double getMonoWeight(Ribonucleotide::RiboNucleotideFragmentType type = Ribonucleotide::Full, Int charge = 0) const;
  bool empty() const;
  EmpiricalFormula getFormula(OpenMS::Ribonucleotide::RiboNucleotideFragmentType type = Ribonucleotide::Full, Int charge = 0) const;

  //TODO:implement
  void set(size_t index, const Ribonucleotide* r);

  bool hasFivePrimeModification() const;
  void setFivePrimeModification(const Ribonucleotide* r);
  Ribonucleotide* getFivePrimeModification();

  bool hasThreePrimeModification() const;
  void setThreePrimeModification(const Ribonucleotide* r);
  Ribonucleotide* getThreePrimeModification();


  typedef std::vector<Ribonucleotide*>::iterator iterator;
  iterator begin() { return s_.begin(); }
  iterator end() { return s_.end(); }

  typedef std::vector<Ribonucleotide*>::const_iterator const_iterator;
  const_iterator cbegin() const { return s_.cbegin(); }
  const_iterator cend() const { return s_.cend(); }

  NASequence getPrefix(Size index) const;
  NASequence getSuffix(Size index) const;

private:
  std::vector<Ribonucleotide*> s_;
  RibonucleotideChainEnd fivePrime_;
  RibonucleotideChainEnd threePrime_;
};
OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const NASequence& nucleotide);

OPENMS_DLLAPI std::istream& operator>>(std::istream& os, const NASequence& nucleotide);
}
#endif // OPENMS_CHEMISTRY_NASEQUENCE_H
