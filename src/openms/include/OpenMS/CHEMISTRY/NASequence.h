#ifndef OPENMS_CHEMISTRY_NASEQUENCE_H
#define OPENMS_CHEMISTRY_NASEQUENCE_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>

// to create iterators more easily
#include <boost/iterator/iterator_adaptor.hpp>

#include <vector>
#include <iosfwd>

namespace OpenMS
{



class OPENMS_DLLAPI NASequence
{

public:
  class ribonucleotide_iterator
    : public boost::iterator_adaptor<
      ribonucleotide_iterator
      , std::vector<Ribonucleotide *>::iterator // Base
      , boost::use_default                      // Value
      , boost::use_default                      // CategoryOrTraversal
      , boost::use_default                      // Reference
      , boost::use_default                      // Difference
    >
  {
  private:
    /*
    struct enabler {};  // a private type avoids misuse

    template <class OtherValue>
    ribonucleotide_iterator(
      ribonucleotide_iterator<OtherValue> const& other
      , typename boost::enable_if<
      boost::is_convertible<OtherValue*, std::vector<Ribonucleotide *>::iterator*>
      , enabler
      >::type = enabler()
    )
      : ribonucleotide_iterator::iterator_adaptor_(other.base()) {}
*/
  public:
    ribonucleotide_iterator()
      : ribonucleotide_iterator::iterator_adaptor_() {}

    explicit ribonucleotide_iterator(const std::vector<Ribonucleotide *>::iterator& p)
      : ribonucleotide_iterator::iterator_adaptor_(p) {}

/*
    explicit ribonucleotide_iterator(const ribonucleotide_iterator::iterator_adapter_::base_type& p)
      : ribonucleotide_iterator::iterator_adaptor_(p) {}
*/
  private:
    friend class boost::iterator_core_access;
  };

public:
  NASequence(); //default constructor
  NASequence(const NASequence& seq); // copy constructor
  NASequence& operator=(const NASequence& rhs); //assignment operator

  NASequence(std::vector<Ribonucleotide*> s, RibonucleotideChainEnd fivePrime, RibonucleotideChainEnd threePrime); //full constructor
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


  NASequence getPrefix(Size index) const;
  NASequence getSuffix(Size index) const;

  ribonucleotide_iterator begin() { return ribonucleotide_iterator(s_.begin()); }
  ribonucleotide_iterator end() { return ribonucleotide_iterator(s_.end()); }
  ribonucleotide_iterator cbegin() const { return ribonucleotide_iterator(s_.cbegin()); }
  ribonucleotide_iterator cend() const { return ribonucleotide_iterator(s_.cend()); }
private:
  std::vector<Ribonucleotide*> s_;
  RibonucleotideChainEnd fivePrime_;
  RibonucleotideChainEnd threePrime_;
};
OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const NASequence& nucleotide);

OPENMS_DLLAPI std::istream& operator>>(std::istream& os, const NASequence& nucleotide);
}
#endif // OPENMS_CHEMISTRY_NASEQUENCE_H
