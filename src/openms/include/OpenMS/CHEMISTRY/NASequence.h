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
  class Iterator;

  /** @brief ConstIterator
*/
  class OPENMS_DLLAPI ConstIterator
  {
  public:
    // TODO Iterator constructor for ConstIterator
    typedef Ribonucleotide value_type;
    typedef const value_type& const_reference;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef std::vector<const value_type*>::difference_type difference_type;
    typedef const value_type* pointer;
    typedef std::random_access_iterator_tag iterator_category;

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    ConstIterator() = default;

    /// detailed constructor with pointer to the vector and offset position
    ConstIterator(const std::vector<const Ribonucleotide*>* vec_ptr, difference_type position)
    {
      vector_ = vec_ptr;
      position_ = position;
    }

    /// copy constructor
    ConstIterator(const ConstIterator& rhs) :
      vector_(rhs.vector_),
      position_(rhs.position_)
    {
    }

    /// copy constructor from Iterator
    ConstIterator(const NASequence::Iterator& rhs) :
      vector_(rhs.vector_),
      position_(rhs.position_)
    {
    }

    /// destructor
    virtual ~ConstIterator() {}

    //@}

    /// assignment operator
    ConstIterator& operator=(const ConstIterator& rhs)
    {
      if (this != &rhs)
      {
        position_ = rhs.position_;
        vector_ = rhs.vector_;
      }
      return *this;
    }

    /** @name Operators
    */
    //@{
    /// dereference operator
    const_reference operator*() const { return *(*vector_)[position_]; }

    /// dereference operator
    const_pointer operator->() const { return (*vector_)[position_]; }

    /// forward jump operator
    const ConstIterator operator+(difference_type diff) const { return ConstIterator(vector_, position_ + diff); }

    difference_type operator-(ConstIterator rhs) const { return position_ - rhs.position_; }

    /// backward jump operator
    const ConstIterator operator-(difference_type diff) const
    {
      return ConstIterator(vector_, position_ - diff);
    }

    /// equality comparator
    bool operator==(const ConstIterator& rhs) const
    {
      return vector_ == rhs.vector_ && position_ == rhs.position_;
    }

    /// inequality operator
    bool operator!=(const ConstIterator& rhs) const
    {
      return vector_ != rhs.vector_ || position_ != rhs.position_;
    }

    /// increment operator
    ConstIterator& operator++()
    {
      ++position_;
      return *this;
    }

    /// decrement operator
    ConstIterator& operator--()
    {
      --position_;
      return *this;
    }

    //@}

  protected:

    // pointer to the vector
    const std::vector<const Ribonucleotide*>* vector_;

    // position in the vector
    difference_type position_;
  };


  /** @brief Iterator class
  */
  class OPENMS_DLLAPI Iterator
  {
  public:

    friend class NASequence::ConstIterator;

    typedef Ribonucleotide value_type;
    typedef const value_type& const_reference;
    typedef value_type& reference;
    typedef const value_type* const_pointer;
    typedef const value_type* pointer;
    typedef std::vector<const value_type*>::difference_type difference_type;

    /** @name Constructors and destructors
    */
    //@{
    Iterator() = default;

    /// detailed constructor with pointer to the vector and offset position
    Iterator(std::vector<const Ribonucleotide*>* vec_ptr, difference_type position)
    {
      vector_ = vec_ptr;
      position_ = position;
    }

    /// copy constructor
    Iterator(const Iterator& rhs) :
      vector_(rhs.vector_),
      position_(rhs.position_)
    {
    }

    /// destructor
    virtual ~Iterator() {}

    //@}

    /// assignment operator
    Iterator& operator=(const Iterator& rhs)
    {
      if (this != &rhs)
      {
        position_ = rhs.position_;
        vector_ = rhs.vector_;
      }
      return *this;
    }

    /** @name Operators
    */
    //@{
    /// dereference operator
    const_reference operator*() const { return *(*vector_)[position_]; }

    /// dereference operator
    const_pointer operator->() const { return (*vector_)[position_]; }

    /// mutable dereference operator
    pointer operator->() { return (*vector_)[position_]; }

    /// forward jump operator
    const Iterator operator+(difference_type diff) const { return Iterator(vector_, position_ + diff); }

    difference_type operator-(Iterator rhs) const { return position_ - rhs.position_; }

    /// backward jump operator
    const Iterator operator-(difference_type diff) const { return Iterator(vector_, position_ - diff); }

    /// equality comparator
    bool operator==(const Iterator& rhs) const { return std::tie(vector_,position_) == std::tie(rhs.vector_, rhs.position_); }

    /// inequality operator
    bool operator!=(const Iterator& rhs) const { return !this->operator==(rhs); }

    /// increment operator
    Iterator& operator++()
    {
      ++position_;
      return *this;
    }

    /// decrement operator
    Iterator& operator--()
    {
      --position_;
      return *this;
    }

    //@}

  protected:

    std::vector<const Ribonucleotide*>* vector_;

    // position in the vector
    difference_type position_;
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

  void set(size_t index, const Ribonucleotide* r);

  bool hasFivePrimeModification() const;

  void setFivePrimeModification(const RibonucleotideChainEnd *r);

  const Ribonucleotide *getFivePrimeModification();

  bool hasThreePrimeModification() const;

  void setThreePrimeModification(const RibonucleotideChainEnd *r);

  const Ribonucleotide *getThreePrimeModification();

  NASequence getPrefix(Size index) const;

  NASequence getSuffix(Size index) const;

  inline Iterator begin() { return Iterator(&s_, 0); }
  inline ConstIterator begin() const { return ConstIterator(&s_, 0); }
  inline Iterator end() { return Iterator(&s_, (Int) s_.size()); }
  inline ConstIterator end()  const { return ConstIterator(&s_, (Int) s_.size()); }
  inline ConstIterator cbegin() const { return ConstIterator(&s_, 0); }
  inline ConstIterator cend() const { return ConstIterator(&s_, (Int) s_.size()); }

  /// return reference to the residue at given position
  const Ribonucleotide& operator[](size_t index) const { return s_[index]; }

private:
  std::vector<const Ribonucleotide*> s_;
  RibonucleotideChainEnd fivePrime_;
  RibonucleotideChainEnd threePrime_;
};
OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const NASequence& nucleotide);

OPENMS_DLLAPI std::istream& operator>>(std::istream& os, const NASequence& nucleotide);
}
#endif // OPENMS_CHEMISTRY_NASEQUENCE_H
