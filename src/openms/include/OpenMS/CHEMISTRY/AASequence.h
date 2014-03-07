// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_AASEQUENCE_H
#define OPENMS_CHEMISTRY_AASEQUENCE_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <vector>
#include <iostream>

namespace OpenMS
{

  //forward declarations
  class ResidueDB;
  class ResidueModification;

  /**
      @brief Representation of a peptide/protein sequence

      This class represents amino acid sequences in %OpenMS. Basically a AASequence instance
      consists of a sequence of residues. The residues are represented as instances of
      Residue. Each amino acid has only one instance which is accessible using the ResidueDB instance (singleton).

      A critical property of amino acid sequence is that they can be modified. Which means that one or more
      amino acids are chemically modified, e.g. oxidized. This is represented via Residue instances which carry
      a ResidueModification object. This is also handled in the ResidueDB.

      If one wants to specify a AASequence the easiest way is simply writing the amino acid sequence. For example
      AASequence seq("DFPIANGER") is sufficient to create a instance of AASequence with DFPIANGER as peptide.

      Modifications are specified using a unique string identifier present in the ModificationsDB in brackets
      after the modified amino acid or by providing the mass of the residue in square brackets. For example 
      AASequence("DFPIAM(Oxidation)GER") creates an instance of the peptide DFPIAMGER with an oxidized methionine 
      (AASequence("DFPIAM[+16]GER") and AASequence("DFPIAM[147]GER") are equivalent). N-terminal modifications 
      are specified by writing the modification as prefix to the sequence. C-terminal modifications are specified by 
      writing the modification as suffix. C-terminal modifications are distinguished from modifications of the last amino
      acid by considering the specificity of the modification as stored in ModificationsDB.

      Note there is a subtle difference between AASequence("DFPIAM[+16]GER") and AASequence("DFPIAM[+15.9949]GER") - while 
      the former will try to find the _first_ modification matching to a mass difference of 16 +/- 0.5, the latter will 
      try to find the closest matching modification to the exact mass. This usually gives the intended results.

      Arbitrary/unknown AA's (usually due to an unknown modification) can be specified using tags preceded by X: 'X[weight]'.
      This indicates a new AA ("X") with the specified weight, e.g. RX[148.5]T. Note that this tag does not alter the 
      AA's to the left (R) or right (T).  Rather, X represents an AA on its own.
      Be careful when converting AASequence to an EmpiricalFormula using .getFormula(), as tags will not be considered
      in this case (there exists no formula for them). However, they have an influence on .getMonoWeight() and .getAverageWeight()!

      If a string cannot be converted into a valid instance of AASequence, the valid flag is false. The flag
      can be read using the isValid() predicate. However, instances of AASequence which are not valid report
      wrong weights, because the weight cannot be calculated then. Also other operations might fail.

      @ingroup Chemistry
  */
  class OPENMS_DLLAPI AASequence
  {

public:
    class Iterator;

    /** @brief ConstIterator for AASequence

            AASequence constant iterator
    */
    class OPENMS_DLLAPI ConstIterator
    {
public:

      // TODO Iterator constructor for ConstIterator

      typedef const Residue & const_reference;
      typedef Residue & reference;
      typedef const Residue * const_pointer;
      typedef std::vector<const Residue *>::difference_type difference_type;
      typedef Residue value_type;
      typedef const Residue * pointer;
      typedef std::random_access_iterator_tag iterator_category;

      /** @name Constructors and destructors
      */
      //@{
      /// default constructor
      ConstIterator()
      {
      }

      /// detailed constructor with pointer to the vector and offset position
      ConstIterator(const std::vector<const Residue *> * vec_ptr, difference_type position)
      {
        vector_ = vec_ptr;
        position_ = position;
      }

      /// copy constructor
      ConstIterator(const ConstIterator & rhs) :
        vector_(rhs.vector_),
        position_(rhs.position_)
      {
      }

      /// copy constructor from Iterator
      ConstIterator(const AASequence::Iterator & rhs) :
        vector_(rhs.vector_),
        position_(rhs.position_)
      {
      }

      /// destructor
      virtual ~ConstIterator()
      {
      }

      //@}

      /// assignment operator
      ConstIterator & operator=(const ConstIterator & rhs)
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
      const_reference operator*() const
      {
        return *(*vector_)[position_];
      }

      /// dereference operator
      const_pointer operator->() const
      {
        return (*vector_)[position_];
      }

      /// forward jump operator
      const ConstIterator operator+(difference_type diff) const
      {
        return ConstIterator(vector_, position_ + diff);
      }

      difference_type operator-(ConstIterator rhs) const
      {
        return position_ - rhs.position_;
      }

      /// backward jump operator
      const ConstIterator operator-(difference_type diff) const
      {
        return ConstIterator(vector_, position_ - diff);
      }

      /// equality comparator
      bool operator==(const ConstIterator & rhs) const
      {
        return vector_ == rhs.vector_ && position_ == rhs.position_;
      }

      /// inequality operator
      bool operator!=(const ConstIterator & rhs) const
      {
        return vector_ != rhs.vector_ || position_ != rhs.position_;
      }

      /// increment operator
      ConstIterator & operator++()
      {
        ++position_;
        return *this;
      }

      /// decrement operator
      ConstIterator & operator--()
      {
        --position_;
        return *this;
      }

      //@}

protected:

      // pointer to the AASequence vector
      const std::vector<const Residue *> * vector_;

      // position in the AASequence vector
      difference_type position_;
    };


    /** @brief Iterator class for AASequence

            Mutable iterator for AASequence
    */
    class OPENMS_DLLAPI Iterator
    {
public:

      friend class AASequence::ConstIterator;

      typedef const Residue & const_reference;
      typedef Residue & reference;
      typedef const Residue * const_pointer;
      typedef const Residue * pointer;
      typedef std::vector<const Residue *>::difference_type difference_type;

      /** @name Constructors and destructors
      */
      //@{
      /// default constructor
      Iterator()
      {
      }

      /// detailed constructor with pointer to the vector and offset position
      Iterator(std::vector<const Residue *> * vec_ptr, difference_type position)
      {
        vector_ = vec_ptr;
        position_ = position;
      }

      /// copy constructor
      Iterator(const Iterator & rhs) :
        vector_(rhs.vector_),
        position_(rhs.position_)
      {
      }

      /// destructor
      virtual ~Iterator()
      {
      }

      //@}

      /// assignment operator
      Iterator & operator=(const Iterator & rhs)
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
      const_reference operator*() const
      {
        return *(*vector_)[position_];
      }

      /// dereference operator
      const_pointer operator->() const
      {
        return (*vector_)[position_];
      }

      /// mutable dereference operator
      pointer operator->()
      {
        return (*vector_)[position_];
      }

      /// forward jump operator
      const Iterator operator+(difference_type diff) const
      {
        return Iterator(vector_, position_ + diff);
      }

      difference_type operator-(Iterator rhs) const
      {
        return position_ - rhs.position_;
      }

      /// backward jump operator
      const Iterator operator-(difference_type diff) const
      {
        return Iterator(vector_, position_ - diff);
      }

      /// equality comparator
      bool operator==(const Iterator & rhs) const
      {
        return vector_ == rhs.vector_ && position_ == rhs.position_;
      }

      /// inequality operator
      bool operator!=(const Iterator & rhs) const
      {
        return vector_ != rhs.vector_ || position_ != rhs.position_;
      }

      /// increment operator
      Iterator & operator++()
      {
        ++position_;
        return *this;
      }

      /// decrement operator
      Iterator & operator--()
      {
        --position_;
        return *this;
      }

      //@}

protected:

      // pointer to the AASequence vector
      std::vector<const Residue *> * vector_;

      // position in the AASequence vector
      difference_type position_;
    };



    /** @name Constructors and Destructors
    */
    //@{
    /// default constructor
    AASequence();

    /// copy constructor
    AASequence(const AASequence & rhs);

    /// copy constructor from a String
    explicit AASequence(const String & rhs);

    /// copy constructor from char* string
    explicit AASequence(const char * rhs);

    /// destructor
    virtual ~AASequence();
    //@}

    /// assignment operator
    AASequence & operator=(const AASequence & rhs);

    /// check if sequence is empty
    bool empty() const;

    /** @name Accessors
    */
    //@{
    /// returns the peptide as string with modifications embedded in brackets
    String toString() const;

    /// returns the peptide as string without any modifications
    String toUnmodifiedString() const;

    /// set the modification of the residue at position index
    void setModification(Size index, const String & modification);

    /// sets the N-terminal modification
    void setNTerminalModification(const String & modification);

    /// returns the Id of the N-term modification; an empty string is returned if none was set
    const String & getNTerminalModification() const;

    /// sets the C-terminal modification
    void setCTerminalModification(const String & modification);

    /// returns the Id of the C-term modification; an empty string is returned if none was set
    const String & getCTerminalModification() const;

    /// sets the string of the sequence; returns true if the conversion to real AASequence was successful, false otherwise
    bool setStringSequence(const String & sequence);

    /// returns a pointer to the residue, which is at position index
    const Residue & getResidue(SignedSize index) const;

    /// returns a pointer to the residue, which is at position index
    const Residue & getResidue(Size index) const;

    /// returns the formula of the peptide
    EmpiricalFormula getFormula(Residue::ResidueType type = Residue::Full, Int charge = 0) const;

    /// returns the average weight of the peptide
    DoubleReal getAverageWeight(Residue::ResidueType type = Residue::Full, Int charge = 0) const;

    /// returns the mono isotopic weight of the peptide
    DoubleReal getMonoWeight(Residue::ResidueType type = Residue::Full, Int charge = 0) const;

    /// returns a pointer to the residue at given position
    const Residue & operator[](SignedSize index) const;

    /// returns a pointer to the residue at given position
    const Residue & operator[](Size index) const;

    /// adds the residues of the peptide
    AASequence operator+(const AASequence & peptide) const;

    /// adds the residues of the peptide, which is given as a string
    AASequence operator+(const String & peptide) const;

    /// adds the residue of the peptide, which is given as string literal
    AASequence operator+(const char * rhs) const;

    /// adds the residue to the peptide; the residue must be a valid residue of the ResidueDB
    AASequence operator+(const Residue * residue) const;

    /// adds the residues of a peptide
    AASequence & operator+=(const AASequence &);

    /// adds the residues of a peptide, which is given as a string
    AASequence & operator+=(const String &);

    /// adds the residues of a peptide, which is given as string literal
    AASequence & operator+=(const char * rhs);

    /// adds the residue to the peptide; the residue must be a valid residue of the ResidueDB
    AASequence & operator+=(const Residue * residue);

    /// returns the number of residues
    Size size() const;

    /// returns a peptide sequence of the first index residues
    AASequence getPrefix(Size index) const;

    /// returns a peptide sequence of the last index residues
    AASequence getSuffix(Size index) const;

    /// returns a peptide sequence of number residues, beginning at position index
    AASequence getSubsequence(Size index, UInt number) const;

    /// counts the number of occurrences of residue given by a string
    Size getNumberOf(const String & residue) const;

    /// compute frequency table of amino acids
    void getAAFrequencies(Map<String, Size> & frequency_table) const;

    //@}

    /** @name Predicates
    */
    //@{
    /** @brief return true if the instance is valid

            Valid means that a possible given sequence as string was successful
            converted into a real amino acid sequence which meaningful amino acids
            and modifications associated with it.
    */
    bool isValid() const;

    /// returns true if the peptide contains the given residue
    bool has(const Residue & residue) const;

    /// returns true if the peptide contains the given residue
    bool has(const String & name) const;

    /// returns true if the peptide contains the given peptide
    /// @note c-term and n-term mods are ignored
    bool hasSubsequence(const AASequence & peptide) const;

    /// returns true if the peptide contains the given peptide
    /// @note c-term and n-term mods are ignored
    bool hasSubsequence(const String & peptide) const;

    /// returns true if the peptide has the given prefix
    /// n-term mod is also checked (c-term as well, if prefix is of same length)
    bool hasPrefix(const AASequence & peptide) const;

    /// returns true if the peptide has the given prefix
    /// n-term mod is also checked (c-term as well, if prefix is of same length)
    bool hasPrefix(const String & peptide) const;

    /// returns true if the peptide has the given suffix
    /// c-term mod is also checked (n-term as well, if suffix is of same length)
    bool hasSuffix(const AASequence & peptide) const;

    /// returns true if the peptide has the given suffix
    /// c-term mod is also checked (n-term as well, if suffix is of same length)
    bool hasSuffix(const String & peptide) const;

    /// predicate which is true if the peptide is N-term modified
    bool hasNTerminalModification() const;

    /// predicate which is true if the peptide is C-term modified
    bool hasCTerminalModification() const;

    // returns true if any of the residues is modified
    bool isModified() const;

    /// returns true if the residue at the position is modified
    bool isModified(Size index) const;

    /// equality operator
    bool operator==(const AASequence & rhs) const;

    /// equality operator given the peptide as a string
    bool operator==(const String & rhs) const;

    /// equality operator given the peptide as string literal
    bool operator==(const char * rhs) const;

    /// lesser than operator which compares the C-term mods, sequence and N-term mods; can be used for maps
    bool operator<(const AASequence & rhs) const;

    /// inequality operator
    bool operator!=(const AASequence & rhs) const;

    /// inequality operator given the peptide as a string
    bool operator!=(const String & rhs) const;

    /// inequality operator given the peptide as string literal
    bool operator!=(const char * rhs) const;
    //@}

    /** @name Iterators
    */
    //@{
    inline Iterator begin() { return Iterator(&peptide_, 0); }

    inline ConstIterator begin() const { return ConstIterator(&peptide_, 0); }

    inline Iterator end() { return Iterator(&peptide_, (Int) peptide_.size()); }

    inline ConstIterator end() const { return ConstIterator(&peptide_, (Int) peptide_.size()); }
    //@}

    /** @name Stream operators
    */
    //@{
    /// writes a peptide to an output stream
    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const AASequence & peptide);

    /// reads a peptide from an input stream
    friend OPENMS_DLLAPI std::istream & operator>>(std::istream & is, const AASequence & peptide);
    //@}

protected:

    std::vector<const Residue *> peptide_;

    String sequence_string_;

    void parseString_(std::vector<const Residue *> & sequence, const String & peptide);

    ResidueDB * getResidueDB_() const;

    bool valid_;

    const ResidueModification * n_term_mod_;

    const ResidueModification * c_term_mod_;
  };

  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const AASequence & peptide);

  OPENMS_DLLAPI std::istream & operator>>(std::istream & os, const AASequence & peptide);

} // namespace OpenMS

#endif
