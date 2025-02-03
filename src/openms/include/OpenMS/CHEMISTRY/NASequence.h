// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein, Timo Sachsenberg, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <iosfwd>
#include <vector>

namespace OpenMS
{
  /**
   * @brief Representation of a nucleic acid sequence
   *
   * NASequence represents nucleic acid sequences (RNA) in %OpenMS. Each NASequence consists
   * of a vector of pointers to Ribonucleotides as well as RibonucleotideChainEnds representing
   * the 5' and 3' ends of the sequence. Each Ribonucleotide has only a single instance.
   * These are accessible through RibonucleotideDB. Modified Ribonucleotides are included in RibonucleotideDB
   * and are expressed as the Modomics Short name surrounded by brackets when converted to string.
   *
   * @ingroup Chemistry
   */


  class OPENMS_DLLAPI NASequence
  {
    /**
      @brief an enum of all possible fragment ion types
      */

  public:
    enum NASFragmentType
    {                //< NB: Not all fragments types are valid for all residue types, this class should probably get split
      Full = 0,      ///< with N-terminus and C-terminus
      Internal,      ///< internal, without any termini
      FivePrime,     ///< only 5' terminus
      ThreePrime,    ///< only 3' terminus
      AIon,          ///< MS:1001229 N-terminus up to the C-alpha/carbonyl carbon bond
      BIon,          ///< MS:1001224 N-terminus up to the peptide bond
      CIon,          ///< MS:1001231 N-terminus up to the amide/C-alpha bond
      XIon,          ///< MS:1001228 amide/C-alpha bond up to the C-terminus
      YIon,          ///< MS:1001220 peptide bond up to the C-terminus
      ZIon,          ///< MS:1001230 C-alpha/carbonyl carbon bond
      Precursor,     ///< MS:1001523 Precursor ion
      BIonMinusH20,  ///< MS:1001222 b ion without water
      YIonMinusH20,  ///< MS:1001223 y ion without water
      BIonMinusNH3,  ///< MS:1001232 b ion without ammonia
      YIonMinusNH3,  ///< MS:1001233 y ion without ammonia
      NonIdentified, ///< MS:1001240 Non-identified ion
      Unannotated,   ///< no stored annotation
      WIon,          ///< W ion, added for nucleic acid support
      AminusB,       ///< A ion with base loss, added for nucleic acid support
      DIon,          ///< D ion, added for nucleic acid support
      SizeOfNASFragmentType
    };

    using ConstRibonucleotidePtr = const Ribonucleotide*;

    class Iterator;

    /** @brief ConstIterator of NASequence class.
     *
     * References to the pointers are returned dereferenced.
     * So we don't need to write (*iterator)->getCode(), but can simply use iterator->getCode().
     */
    class OPENMS_DLLAPI ConstIterator
    {
    public:
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
      ConstIterator(const ConstIterator& rhs) : vector_(rhs.vector_), position_(rhs.position_)
      {
      }

      /// copy constructor from Iterator
      ConstIterator(const NASequence::Iterator& rhs) : vector_(rhs.vector_), position_(rhs.position_)
      {
      }

      /// destructor
      virtual ~ConstIterator()
      {
      }

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
      bool operator==(const ConstIterator& rhs) const
      {
        return (std::tie(vector_, position_) == std::tie(rhs.vector_, rhs.position_));
      }

      /// inequality operator
      bool operator!=(const ConstIterator& rhs) const
      {
        return !(operator==(rhs));
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


    /** @brief Iterator of NASequence class.
     *
     * References to the pointers are returned dereferenced.
     * So we don't need to write (*iterator)->getCode(), but can simply use iterator->getCode().
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
      Iterator(const Iterator& rhs) : vector_(rhs.vector_), position_(rhs.position_)
      {
      }

      /// destructor
      virtual ~Iterator()
      {
      }

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
      bool operator==(const Iterator& rhs) const
      {
        return (std::tie(vector_, position_) == std::tie(rhs.vector_, rhs.position_));
      }

      /// inequality operator
      bool operator!=(const Iterator& rhs) const
      {
        return !this->operator==(rhs);
      }

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
    /*
     * Default constructors and assignment operators.
     */
    NASequence() = default;                               /// default constructor
    NASequence(const NASequence&) = default;              ///< Copy constructor
    NASequence(NASequence&&) = default;                   ///< Move constructor
    NASequence& operator=(const NASequence&) & = default; ///< Copy assignment operator
    NASequence& operator=(NASequence&&) & = default;      ///< Move assignment operator

    /// full constructor
    NASequence(std::vector<const Ribonucleotide*> s, const RibonucleotideChainEnd* five_prime, const RibonucleotideChainEnd* three_prime);

    virtual ~NASequence() = default; /// destructor

    bool operator==(const NASequence& rhs) const; ///< element-wise equality
    bool operator!=(const NASequence& rhs) const; ///< not quality
    bool operator<(const NASequence& rhs) const;  ///< less operator

    /// getter / setter for sequence
    void setSequence(const std::vector<const Ribonucleotide*>& seq);

    const std::vector<const Ribonucleotide*>& getSequence() const
    {
      return seq_;
    }

    std::vector<const Ribonucleotide*>& getSequence()
    {
      return seq_;
    }

    /// getter / setter for ribonucleotide elements (easily wrapped using pyOpenMS)
    void set(size_t index, const Ribonucleotide* r);

    const Ribonucleotide* get(size_t index)
    {
      return seq_[index];
    }

    /// getter / setter for sequence elements (C++ container style)
    inline const Ribonucleotide*& operator[](size_t index)
    {
      return seq_[index];
    }

    inline const Ribonucleotide* const& operator[](size_t index) const
    {
      return seq_[index];
    }

    bool empty() const;
    size_t size() const;
    void clear();

    /// 5' and 3' modifications
    bool hasFivePrimeMod() const;
    void setFivePrimeMod(const RibonucleotideChainEnd* r);
    const RibonucleotideChainEnd* getFivePrimeMod() const;
    bool hasThreePrimeMod() const;
    void setThreePrimeMod(const RibonucleotideChainEnd* r);
    const RibonucleotideChainEnd* getThreePrimeMod() const;

    /// iterators
    inline Iterator begin()
    {
      return Iterator(&seq_, 0);
    }

    inline ConstIterator begin() const
    {
      return ConstIterator(&seq_, 0);
    }

    inline Iterator end()
    {
      return Iterator(&seq_, (Int)seq_.size());
    }

    inline ConstIterator end() const
    {
      return ConstIterator(&seq_, (Int)seq_.size());
    }

    inline ConstIterator cbegin() const
    {
      return ConstIterator(&seq_, 0);
    }

    inline ConstIterator cend() const
    {
      return ConstIterator(&seq_, (Int)seq_.size());
    }

    /// utility functions

    /**
     * @brief Get the Monoisotopic Weight of a NASequence. NB returns the uncharged mass + or - proton masses to match the charge param
     *
     * @param type fragment type to return
     * @param charge how many protons to add or subtract, NB the mass returned is the UNCHARGED MASS
     *
     * @return double
     */
    double getMonoWeight(NASFragmentType type = Full, Int charge = 0) const;

    /**
     * @brief Get the Average Weight of a NASequence. NB returns the uncharged mass + or - proton masses to match the charge param
     *
     * @param type fragment type to return
     * @param charge how many protons to add or subtract, NB the mass returned is the UNCHARGED MASS
     *
     * @return double
     */
    double getAverageWeight(NASFragmentType type = Full, Int charge = 0) const;

    /**
     * @brief Get the formula for a NASequence
     *
     * @param type fragment type for formula
     * @param charge how many H to add or subtract
     *
     * @return EmpiricalFormula
     */
    EmpiricalFormula getFormula(NASFragmentType type = Full, Int charge = 0) const;

    /**
     * @brief Return sequence prefix of the given length (not end index!)
     *
     * @param length

     * @return NASequence
     */
    NASequence getPrefix(Size length) const;

    /**
     * @brief Return sequence suffix of the given length (not start index!)
     *
     * @param length

     * @return NASequence
     */
    NASequence getSuffix(Size length) const;

    /**
     * @brief Return subsequence with given starting position and length
     *
     * @param start
     * @param length
     *
     * @return NASequence
     */
    NASequence getSubsequence(Size start = 0, Size length = Size(-1)) const;

    /**
       @brief create NASequence object by parsing an OpenMS string

       @param s Input string

       @throws Exception::ParseError if an invalid string representation of a nucleic acid sequence is passed
    */
    static NASequence fromString(const String& s);

    /** @name Stream operators
        writes a NASequence to an output stream
    */
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const NASequence& seq);

    /**
       @brief create NASequence object by parsing a C string (character array)

       @param s Input string

       @throws Exception::ParseError if an invalid string representation of a nucleic acid sequence is passed
    */
    static NASequence fromString(const char* s);

    std::string toString() const;

  private:
    // TODO: query RNA / DNA depending on type
    static void parseString_(const String& s, NASequence& nas);

    /**
       @brief Parses modifications in square brackets

       @param str_it Current position in the string to be parsed
       @param str Full input string
       @param nas Current AASequence object (will be modified with the correct ribo added)

       @return Position at which to continue parsing
    */
    // TODO: query RNA / DNA depending on type
    static String::ConstIterator parseMod_(const String::ConstIterator str_it, const String& str, NASequence& nas);

    std::vector<const Ribonucleotide*> seq_;

    const RibonucleotideChainEnd* five_prime_ = nullptr;
    const RibonucleotideChainEnd* three_prime_ = nullptr;
  };

} // namespace OpenMS
