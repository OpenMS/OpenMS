// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein, Timo Sachsenberg $
// --------------------------------------------------------------------------

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
      ConstIterator(const std::vector<const Ribonucleotide*>* vec_ptr,
                    difference_type position)
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
        return (std::tie(vector_, position_) ==
                std::tie(rhs.vector_, rhs.position_));
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
      Iterator(std::vector<const Ribonucleotide*>* vec_ptr,
               difference_type position)
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
        return (std::tie(vector_,position_) ==
                std::tie(rhs.vector_, rhs.position_));
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
    NASequence() = default; /// default constructor
    NASequence(const NASequence&) = default; // Copy constructor
    NASequence(NASequence&&) = default; // Move constructor
    NASequence& operator=(const NASequence&) & = default;  // Copy assignment operator
    NASequence& operator=(NASequence&&) & = default; // Move assignment operator

    /// full constructor
    NASequence(std::vector<const Ribonucleotide*> s,
               const RibonucleotideChainEnd* five_prime,
               const RibonucleotideChainEnd* three_prime);

    virtual ~NASequence() = default; // destructor

    bool operator==(const NASequence& rhs) const; // element-wise equality
    bool operator!=(const NASequence& rhs) const; // not quality
    bool operator<(const NASequence& rhs) const; // less operator

    // getter / setter for sequence
    void setSequence(const std::vector<const Ribonucleotide*>& s);

    const std::vector<const Ribonucleotide*>& getSequence() const
    {
      return s_;
    }

    std::vector<const Ribonucleotide*>& getSequence()
    {
      return s_;
    }

    // getter / setter for ribonucleotide elements (easy wrapped using pyOpenMS)
    void set(size_t index, const Ribonucleotide* r);

    const Ribonucleotide* get(size_t index)
    {
      return s_[index];
    }

    // getter / setter for sequence elements (C++ container style)
    inline const Ribonucleotide*& operator[](size_t index)
    {
      return s_[index];
    }

    inline const Ribonucleotide* const& operator[](size_t index) const
    {
      return s_[index];
    }

    bool empty() const;
    size_t size() const;
    void clear();

    // 5' and 3' modifications
    bool hasFivePrimeModification() const;
    void setFivePrimeModification(const RibonucleotideChainEnd* r);
    const RibonucleotideChainEnd* getFivePrimeModification() const;
    bool hasThreePrimeModification() const;
    void setThreePrimeModification(const RibonucleotideChainEnd* r);
    const RibonucleotideChainEnd* getThreePrimeModification() const;

    // iterators
    inline Iterator begin()
    {
      return Iterator(&s_, 0);
    }

    inline ConstIterator begin() const
    {
      return ConstIterator(&s_, 0);
    }

    inline Iterator end()
    {
      return Iterator(&s_, (Int) s_.size());
    }

    inline ConstIterator end() const
    {
      return ConstIterator(&s_, (Int) s_.size());
    }

    inline ConstIterator cbegin() const
    {
      return ConstIterator(&s_, 0);
    }

    inline ConstIterator cend() const
    {
      return ConstIterator(&s_, (Int) s_.size());
    }

    // utility functions
    double getMonoWeight(
      Ribonucleotide::RibonucleotideFragmentType type = Ribonucleotide::Full,
      Int charge = 0) const;
    EmpiricalFormula getFormula(
      Ribonucleotide::RibonucleotideFragmentType type = Ribonucleotide::Full,
      Int charge = 0) const;
    NASequence getPrefix(Size index) const;
    NASequence getSuffix(Size index) const;

    /**
       @brief create NASequence object by parsing an OpenMS string

       @param s Input string

       @throws Exception::ParseError if an invalid string representation of an AA sequence is passed
    */
    static NASequence fromString(const String& s,
                                 Ribonucleotide::NucleicAcidType type);

    /** @name Stream operators
        writes a NASequence to an output stream
    */
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os,
                                                  const NASequence& seq);

    /**
       @brief create NASequence object by parsing a C string (character array)

       @param s Input string

       @throws Exception::ParseError if an invalid string representation of an AA sequence is passed
    */
    static NASequence fromString(const char* s,
                                 Ribonucleotide::NucleicAcidType type);

    std::string toString() const ;

  private:
    //TODO: query RNA / DNA depending on type
    static void parseString_(const String& s, NASequence& nss,
                             Ribonucleotide::NucleicAcidType type);

    /**
       @brief Parses modifications in square brackets

       @param str_it Current position in the string to be parsed
       @param str Full input string
       @param aas Current AASequence object (will be modified with the correct ribo added)

       @return Position at which to continue parsing
    */
    //TODO: query RNA / DNA depending on type
    static String::ConstIterator parseModSquareBrackets_(
      const String::ConstIterator str_it, const String& str, NASequence& nss,
      Ribonucleotide::NucleicAcidType type);

    std::vector<const Ribonucleotide*> s_;

    const RibonucleotideChainEnd* five_prime_ = nullptr;
    const RibonucleotideChainEnd* three_prime_ = nullptr;
  };

}

#endif // OPENMS_CHEMISTRY_NASEQUENCE_H
