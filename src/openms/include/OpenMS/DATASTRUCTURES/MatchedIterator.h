// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <iterator>
#include <type_traits>

namespace OpenMS
{
  
  /**
    @brief For each element in the reference container the closest peak in the target will be searched. If no match is found within the tolerance window, the peak will be skipped over.

    This class can be used for example to iterate through the matching peaks in two spectra (e.g. experimental spectrum and reference spectrum) that are 
    within a given tolerance (in m/z, RT, or something user-defined).

    The iterator always chooses the closest matching peak in the target container, if more than one candidate is found in the
    match-window. If two peaks have equal distance, the smaller value is preferred.
    If no peak is found within the given tolerance (distance), the reference peak does not yield a result and the next reference peak is tested.
    This means the operator++ can be called at most(!) ref.size()-1 times before it == is.end().
    
    The TRAIT template argument (e.g., ValueTrait, DaTrait or PpmTrait) encodes the distance metric (on the value directly, or a member of the value_type, 
    e.g. ppm or Da for m/z, or RT or any other metric you like).
    The simplest use case would be a vector<float> or similar.
    The TRAIT struct just defines two functions which return some distance metrics and accept elements of the container CONT_T as arguments.
    Both containers must be sorted with respect to the comparator used in TRAIT.

    The CONST_T template argument switches between a 'const_iterator' and an 'iterator'.
        
    This iterator is much more efficient than iterating over the reference container and calling findNearest(), i.e. binary search on the target container,
    i.e. O(n+m) vs. O(n*log(m)). Since this container is much more cache-friendly, the actual speedups are even larger.

  */
  template <typename CONT_T, typename TRAIT, bool CONST_T = true >
  class MatchedIterator
  {
    public:
      // define the 5 types required for an iterator. Deriving from std::iterator is deprecated in C++17.
      using iterator_category = std::forward_iterator_tag;
      using value_type = typename CONT_T::value_type; //< dereferences to an element in the target container
      using difference_type = std::ptrdiff_t;
      using pointer = typename std::conditional<CONST_T, typename CONT_T::value_type const*, typename CONT_T::value_type*>::type;
      using reference = typename std::conditional<CONST_T, typename CONT_T::value_type const&, typename CONT_T::value_type&>::type;
      
      using CONT_IT = typename std::conditional<CONST_T, typename CONT_T::const_iterator, typename CONT_T::iterator>::type; // for dereferencing etc
      using CONST_CONT_IT = typename CONT_T::const_iterator; // for input containers
      
      /**
        @brief Constructs a MatchedIterator on two containers. The way a match is found, depends on the TRAIT type (ppm or Da tolerance)

        For each element in the reference container the closest peak in the target will be searched. If no match is found within the tolerance window, the peak will be skipped over.

        @param ref For each element in this reference container the closest peak in the target will be searched
        @param target Target container
        @param tolerance Maximal distance between a valid matching pair in reference and target (unit is according to TRAIT::getDiffAbsolute(), i.e. could be ppm, Da, seconds, ...)
      */
      explicit MatchedIterator(const CONT_T& ref, const CONT_T& target, float tolerance)
        : MatchedIterator(ref.begin(), ref.end(), target.begin(), target.end(), tolerance)
      {
      }
    
      /**
      @brief Constructs a MatchedIterator on two containers. The way a match is found, depends on the TRAIT type (ppm or Da tolerance)

      For each element in the reference container the closest peak in the target will be searched. If no match is found within the tolerance window, the peak will be skipped over.

      @param ref_begin Begin range of reference container
      @param ref_end End range of reference container
      @param tgt_begin Begin range of reference container
      @param tgt_end End range of reference container
      @param tolerance Maximal distance between a valid matching pair in reference and target (unit is according to TRAIT::getDiffAbsolute(), i.e. could be ppm, Da, seconds, ...)
      */
      explicit MatchedIterator(const CONST_CONT_IT ref_begin, const CONST_CONT_IT ref_end,
                               const CONST_CONT_IT tgt_begin, const CONST_CONT_IT tgt_end,
                               float tolerance)
        : ref_begin_(ref_begin), ref_end_(ref_end), tgt_begin_(tgt_begin), tgt_end_(tgt_end), it_ref_(ref_begin), it_tgt_(tgt_begin), tol_(tolerance), is_end_(false)
      {
        if (tgt_begin_ == tgt_end_)
        { // nothing to iterate over in target (if ref_ were empty, isEnd_() is automatically true)
          setToEnd_();
          return;
        }
        advanceTarget_();
      }

      /**
          @brief Default CTor; do not use this for anything other than assigning to it;
      */ 
      explicit MatchedIterator()
        : ref_begin_(), ref_end_(), tgt_begin_(), tgt_end_(), it_ref_(), it_tgt_(), tol_(), is_end_(false)
      {
      }

      /// Copy CTor (default)
      MatchedIterator(const MatchedIterator& rhs) = default;
      
      /// Assignment operator (default)
      MatchedIterator& operator=(const MatchedIterator& rhs) = default;

      bool operator==(const MatchedIterator& rhs) const
      {
        if (this == &rhs) return true;

        if (isEnd_() || rhs.isEnd_())
        {
          return isEnd_() == rhs.isEnd_();
        }

        return (it_ref_ == rhs.it_ref_ &&
                it_tgt_ == rhs.it_tgt_ &&
                ref_begin_ == rhs.ref_begin_ &&
                ref_end_ == rhs.ref_end_ &&
                tgt_begin_ == rhs.tgt_begin_ &&
                tgt_end_ == rhs.tgt_end_);
      }
      bool operator!=(const MatchedIterator& rhs) const
      {
        return !(*this == rhs);
      }

      /// dereference current target element
      template< bool _CONST = CONST_T >
      typename std::enable_if< _CONST, reference >::type operator*() const
      {
        return *it_tgt_;
      }
      template< bool _CONST = CONST_T >
      typename std::enable_if< !_CONST, reference >::type operator*()
      {
        return *it_tgt_;
      }

      /// pointer to current target element
      template< bool _CONST = CONST_T >
      typename std::enable_if< _CONST, pointer >::type operator->() const
      {
        return &(*it_tgt_);
      }
      template< bool _CONST = CONST_T >
      typename std::enable_if< !_CONST, pointer >::type operator->()
      {
        return &(*it_tgt_);
      }

      /// current element in reference container
      const value_type& ref() const
      {
        return *it_ref_;
      }

      /// index into reference container
      size_t refIdx() const
      {
        return it_ref_ - ref_begin_;
      }
      
      /// index into target container
      size_t tgtIdx() const
      {
        return it_tgt_ - tgt_begin_;
      }

      /**
        @brief Advances to the next valid pair

        @exception Exception::InvalidIterator If iterator is already at end
      */
      MatchedIterator& operator++()
      {
        // are we at end already? --> wrong usage
        OPENMS_PRECONDITION(!isEnd_(), "Tried to advance beyond end iterator!");
        ++it_ref_;
        advanceTarget_();
        return *this;
      }

      /// post-increment
      MatchedIterator operator++(int)
      {
        MatchedIterator n(*this);
        ++(*this);
        return n;
      }
      
      /// the end iterator
      static MatchedIterator end()
      {
        return MatchedIterator(true);
      }

    protected:

      /// protected Ctor which creates and end() iterator
      /// the bool argument is only used to call the correct target (its value is ignored)
      MatchedIterator(bool /*is_end*/)
        : ref_begin_(), ref_end_(), tgt_begin_(), tgt_end_(), it_ref_(), it_tgt_(), tol_(), is_end_(true)
      {
      }

      void setToEnd_()
      {
        is_end_ = true;
      }

      bool isEnd_() const
      {
        return is_end_;
      }

      void advanceTarget_()
      {
        while (it_ref_ != ref_end_)
        { // note: it_tgt_ always points to a valid element (unless the whole container was empty -- see CTor)

          double max_dist = TRAIT::allowedTol(tol_, *it_ref_);

          // forward iterate over elements in target data until distance gets worse
          float diff = std::numeric_limits<float>::max();
          do
          {
            auto d = TRAIT::getDiffAbsolute(*it_ref_, *it_tgt_);
            if (diff > d) // getting better
            {
              diff = d;
            }
            else   // getting worse (overshot)
            {
              --it_tgt_;
              break;
            }
            ++it_tgt_;
          } while (it_tgt_ != tgt_end_);

          if (it_tgt_ == tgt_end_)
          { // reset to last valid entry
            --it_tgt_;
          }
          if (diff <= max_dist) return; // ok, found match

          // try next ref peak
          ++it_ref_;
        }
        // reached end of ref container
        setToEnd_();
        // i.e. isEnd() is true now
      }
      
      CONST_CONT_IT ref_begin_, ref_end_;
      CONST_CONT_IT tgt_begin_, tgt_end_; 
      CONT_IT it_ref_, it_tgt_;
      float tol_;
      bool is_end_ = false;
  };
  
  /// Trait for MatchedIterator to find pairs with a certain distance, which is computed directly on the value_type of the container
  struct ValueTrait
  {
    template <typename T>
    static float allowedTol(float tol, const T& /*mz_ref*/)
    {
      return tol;
    }
    /// just use fabs on the value directly
    template <typename T>
    static T getDiffAbsolute(const T& elem_ref, const T& elem_tgt)
    {
      return fabs(elem_ref - elem_tgt);
    }
  };

  /// Trait for MatchedIterator to find pairs with a certain ppm distance in m/z.
  /// Requires container elements to support .getMZ() as member function
  struct PpmTrait
  {
    template <typename T>
    static float allowedTol(float tol, const T& elem_ref)
    {
      return Math::ppmToMass(tol, (float)elem_ref.getMZ());
    }
    /// for Peak1D & Co
    template <typename T>
    static float getDiffAbsolute(const T& elem_ref, const T& elem_tgt)
    {
      return fabs(elem_ref.getMZ() - elem_tgt.getMZ());
    }
  };

  /// Trait for MatchedIterator to find pairs with a certain Th/Da distance in m/z.
  /// Requires container elements to support .getMZ() as member function
  struct DaTrait
  {
    template <typename T>
    static float allowedTol(float tol, const T& /*mz_ref*/)
    {
      return tol;
    }
    /// for Peak1D & Co
    template <typename T>
    static float getDiffAbsolute(const T& elem_ref, const T& elem_tgt)
    {
      return fabs(elem_ref.getMZ() - elem_tgt.getMZ());
    }
  };

} // namespace OpenMS
