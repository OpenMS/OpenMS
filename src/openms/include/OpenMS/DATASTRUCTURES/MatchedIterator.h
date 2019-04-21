// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <iterator>

namespace OpenMS
{
  
  /**
    @brief Iterates over all elements of a reference container, and finding the closest match in another target container.

    Both containers must be sorted with respect to the comparator.
    The TRAIT template argument (DaTrait or PpmTrait) encodes the distance metric (ppm or Da) and the threshold.
    The iterator always choses the closest matching peak in the target container, if more than one candidate is found in the
    match-window.


    This iterator is much more efficient than iterating over the reference container and calling findNearest() on the target container
    O(n+m) vs. O(n*log(m). Since this container is much more cache-friendly, the actual speedups are even larger.

  */
  template <typename CONT, typename TRAIT>
  class MatchedIterator
  {
    public:
      // define the 5 types required for an iterator. Deriving from std::iterator is deprecated in C++17.
      using iterator_category = std::forward_iterator_tag;
      using value_type = typename CONT::value_type; //< dereferences to an element in the target container
      using difference_type = std::ptrdiff_t;
      using pointer = typename CONT::value_type*;
      using reference = typename CONT::value_type&;

      typedef typename TRAIT DiffType;
      /**
        @brief Constructs a MatchedIterator using two containers. The way a match is found, depends on the TRAIT type (ppm or Da tolerance)

      */
      explicit MatchedIterator(const CONT& ref, const CONT& target, float tolerance)
        : ref_(ref), target_(target), it_ref_(ref.begin()), it_tgt_(target.begin()), tol_(tolerance), is_end_(false)
      {
        advanceTarget_();
      }
    
      /// same as end()
      explicit MatchedIterator()
        : ref_(empty_), target_(empty_), it_ref_(empty_.end()), it_tgt_(empty_.end()), is_end_(true)
      {
      }

      /// assignment operator (default)
      MatchedIterator(const MatchedIterator& rhs) = default;

      bool operator==(const MatchedIterator& rhs) const
      {
        if (this == &rhs) return true;

        if (isEnd_() || rhs.isEnd_())
        {
          return isEnd_() == rhs.isEnd_();
        }

        return (it_ref_ == rhs.it_ref_ &&
                it_tgt_ == rhs.it_tgt_ &&
                &ref_ == &rhs.ref_ &&
                &target_ == &rhs.target_);
      }
      bool operator!=(const MatchedIterator& rhs) const
      {
        return !(*this == rhs);
      }

      const value_type& operator*() const
      {
        return *it_tgt_;
      }
      const value_type& operator->() const
      {
        return *it_tgt_;
      }

      const value_type& curRef() const
      {
        return *it_ref_;
      }

      size_t refIdx() const
      {
        return it_ref_ - ref_.begin();
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
        while (it_ref_ != ref_.end() && it_tgt_ != target_.end())
        {

          double max_dist_dalton = DiffType::allowedTol(tol_, *it_ref_);
          
          // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
          float diff = std::numeric_limits<float>::max();
          do
          {
            auto d = getMZDiffAbs_(*it_ref_, *it_tgt_);
            if (diff >= d) // we need >=, so the very first comparison is true sets 'diff'
            {
              diff = d; // getting better
            }
            else {
              // getting worse (overshot)
              --it_tgt_;
              break;
            }
            ++it_tgt_;
          } while (it_tgt_ != target_.end());

          if (it_tgt_ == target_.end())
          { // reset to last valid entry
            --it_tgt_;
          }
          if (diff < max_dist_dalton) return; // ok, found match
          
          // try next ref peak
          ++it_ref_;
        }
        // reached end of ref or target container
        setToEnd_();
      }

      MatchedIterator& operator++()
      {
        // are we at end already? --> wrong usage
        if (isEnd_()) throw Exception::InvalidIterator(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        ++it_ref_;
        advanceTarget_();
        return *this;
      }

      MatchedIterator operator++(int) const
      {
        MatchedIterator n(*this);
        ++n;
        return n;
      }

      static const MatchedIterator& end()
      {
        const static MatchedIterator it_end;
        return it_end;
      }

    protected:
      const CONT& ref_;
      const CONT& target_;
      typename CONT::const_iterator it_ref_, it_tgt_;
      float tol_;

      static const CONT empty_;
      bool is_end_;

      /// for Peak1D & Co
      template <typename T>
      static float getMZDiffAbs_(const T& mz_ref, const T& mz_obs)
      {
        return fabs(getMZ(mz_ref) - getMZ(mz_obs));
      }
  };

  template <typename CONT, typename TRAIT>
  const CONT MatchedIterator<CONT, TRAIT>::empty_;

  inline float getMZ(float m) { return m;}
  template <typename T>
  float getMZ(const T& o) { return o.getMZ(); }

  struct PpmTrait
  {
    template <typename T>
    static float allowedTol(float tol, const T& mz_ref)
    {
      return Math::ppmToMass(tol, getMZ(mz_ref));
    }
  };

  struct DaTrait
  {
    template <typename T>
    static float allowedTol(float tol, const T& /*mz_ref*/)
    {
      return tol;
    }
  };

} // namespace OpenMS
