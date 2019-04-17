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
#include <iterator>

namespace OpenMS
{

  /**
    @brief Iterates over all elements of a reference container, and finding the closest match in another target container.

    Both containers must be sorted with respect to the comparator.
    This iterator is much more efficient than iterating over the reference container and calling findNearest() on the target container
    O(n+m) vs. O(n*log(m). Since this container is much more cache-friendly, the actual speedups are even larger.
    
    

  */
  template <typename CONT, typename L_DIFF, typename L_MATCH>
  class MatchedIterator
  {
    public:
      // define the 5 types required for an iterator. Deriving from std::iterator is deprecated in C++17.
      using iterator_category = std::forward_iterator_tag;
      using value_type = CONT::value_type; //< dereferences to an element in the target container
      using difference_type = std::ptrdiff_t;
      using pointer = int*;
      using reference = int&;

      /**
        @brief Constructs a MatchedIterator using two containers, a distance-lambda and a match-lambda function.

        The distance lambda must return an absolute value, the match-lambda a bool (true if the two items are close enough)

      */
      MatchedIterator(const CONT& ref, const CONT& target, const L_DIFF diff, const L_MATCH)
      {
        
      }


    protected:
  };

} // namespace OpenMS
