// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_EUCLIDEANSIMILARITY_H
#define OPENMS_COMPARISON_CLUSTERING_EUCLIDEANSIMILARITY_H

#include <cmath>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

  /**
      @brief CompareFunctor for 2Dpoints

      each 2D point as a pair of Real holds a Real coordinate for each Dimension

      @ingroup DummyComparison
  */

  class OPENMS_DLLAPI EuclideanSimilarity
  {
private:

    Real scale_;

public:

    /// default constructor
    EuclideanSimilarity();

    /// copy constructor
    EuclideanSimilarity(const EuclideanSimilarity & source);

    /// destructor
    virtual ~EuclideanSimilarity();

    /// assignment operator
    EuclideanSimilarity & operator=(const EuclideanSimilarity & source);


    /**
        @brief calculates similarity between two points in euclidean space

        @param a a pair of Real, giving the x and the y coordinates of the first point
        @param b a pair of Real, giving the x and the y coordinates of the second point

        calculates similarity from the euclidean distance between given 2D points, scaled in [0,1] @see setScale
    */
    Real operator()(const std::pair<Real, Real> & a, const std::pair<Real, Real> & b) const;

    /**
        @brief calculates self similarity, will yield 0

        @param c a pair of Real, giving the x and the y coordinates

    */
    Real operator()(const std::pair<Real, Real> & c) const;

    /**
        @brief clusters the indices according to their respective element distances

        @param x Real value to scale the result
        @throw Exception::DivisionByZero if scaling is unapplicable because it is 0

        sets the scale so that similarities can be correctly calculated from distances. Should be set so that the greatest distance in a chosen set will be scales to 1 (i.e. @p x = greatest possible distance in the set)
*/
    void setScale(Real x);

  };

}
#endif //OPENMS_COMPARISON_CLUSTERING_EUCLIDEANSIMILARITY_H
