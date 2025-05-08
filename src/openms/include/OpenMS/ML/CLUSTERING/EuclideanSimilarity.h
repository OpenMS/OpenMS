// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <cmath>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

  /**
      @brief CompareFunctor for 2Dpoints

      each 2D point as a pair of float holds a float coordinate for each Dimension

  */

  class OPENMS_DLLAPI EuclideanSimilarity
  {
private:

    float scale_;

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

        @param a a pair of float, giving the x and the y coordinates of the first point
        @param b a pair of float, giving the x and the y coordinates of the second point

        calculates similarity from the euclidean distance between given 2D points, scaled in [0,1] @see setScale
    */
    float operator()(const std::pair<float, float> & a, const std::pair<float, float> & b) const;

    /**
        @brief calculates self similarity, will yield 0

        @param c a pair of float, giving the x and the y coordinates

    */
    float operator()(const std::pair<float, float> & c) const;

    /**
        @brief clusters the indices according to their respective element distances

        @param x float value to scale the result
        @throw Exception::DivisionByZero if scaling is inapplicable because it is 0

        sets the scale so that similarities can be correctly calculated from distances. Should be set so that the greatest distance in a chosen set will be scales to 1 (i.e. @p x = greatest possible distance in the set)
*/
    void setScale(float x);

  };

}
