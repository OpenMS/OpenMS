// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#pragma once

#include <utility>
#include <map>
#include <memory>

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IntegerMassDecomposer.h>

namespace OpenMS
{
  namespace ims
  {

    /**
      @brief Handles decomposing of non-integer values/masses
      over a set of non-integer weights with an error allowed.

      Implements a decomposition of non-integer values with a certain error
      allowed. Exactness of decomposition can also be tuned by setting a
      precision factor for weights defining their scaling magnitude.

      Works in fact as a wrapper for classes that handle exact mass decomposing
      using integer arithmetics. Instead of decomposing a single value as done
      by integer mass decomposers, @c RealMassDecomposer defines a set of values
      that lie in the allowed range (defined by error and false negatives
      appeared due to rounding), scales those to integers, decomposes
      them using @c IntegerMassDecomposer, does some checks (i.e. on false
      positives appeared due to rounding) and collects decompositions together.

      @author Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE>
    */
    class OPENMS_DLLAPI RealMassDecomposer
    {
public:

      /// Type of integer decomposer.
      typedef IntegerMassDecomposer<> integer_decomposer_type;

      /// Type of integer values that are decomposed.
      typedef integer_decomposer_type::value_type integer_value_type;

      /// Type of result decompositions from integer decomposer.
      typedef integer_decomposer_type::decompositions_type decompositions_type;

      /// Type of the number of decompositions.
      typedef unsigned long long number_of_decompositions_type;

      typedef std::map<unsigned int, std::pair<unsigned int, unsigned int> > constraints_type;

      /**
        Constructor with weights.

        @param weights Weights over which values/masses to be decomposed.
      */
      explicit RealMassDecomposer(const Weights & weights);

      /**
        Gets all decompositions for a @c mass with an @c error allowed.

        @param mass Mass to be decomposed.
        @param error Error allowed between given and result decomposition.
        @return All possible decompositions for a given mass and error.
      */
      decompositions_type getDecompositions(double mass, double error);

      decompositions_type getDecompositions(double mass, double error, const constraints_type & constraints);

      /**
       Gets a number of all decompositions for a @c mass with an @c error
       allowed. It's similar to the @c getDecompositions(double,double) function
       but less space consuming, since doesn't use container to store decompositions.

       @param mass Mass to be decomposed.
       @param error Error allowed between given and result decomposition.
       @return Number of all decompositions for a given mass and error.
      */
      number_of_decompositions_type getNumberOfDecompositions(double mass, double error);

private:
      /// Weights over which values/masses to be decomposed.
      Weights weights_;

      /// Minimal and maximal rounding errors.
      std::pair<double, double> rounding_errors_;

      /// Precision to scale double values to integer
      double precision_;

      /**
        Decomposer to be used for exact decomposing using
        integer arithmetic.
      */
      std::shared_ptr<integer_decomposer_type> decomposer_;
    };

  } // namespace ims
} // namespace OpenMS

