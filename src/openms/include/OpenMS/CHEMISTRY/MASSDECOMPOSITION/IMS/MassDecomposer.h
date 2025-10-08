// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#pragma once

#include <vector>

namespace OpenMS
{

  namespace ims
  {
    /**
      @brief An interface to handle decomposing of integer values/masses
      over a set of integer weights (alphabet).

      An interface that addresses the following "mass decomposition" problems:

      - Existence Problem (whether the decomposition of the given mass exists),
      - One Decomposition Problem (returns one possible decomposition),
      - All Decompositions Problem (returns all possible decompositions),
      - Number of Decompositions Problem (returns the number of possible
       decompositions).

      Those problems are solved in integer arithmetic, i.e. only exact
      solutions are found with no error allowed.

      @param ValueType Type of values to be decomposed.
      @param DecompositionValueType Type of decomposition elements.

      @author Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE>
    */
    template <typename ValueType,
              typename DecompositionValueType>
    class MassDecomposer
    {
public:
      /**
        Type of value to be decomposed.
      */
      typedef ValueType value_type;

      /**
        Type of decomposition value.
      */
      typedef DecompositionValueType decomposition_value_type;

      /**
        Type of decomposition container.
      */
      typedef std::vector<decomposition_value_type> decomposition_type;

      /**
        Type of container for many decompositions.
      */
      typedef std::vector<decomposition_type> decompositions_type;

      /**
        A virtual destructor.
      */
      virtual ~MassDecomposer(){}

      /**
        Returns true if the decomposition for the given @c mass exists, otherwise - false.

        @param mass Mass to be checked on decomposing.
        @return true, if the decomposition for @c mass exist, otherwise - false.
      */
      virtual bool exist(value_type mass) = 0;

      /**
        Returns one possible decomposition of the given @c mass.

        @param mass Mass to be decomposed.
        @return The decomposition of the @c mass, if one exists, otherwise - an empty container.
      */
      virtual decomposition_type getDecomposition(value_type mass) = 0;

      /**
        Returns all possible decompositions for the given @c mass.

        @param mass Mass to be decomposed.
        @return All possible decompositions of the @c mass, if there are any exist,
        otherwise - an empty container.
      */
      virtual decompositions_type getAllDecompositions(value_type mass) = 0;

      /**
        Returns the number of possible decompositions for the given @c mass.
        *
        @param mass Mass to be decomposed.
        @return The number of possible decompositions for the @c mass.
      */
      virtual decomposition_value_type getNumberOfDecompositions(value_type mass) = 0;

    };

  } // namespace ims
} // namespace OpenMS

