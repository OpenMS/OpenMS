// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/RealMassDecomposer.h>
#include <iostream>
#include <memory>

namespace OpenMS::ims
{
  RealMassDecomposer::RealMassDecomposer(const Weights & weights) :
    weights_(weights)
  {

    rounding_errors_ = std::make_pair(weights.getMinRoundingError(), weights.getMaxRoundingError());
    precision_ = weights.getPrecision();
    decomposer_ = std::make_shared<integer_decomposer_type>(
      weights);
  }

  RealMassDecomposer::decompositions_type RealMassDecomposer::getDecompositions(double mass, double error)
  {
    // defines the range of integers to be decomposed
    integer_value_type start_integer_mass = static_cast<integer_value_type>(
      ceil((1 + rounding_errors_.first) * (mass - error) / precision_));
    integer_value_type end_integer_mass = static_cast<integer_value_type>(
      floor((1 + rounding_errors_.second) * (mass + error) / precision_));

    decompositions_type all_decompositions_from_range;

    // loops and finds decompositions for every integer mass,
    // then checks if real mass of decomposition lays in the allowed
    // error interval [mass-error; mass+error]
    for (integer_value_type integer_mass = start_integer_mass;
         integer_mass < end_integer_mass; ++integer_mass)
    {
      decompositions_type decompositions =
        decomposer_->getAllDecompositions(integer_mass);
      for (decompositions_type::iterator pos = decompositions.begin();
           pos != decompositions.end(); )
      {
        double parent_mass = weights_.getParentMass(*pos);
        if (fabs(parent_mass - mass) > error)
        {
          pos = decompositions.erase(pos);
        }
        else
        {
          ++pos;
        }
      }
      all_decompositions_from_range.insert(all_decompositions_from_range.end(),
                                           decompositions.begin(), decompositions.end());
    }
    return all_decompositions_from_range;
  }

  RealMassDecomposer::decompositions_type RealMassDecomposer::getDecompositions(double mass, double error,
                                                                                const constraints_type & constraints)
  {

    // defines the range of integers to be decomposed
    integer_value_type start_integer_mass = static_cast<integer_value_type>(
      ceil((1 + rounding_errors_.first) * (mass - error) / precision_));
    integer_value_type end_integer_mass = static_cast<integer_value_type>(
      floor((1 + rounding_errors_.second) * (mass + error) / precision_));

    decompositions_type all_decompositions_from_range;

    // loops and finds decompositions for every integer mass,
    // then checks if real mass of decomposition lays in the allowed
    // error interval [mass-error; mass+error]
    for (integer_value_type integer_mass = start_integer_mass;
         integer_mass < end_integer_mass; ++integer_mass)
    {
      decompositions_type decompositions =
        decomposer_->getAllDecompositions(integer_mass);
      for (decompositions_type::iterator pos = decompositions.begin();
           pos != decompositions.end(); )
        {
        double parent_mass = weights_.getParentMass(*pos);
        if (fabs(parent_mass - mass) > error)
        {
          pos = decompositions.erase(pos);
        }
        else
        {
          bool to_erase = false;
          if (!constraints.empty())
          {
            for (constraints_type::const_iterator it =
                   constraints.begin(); it != constraints.end(); ++it)
              {
              if ((*pos)[it->first] < it->second.first ||
                  (*pos)[it->first] > it->second.second)
              {
                to_erase = true;
                break;
              }
            }
          }
          if (to_erase)
          {
            pos = decompositions.erase(pos);
          }
          else
          {
            ++pos;
          }
        }
      }
      all_decompositions_from_range.insert(all_decompositions_from_range.end(),
                                           decompositions.begin(), decompositions.end());
    }
    return all_decompositions_from_range;
  }

  RealMassDecomposer::number_of_decompositions_type RealMassDecomposer::getNumberOfDecompositions(double mass, double error)
  {
    // defines the range of integers to be decomposed
    integer_value_type start_integer_mass = static_cast<integer_value_type>(1);
    if (mass - error > 0)
    {
      start_integer_mass = static_cast<integer_value_type>(
        ceil((1 + rounding_errors_.first) * (mass - error) / precision_));
    }
    integer_value_type end_integer_mass = static_cast<integer_value_type>(
      floor((1 + rounding_errors_.second) * (mass + error) / precision_));

    number_of_decompositions_type number_of_decompositions = static_cast<number_of_decompositions_type>(0);

    // loops and finds decompositions for every integer mass,
    // then checks if real mass of decomposition lays in the allowed
    // error interval [mass-error; mass+error]
    for (integer_value_type integer_mass = start_integer_mass;
         integer_mass < end_integer_mass; ++integer_mass)
    {
      decompositions_type decompositions =
        decomposer_->getAllDecompositions(integer_mass);
      for (decompositions_type::iterator pos = decompositions.begin();
           pos != decompositions.end(); ++pos)
      {
        double parent_mass = weights_.getParentMass(*pos);
        if (fabs(parent_mass - mass) <= error)
        {
          ++number_of_decompositions;
        }
      }
    }
    return number_of_decompositions;
  }

} // namespace OpenMS  // namespace ims
