// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/WeightWrapper.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

  WeightWrapper::WeightWrapper() :
    weight_mode_(MONO)
  {
  }

  WeightWrapper::WeightWrapper(const WEIGHTMODE weight_mode) :
    weight_mode_(weight_mode)
  {
  }

  WeightWrapper::WeightWrapper(const WeightWrapper & source) = default;

  WeightWrapper::~WeightWrapper() = default;

  void WeightWrapper::setWeightMode(const WEIGHTMODE mode)
  {
    if (mode >= WeightWrapper::SIZE_OF_WEIGHTMODE)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "setWeightMode() received illegal 'mode' value!");
    }
    weight_mode_ = mode;
  }

  WeightWrapper::WEIGHTMODE WeightWrapper::getWeightMode() const
  {
    return weight_mode_;
  }

  double WeightWrapper::getWeight(const AASequence & aa) const
  {
    if (weight_mode_ == WeightWrapper::MONO)
    {
      return aa.getMonoWeight();
    }
    else
    {
      return aa.getAverageWeight();
    }
  }

  double WeightWrapper::getWeight(const EmpiricalFormula & ef) const
  {
    if (weight_mode_ == WeightWrapper::MONO)
    {
      return ef.getMonoWeight();
    }
    else
    {
      return ef.getAverageWeight();
    }
  }

  double WeightWrapper::getWeight(const Residue & r, Residue::ResidueType res_type) const
  {
    if (weight_mode_ == WeightWrapper::MONO)
    {
      return r.getMonoWeight(res_type);
    }
    else
    {
      return r.getAverageWeight(res_type);
    }
  }

}
