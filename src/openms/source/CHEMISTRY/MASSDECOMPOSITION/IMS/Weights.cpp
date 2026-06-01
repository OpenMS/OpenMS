// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/Weights.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS::ims
{

  Weights & Weights::operator=(const Weights & other)
  {
    if (this != &other)
    {
      alphabet_masses_ = other.alphabet_masses_;
      precision_ = other.precision_;
      weights_ = other.weights_;
    }
    return *this;
  }

  void Weights::setPrecision(Weights::alphabet_mass_type precision)
  {
    this->precision_ = precision;
    weights_.clear();
    // convert alphabet masses (double) to integer masses (weights) with the given precision
    for (alphabet_masses_type::size_type i = 0; i < alphabet_masses_.size(); ++i)
    {
      weights_.push_back(static_cast<weight_type>(floor((alphabet_masses_[i] / precision) + 0.5)));
    }
  }

  void Weights::swap(size_type index1, size_type index2)
  {
    weight_type weight = weights_[index1];
    weights_[index1] = weights_[index2];
    weights_[index2] = weight;

    alphabet_mass_type mass = alphabet_masses_[index1];
    alphabet_masses_[index1] = alphabet_masses_[index2];
    alphabet_masses_[index2] = mass;
  }

  Weights::alphabet_mass_type Weights::getParentMass(const std::vector<unsigned int> & decomposition) const
  {
    // checker whether the passed decomposition is applicable
    if (alphabet_masses_.size() != decomposition.size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("The passed decomposition has the wrong size. Expected ") + String(alphabet_masses_.size()) + String(" but got ") + String(decomposition.size()) + String("."));
    }

    alphabet_mass_type parent_mass = 0;

    for (std::vector<unsigned int>::size_type i = 0; i < decomposition.size(); ++i)
    {
      parent_mass += alphabet_masses_[i] * decomposition[i];
    }
    return parent_mass;
  }

  bool Weights::divideByGCD()
  {
    if (weights_.size() < 2)
    {
      return false;
    }
    weight_type d = Math::gcd(weights_[0], weights_[1]);
    for (weights_type::size_type i = 2; i < weights_.size(); ++i)
    {
      d = Math::gcd(d, weights_[i]);
      if (d == 1)
      {
        return false;
      }
    }
    // if we're here: d != 1
    precision_ *= d;

    // rescales the integer weights. Don't use setPrecision() here since
    // the result could be different due to rounding errors.
    for (weights_type::size_type i = 0; i < weights_.size(); ++i)
    {
      weights_[i] /= d;
    }
    return true;
  }

  Weights::alphabet_mass_type Weights::getMinRoundingError() const
  {
    alphabet_mass_type min_error = 0;
    for (size_type i = 0; i < weights_.size(); ++i)
    {
      alphabet_mass_type error = (precision_ * static_cast<alphabet_mass_type>(weights_[i]) - alphabet_masses_[i]) / alphabet_masses_[i];
      if (error < 0 && error < min_error)
      {
        min_error = error;
      }
    }
    return min_error;
  }

  Weights::alphabet_mass_type Weights::getMaxRoundingError() const
  {
    alphabet_mass_type max_error = 0;
    for (size_type i = 0; i < weights_.size(); ++i)
    {
      alphabet_mass_type error = (precision_ * static_cast<alphabet_mass_type>(weights_[i]) - alphabet_masses_[i]) / alphabet_masses_[i];
      if (error > 0 && error > max_error)
      {
        max_error = error;
      }
    }
    return max_error;
  }

  std::ostream & operator<<(std::ostream & os, const Weights & weights)
  {
    for (Weights::size_type i = 0; i < weights.size(); ++i)
    {
      os << weights.getWeight(i) << std::endl;
    }
    return os;
  }

} // namespace OpenMS // namespace ims
