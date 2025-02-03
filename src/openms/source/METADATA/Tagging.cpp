// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Tagging.h>

using namespace std;

namespace OpenMS
{
  const std::string Tagging::NamesOfIsotopeVariant[] = {"LIGHT", "HEAVY"};

  Tagging::Tagging() :
    Modification(),
    mass_shift_(0.0),
    variant_(LIGHT)
  {
    type_ = "Tagging";
  }

  Tagging::~Tagging() = default;

  bool Tagging::operator==(const SampleTreatment & rhs) const
  {
    if (type_ != rhs.getType())
    {
      return false;
    }

    const Tagging * tmp = dynamic_cast<const Tagging *>(&rhs);
    return Modification::operator==(rhs)
    && mass_shift_ == tmp->mass_shift_
    && variant_ == tmp->variant_;
  }

  SampleTreatment * Tagging::clone() const
  {
    SampleTreatment * tmp = new Tagging(*this);
    return tmp;
  }

  double Tagging::getMassShift() const
  {
    return mass_shift_;
  }

  void Tagging::setMassShift(double mass_shift)
  {
    mass_shift_ = mass_shift;
  }

  const Tagging::IsotopeVariant & Tagging::getVariant() const
  {
    return variant_;
  }

  void Tagging::setVariant(const Tagging::IsotopeVariant & variant)
  {
    variant_ = variant;
  }

}

