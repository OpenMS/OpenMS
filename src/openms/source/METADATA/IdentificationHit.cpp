// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/IdentificationHit.h>

using namespace std;

namespace OpenMS
{

  IdentificationHit::IdentificationHit() :
    MetaInfoInterface(),
    id_(),
    charge_(0),
    calculated_mass_to_charge_(0.0),
    name_(""),
    pass_threshold_(true),
    rank_(0)
  {
  }

  IdentificationHit::~IdentificationHit() = default;

  // Equality operator
  bool IdentificationHit::operator==(const IdentificationHit & rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && id_ == rhs.id_
           && charge_ == rhs.charge_
           && calculated_mass_to_charge_ == rhs.calculated_mass_to_charge_
           && experimental_mass_to_charge_ == rhs.experimental_mass_to_charge_
           && name_ == rhs.name_
           && pass_threshold_ == rhs.pass_threshold_
           && rank_ == rhs.rank_;
  }

  // Inequality operator
  bool IdentificationHit::operator!=(const IdentificationHit & rhs) const
  {
    return !(*this == rhs);
  }

  void IdentificationHit::setId(const String & id)
  {
    id_ = id;
  }

  const String & IdentificationHit::getId() const
  {
    return id_;
  }

  void IdentificationHit::setCharge(Int charge)
  {
    charge_ = charge;
  }

  Int IdentificationHit::getCharge() const
  {
    return charge_;
  }

  void IdentificationHit::setCalculatedMassToCharge(double mz)
  {
    calculated_mass_to_charge_ = mz;
  }

  double IdentificationHit::getCalculatedMassToCharge() const
  {
    return calculated_mass_to_charge_;
  }

  void IdentificationHit::setExperimentalMassToCharge(double mz)
  {
    experimental_mass_to_charge_ = mz;
  }

  double IdentificationHit::getExperimentalMassToCharge() const
  {
    return experimental_mass_to_charge_;
  }

  void IdentificationHit::setName(const String & name)
  {
    name_ = name;
  }

  const String & IdentificationHit::getName() const
  {
    return name_;
  }

  void IdentificationHit::setPassThreshold(bool pass)
  {
    pass_threshold_ = pass;
  }

  bool IdentificationHit::getPassThreshold() const
  {
    return pass_threshold_;
  }

  void IdentificationHit::setRank(Int rank)
  {
    rank_ = rank;
  }

  Int IdentificationHit::getRank() const
  {
    return rank_;
  }

} // namespace OpenMS

