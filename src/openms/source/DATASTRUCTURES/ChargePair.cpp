// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/ChargePair.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>

#include <ostream>

namespace OpenMS
{

  ChargePair::ChargePair() :
    feature0_index_(0),
    feature1_index_(0),
    feature0_charge_(0),
    feature1_charge_(0),
    compomer_(),
    mass_diff_(0),
    score_(1),
    is_active_(false)
  {
  }

  /// Constructor from map index, element index and Feature
  ChargePair::ChargePair(const Size& index0,
                         const Size& index1,
                         const Int& charge0,
                         const Int& charge1,
                         const Compomer& compomer,
                         const double& mass_diff,
                         const bool active) :
    feature0_index_(index0),
    feature1_index_(index1),
    feature0_charge_(charge0),
    feature1_charge_(charge1),
    compomer_(compomer),
    mass_diff_(mass_diff),
    score_(1),
    is_active_(active)
  {
  }

  /// Copy constructor
  ChargePair::ChargePair(const ChargePair& rhs) = default;

  /// Assignment operator
  ChargePair& ChargePair::operator=(const ChargePair& rhs)
  {
    if (&rhs == this)
    {
      return *this;
    }
    feature0_index_ = rhs.feature0_index_;
    feature1_index_ = rhs.feature1_index_;
    feature0_charge_ = rhs.feature0_charge_;
    feature1_charge_ = rhs.feature1_charge_;
    compomer_ = rhs.compomer_;
    mass_diff_ = rhs.mass_diff_;
    score_ = rhs.score_;
    is_active_ = rhs.is_active_;

    return *this;
  }

  //@}

  //@name Accessors
  //@{
  /// Returns the charge (for element 0 or 1)
  Int ChargePair::getCharge(UInt pairID) const
  {
    if (pairID == 0)
    {
      return feature0_charge_;
    }
    else
    {
      return feature1_charge_;
    }
  }

  /// Set the charge (for element 0 or 1)
  void ChargePair::setCharge(UInt pairID, Int e)
  {
    if (pairID == 0)
    {
      feature0_charge_ = e;
    }
    else
    {
      feature1_charge_ = e;
    }
  }

  /// Returns the element index (for element 0 or 1)
  Size ChargePair::getElementIndex(UInt pairID) const
  {
    if (pairID == 0)
    {
      return feature0_index_;
    }
    else
    {
      return feature1_index_;
    }
  }

  /// Set the element index (for element 0 or 1)
  void ChargePair::setElementIndex(UInt pairID, Size e)
  {
    if (pairID == 0)
    {
      feature0_index_ = e;
    }
    else
    {
      feature1_index_ = e;
    }
  }

  /// Returns the Id of the compomer that explains the mass difference
  const Compomer& ChargePair::getCompomer() const
  {
    return compomer_;
  }

  /// Set the compomer id
  void ChargePair::setCompomer(const Compomer& compomer)
  {
    compomer_ = compomer;
  }

  /// Returns the mass difference
  double ChargePair::getMassDiff() const
  {
    return mass_diff_;
  }

  /// Sets the mass difference
  void ChargePair::setMassDiff(double mass_diff)
  {
    mass_diff_ = mass_diff;
  }

  /// Returns the ILP edge score
  double ChargePair::getEdgeScore() const
  {
    return score_;
  }

  /// Sets the ILP edge score
  void ChargePair::setEdgeScore(double score)
  {
    score_ = score;
  }

  /// is this pair realized?
  bool ChargePair::isActive() const
  {
    return is_active_;
  }

  void ChargePair::setActive(const bool active)
  {
    is_active_ = active;
  }

  //@}

  /// Equality operator
  bool ChargePair::operator==(const ChargePair& i) const
  {
    return (feature0_index_ == i.feature0_index_) &&
           (feature1_index_ == i.feature1_index_) &&
           (feature0_charge_ == i.feature0_charge_) &&
           (feature1_charge_ == i.feature1_charge_) &&
           (compomer_ == i.compomer_) &&
           (mass_diff_ == i.mass_diff_) &&
           (is_active_ == i.is_active_);
  }

  /// Equality operator
  bool ChargePair::operator!=(const ChargePair& i) const
  {
    return !(this->operator==(i));
  }

  std::ostream& operator<<(std::ostream& os, const ChargePair& cp)
  {
    os << "---------- ChargePair -----------------\n"
       << "Mass Diff: " << cp.getMassDiff() << "\n"
       << "Compomer: " << cp.getCompomer() << "\n"
       << "Charge: " << cp.getCharge(0) << " : " << cp.getCharge(1) << "\n"
       << "Element Index: " << cp.getElementIndex(0) << " : " << cp.getElementIndex(1) << "\n";
    return os;
  }

}
