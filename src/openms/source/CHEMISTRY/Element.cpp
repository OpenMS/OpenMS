// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <ostream>
#include <algorithm>
#include <cassert>

using namespace std;

namespace OpenMS
{
  Element::Element() :
    name_(OPENMS_CHEMISTRY_ELEMENT_NAME_DEFAULT),
    symbol_(OPENMS_CHEMISTRY_ELEMENT_SYMBOL_DEFAULT),
    atomic_number_(OPENMS_CHEMISTRY_ELEMENT_ATOMICNUMBER_DEFAULT),
    average_weight_(OPENMS_CHEMISTRY_ELEMENT_WEIGHT_DEFAULT),
    mono_weight_(OPENMS_CHEMISTRY_ELEMENT_WEIGHT_DEFAULT)
  {
  }

  Element::Element(const Element & e) = default;

  Element::Element(const string & name,
                   const string & symbol,
                   unsigned int atomic_number,
                   double average_weight,
                   double mono_weight,
                   const IsotopeDistribution & isotopes) :
    name_(name),
    symbol_(symbol),
    atomic_number_(atomic_number),
    average_weight_(average_weight),
    mono_weight_(mono_weight)
  {
    this->setIsotopeDistribution(isotopes);
  }

  Element::~Element() = default;

  void Element::setAtomicNumber(unsigned int atomic_number)
  {
    atomic_number_ = atomic_number;
  }

  unsigned int Element::getAtomicNumber() const
  {
    return atomic_number_;
  }

  void Element::setAverageWeight(double weight)
  {
    average_weight_ = weight;
  }

  double Element::getAverageWeight() const
  {
    return average_weight_;
  }

  void Element::setMonoWeight(double weight)
  {
    mono_weight_ = weight;
  }

  double Element::getMonoWeight() const
  {
    return mono_weight_;
  }

  void Element::setIsotopeDistribution(const IsotopeDistribution & distribution)
  {
    //force sortedness by mz. A lot of code relies on this.
    assert(std::is_sorted(distribution.begin(), distribution.end(), Peak1D::MZLess()));
    isotopes_ = distribution;
  }

  const IsotopeDistribution & Element::getIsotopeDistribution() const
  {
    return isotopes_;
  }

  void Element::setName(const string & name)
  {
    name_ = name;
  }

  const string & Element::getName() const
  {
    return name_;
  }

  void Element::setSymbol(const string & symbol)
  {
    symbol_ = symbol;
  }

  const string & Element::getSymbol() const
  {
    return symbol_;
  }

  Element & Element::operator=(const Element & element) = default;

  bool Element::operator==(const Element & element) const
  {
    return name_ == element.name_ &&
           symbol_ == element.symbol_ &&
           atomic_number_ == element.atomic_number_ &&
           average_weight_ == element.average_weight_ &&
           mono_weight_ == element.mono_weight_ &&
           isotopes_ == element.isotopes_;
  }

  bool Element::operator<(const Element & rhs) const
  {
    return std::tie(
     atomic_number_, 
     mono_weight_, 
     symbol_, 
     name_, 
     average_weight_, 
     isotopes_) 
     < 
     std::tie(
      rhs.atomic_number_, 
      rhs.mono_weight_, 
      rhs.symbol_, 
      rhs.name_, 
      rhs.average_weight_, 
      rhs.isotopes_);
  }

  bool Element::operator!=(const Element & element) const
  {
    return !(*this == element);
  }

  std::ostream & operator<<(std::ostream & os, const Element & element)
  {
    os  << element.name_ << " "
    << element.symbol_ << " "
    << element.atomic_number_ << " "
    << element.average_weight_ << " "
    << element.mono_weight_;

    for (const auto& isotope : element.isotopes_)
    {
      if (isotope.getIntensity() > 0.0f)
      {
        os << " " << isotope.getPosition() << "=" << isotope.getIntensity() * 100 << "%";
      }
    }
    return os;
  }

} // namespace OpenMS
