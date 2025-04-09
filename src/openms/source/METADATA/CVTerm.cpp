// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/CVTerm.h>

using namespace std;

namespace OpenMS
{

  CVTerm::CVTerm(const String& accession, const String& name, const String& cv_identifier_ref, const String& value, const Unit& unit) :
    accession_(accession),
    name_(name),
    cv_identifier_ref_(cv_identifier_ref),
    unit_(unit),
    value_(value)
  {
  }

  CVTerm::~CVTerm() = default;

  bool CVTerm::operator==(const CVTerm & rhs) const
  {
    return accession_ == rhs.accession_ &&
           name_ == rhs.name_ &&
           cv_identifier_ref_ == rhs.cv_identifier_ref_ &&
           unit_ == rhs.unit_ &&
           value_ == rhs.value_;
  }

  bool CVTerm::operator!=(const CVTerm& rhs) const
  {
    return !(*this == rhs);
  }

  void CVTerm::setAccession(const String& accession)
  {
    accession_ = accession;
  }

  const String& CVTerm::getAccession() const
  {
    return accession_;
  }

  void CVTerm::setName(const String& name)
  {
    name_ = name;
  }

  const String& CVTerm::getName() const
  {
    return name_;
  }

  void CVTerm::setCVIdentifierRef(const String& cv_identifier_ref)
  {
    cv_identifier_ref_ = cv_identifier_ref;
  }

  const String& CVTerm::getCVIdentifierRef() const
  {
    return cv_identifier_ref_;
  }

  void CVTerm::setUnit(const Unit& unit)
  {
    unit_ = unit;
  }

  const CVTerm::Unit& CVTerm::getUnit() const
  {
    return unit_;
  }

  void CVTerm::setValue(const DataValue& value)
  {
    value_ = value;
  }

  const DataValue& CVTerm::getValue() const
  {
    return value_;
  }

  bool CVTerm::hasUnit() const
  {
    return !unit_.accession.empty();
  }

  bool CVTerm::hasValue() const
  {
    return value_ != DataValue::EMPTY;
  }

} // namespace OpenMS

