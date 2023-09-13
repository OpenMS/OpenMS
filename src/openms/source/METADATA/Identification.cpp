// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Identification.h>

using namespace std;

namespace OpenMS
{

  Identification::~Identification() = default;

  // Equality operator
  bool Identification::operator==(const Identification & rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && id_ == rhs.id_
           && creation_date_ == rhs.creation_date_
           && spectrum_identifications_ == rhs.spectrum_identifications_;
  }

  // Inequality operator
  bool Identification::operator!=(const Identification & rhs) const
  {
    return !(*this == rhs);
  }

  void Identification::setCreationDate(const DateTime & date)
  {
    creation_date_ = date;
  }

  const DateTime & Identification::getCreationDate() const
  {
    return creation_date_;
  }

  void Identification::setSpectrumIdentifications(const vector<SpectrumIdentification> & ids)
  {
    spectrum_identifications_ = ids;
  }

  void Identification::addSpectrumIdentification(const SpectrumIdentification & id)
  {
    spectrum_identifications_.push_back(id);
  }

  const vector<SpectrumIdentification> & Identification::getSpectrumIdentifications() const
  {
    return spectrum_identifications_;
  }

} // namespace OpenMS

