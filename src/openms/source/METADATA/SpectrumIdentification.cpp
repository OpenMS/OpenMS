// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SpectrumIdentification.h>

using namespace std;

namespace OpenMS
{

  SpectrumIdentification::~SpectrumIdentification() = default;

  // Equality operator
  bool SpectrumIdentification::operator==(const SpectrumIdentification & rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && id_ == rhs.id_
           && hits_ == rhs.hits_;
  }

  // Inequality operator
  bool SpectrumIdentification::operator!=(const SpectrumIdentification & rhs) const
  {
    return !(*this == rhs);
  }

  void SpectrumIdentification::setHits(const vector<IdentificationHit> & hits)
  {
    hits_ = hits;
  }

  void SpectrumIdentification::addHit(const IdentificationHit & hit)
  {
    hits_.push_back(hit);
  }

  const vector<IdentificationHit> & SpectrumIdentification::getHits() const
  {
    return hits_;
  }

} // namespace OpenMS

