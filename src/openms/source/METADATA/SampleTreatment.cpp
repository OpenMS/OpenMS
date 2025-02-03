// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SampleTreatment.h>

using namespace std;

namespace OpenMS
{

  SampleTreatment::SampleTreatment(const String & type) :
    MetaInfoInterface(),
    type_(type),
    comment_()
  {
  }

  SampleTreatment::~SampleTreatment() = default;

  const String & SampleTreatment::getType() const
  {
    return type_;
  }

  const String & SampleTreatment::getComment() const
  {
    return comment_;
  }

  void SampleTreatment::setComment(const String & comment)
  {
    comment_ = comment;
  }

  bool SampleTreatment::operator==(const SampleTreatment & rhs) const
  {
    return MetaInfoInterface::operator==(rhs) &&
           comment_ == rhs.comment_;
  }

}
