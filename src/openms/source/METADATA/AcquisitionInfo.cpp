// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/AcquisitionInfo.h>

using namespace std;

namespace OpenMS
{

  bool AcquisitionInfo::operator==(const AcquisitionInfo & rhs) const
  {
    return method_of_combination_ == rhs.method_of_combination_ &&
           MetaInfoInterface::operator==(rhs) &&
           std::operator==(*this, rhs);
  }

  bool AcquisitionInfo::operator!=(const AcquisitionInfo & rhs) const
  {
    return !(operator==(rhs));
  }

  const String & AcquisitionInfo::getMethodOfCombination() const
  {
    return method_of_combination_;
  }

  void AcquisitionInfo::setMethodOfCombination(const String & method_of_combination)
  {
    method_of_combination_ = method_of_combination;
  }

}

