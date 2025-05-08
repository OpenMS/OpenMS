// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Acquisition.h>

using namespace std;

namespace OpenMS
{

  bool Acquisition::operator==(const Acquisition & rhs) const
  {
    return identifier_ == rhs.identifier_ && MetaInfoInterface::operator==(rhs);
  }

  bool Acquisition::operator!=(const Acquisition & rhs) const
  {
    return !(operator==(rhs));
  }

  const String & Acquisition::getIdentifier() const
  {
    return identifier_;
  }

  void Acquisition::setIdentifier(const String & identifier)
  {
    identifier_ = identifier;
  }

}

