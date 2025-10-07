// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/CVReference.h>

using namespace std;

namespace OpenMS
{
  // CV reference implementation
  CVReference::CVReference() = default;

  CVReference::~CVReference() = default;

  CVReference::CVReference(const CVReference& rhs) = default;

  CVReference& CVReference::operator=(const CVReference& rhs)
  {
    if (this != &rhs)
    {
      name_ = rhs.name_;
      identifier_ = rhs.identifier_;
    }
    return *this;
  }

  bool CVReference::operator==(const CVReference& rhs) const
  {
    return name_ == rhs.name_ && identifier_ == rhs.identifier_;
  }

  bool CVReference::operator!=(const CVReference& rhs) const
  {
    return !(*this == rhs);
  }

  void CVReference::setName(const String& name)
  {
    name_ = name;
  }

  const String& CVReference::getName() const
  {
    return name_;
  }

  void CVReference::setIdentifier(const String& identifier)
  {
    identifier_ = identifier;
  }

  const String& CVReference::getIdentifier() const
  {
    return identifier_;
  }

} // namespace OpenMS
