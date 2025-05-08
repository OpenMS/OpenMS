// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Software.h>

using namespace std;

namespace OpenMS
{

  Software::Software(const String& name, const String& version) :
    CVTermList(),
    name_(name),
    version_(version)
  {
  }

  Software::~Software() = default;

  bool Software::operator==(const Software& rhs) const
  {
    return CVTermList::operator==(rhs) &&
           name_ == rhs.name_ &&
           version_ == rhs.version_;
  }

  bool Software::operator!=(const Software& rhs) const
  {
    return !(operator==(rhs));
  }

  bool Software::operator<(const Software& rhs) const
  {
    return tie(name_, version_) < tie(rhs.name_, rhs.version_);
  }

  const String& Software::getName() const
  {
    return name_;
  }

  void Software::setName(const String& name)
  {
    name_ = name;
  }

  const String& Software::getVersion() const
  {
    return version_;
  }

  void Software::setVersion(const String& version)
  {
    version_ = version;
  }

}
