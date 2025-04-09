// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MetaInfoDescription.h>

#include <OpenMS/CONCEPT/Helpers.h>

using namespace std;

namespace OpenMS
{

  MetaInfoDescription::~MetaInfoDescription() = default;

  bool MetaInfoDescription::operator==(const MetaInfoDescription & rhs) const
  {
    return MetaInfoInterface::operator==(rhs) &&
           comment_ == rhs.comment_ &&
           name_ == rhs.name_ &&
           ( data_processing_.size() == rhs.data_processing_.size() &&
           std::equal(data_processing_.begin(),
                      data_processing_.end(),
                      rhs.data_processing_.begin(),
                      OpenMS::Helpers::cmpPtrSafe<DataProcessingPtr>) );
  }

  void MetaInfoDescription::setName(const String & name)
  {
    name_ = name;
  }

  const String & MetaInfoDescription::getName() const
  {
    return name_;
  }

  const vector<ConstDataProcessingPtr> & MetaInfoDescription::getDataProcessing() const
  {
    return OpenMS::Helpers::constifyPointerVector(data_processing_);
  }

  vector<DataProcessingPtr> & MetaInfoDescription::getDataProcessing()
  {
    return data_processing_;
  }

  void MetaInfoDescription::setDataProcessing(const vector<DataProcessingPtr> & processing_method)
  {
    data_processing_ = processing_method;
  }

}

