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
           name_ == rhs.name_ &&
           ( data_processing_.size() == rhs.data_processing_.size() &&
           std::equal(data_processing_.begin(),
                      data_processing_.end(),
                      rhs.data_processing_.begin(),
                      OpenMS::Helpers::cmpPtrSafe<DataProcessingPtr>) );
  }

  auto MetaInfoDescription::operator<=>(const MetaInfoDescription & rhs) const
  {
    // First compare the MetaInfoInterface base
    if (MetaInfoInterface::operator!=(rhs))
    {
      // For MetaInfoInterface comparison, we need to create a deterministic ordering
      // Compare by checking if one is empty and the other is not
      bool lhs_empty = isMetaEmpty();
      bool rhs_empty = rhs.isMetaEmpty();
      if (lhs_empty != rhs_empty)
      {
        return lhs_empty <=> rhs_empty;
      }
      
      // If both non-empty, compare by getting keys and comparing them
      std::vector<UInt> lhs_keys, rhs_keys;
      getKeys(lhs_keys);
      rhs.getKeys(rhs_keys);
      
      if (auto cmp = lhs_keys <=> rhs_keys; cmp != 0)
        return cmp;
      
      // If keys are same, compare values for each key
      for (UInt key : lhs_keys)
      {
        if (auto cmp = getMetaValue(key).toString() <=> rhs.getMetaValue(key).toString(); cmp != 0)
          return cmp;
      }
    }
    
    // Compare name
    if (auto cmp = name_ <=> rhs.name_; cmp != 0)
      return cmp;
    
    // Compare data processing by size first
    if (auto cmp = data_processing_.size() <=> rhs.data_processing_.size(); cmp != 0)
      return cmp;
    
    // Compare data processing elements (simplified comparison by comparing names)
    for (size_t i = 0; i < data_processing_.size(); ++i)
    {
      if (data_processing_[i] && rhs.data_processing_[i])
      {
        if (auto cmp = data_processing_[i]->getProcessingActions().size() <=> rhs.data_processing_[i]->getProcessingActions().size(); cmp != 0)
          return cmp;
      }
      else if (data_processing_[i] != rhs.data_processing_[i])
      {
        return (data_processing_[i] != nullptr) <=> (rhs.data_processing_[i] != nullptr);
      }
    }
    
    return std::strong_ordering::equal;
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

