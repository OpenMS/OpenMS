// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <compare>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
  /**
      @brief Description of the meta data arrays of MSSpectrum.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI MetaInfoDescription :
    public MetaInfoInterface
  {
public:
    /// Constructor
    MetaInfoDescription() = default;
    /// Copy constructor
    MetaInfoDescription(const MetaInfoDescription &) = default;
    /// Move constructor
    MetaInfoDescription(MetaInfoDescription&&) = default;
    /// Destructor
    ~MetaInfoDescription();

    /// Assignment operator
    MetaInfoDescription & operator=(const MetaInfoDescription &) = default;
    /// Move assignment operator
    MetaInfoDescription& operator=(MetaInfoDescription&&) & = default;

    /// Equality operator
    bool operator==(const MetaInfoDescription & rhs) const;
    /// Three-way comparison operator
    auto operator<=>(const MetaInfoDescription & rhs) const
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

    /// returns the name of the peak annotations
    const String & getName() const;
    /// sets the name of the peak annotations
    void setName(const String & name);

    /// returns a const reference to the description of the applied processing
    const std::vector<ConstDataProcessingPtr> & getDataProcessing() const;
    /// returns a mutable reference to the description of the applied processing
    std::vector<DataProcessingPtr> & getDataProcessing();
    /// sets the description of the applied processing
    void setDataProcessing(const std::vector<DataProcessingPtr> & data_processing);

protected:
    String name_;
    std::vector<DataProcessingPtr> data_processing_;
  };
} // namespace OpenMS

