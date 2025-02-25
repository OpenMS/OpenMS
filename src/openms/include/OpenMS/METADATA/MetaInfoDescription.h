// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

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
    String comment_;
    String name_;
    std::vector<DataProcessingPtr> data_processing_;
  };
} // namespace OpenMS

