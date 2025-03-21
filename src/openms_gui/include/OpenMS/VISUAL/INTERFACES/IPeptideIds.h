// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h> 

#include <vector>

namespace OpenMS
{

  /**
  @brief Abstract base class which defines an interface for PeptideIdentifications
  */
  class OPENMS_GUI_DLLAPI IPeptideIds
  {
  public:
    using PepIds = std::vector<PeptideIdentification>;
    
    /// get the peptide IDs for this layer
    virtual const PepIds& getPeptideIds() const = 0;
    virtual PepIds& getPeptideIds() = 0;

    /// overwrite the peptide IDs for this layer
    virtual void setPeptideIds(const PepIds& ids) = 0;
    virtual void setPeptideIds(PepIds&& ids) = 0;
  };

}// namespace OpenMS
