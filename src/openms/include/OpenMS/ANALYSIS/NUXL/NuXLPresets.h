// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  // MS1 (modifications) and MS2 (fragments) for the different protocols
  namespace NuXLPresets
  {
    /**
      @brief Get all available presets names from JSON file
      @return StringList containing all available preset names
    */
    OPENMS_DLLAPI StringList getAllPresetsNames();

    /**
      @brief Get preset parameters for a given preset name
      @param p The preset name
      @param nucleotides Output parameter for nucleotides
      @param mapping Output parameter for mapping
      @param modifications Output parameter for modifications
      @param fragment_adducts Output parameter for fragment adducts
      @param can_cross_link Output parameter for can_cross_link
    */
   OPENMS_DLLAPI void getPresets(const String& p, 
    StringList& nucleotides, 
    StringList& mapping, 
    StringList& modifications, 
    StringList& fragment_adducts, 
    String& can_cross_link); 
  }

}
