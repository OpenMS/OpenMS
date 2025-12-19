// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
      @param[in] custom_presets_file Optional path to a custom presets file
      @return StringList containing all available preset names
    */
    OPENMS_DLLAPI StringList getAllPresetsNames(const String& custom_presets_file = "");

    /**
      @brief Get preset parameters for a given preset name
      @param[in] p The preset name
      @param[in] custom_presets_file Optional path to a custom presets file
      @param[out] nucleotides Output parameter for nucleotides
      @param[out] mapping Output parameter for mapping
      @param[out] modifications Output parameter for modifications
      @param[out] fragment_adducts Output parameter for fragment adducts
      @param[out] can_cross_link Output parameter for can_cross_link
    */
   OPENMS_DLLAPI void getPresets(const String& p, 
    const String& custom_presets_file,
    StringList& nucleotides, 
    StringList& mapping, 
    StringList& modifications, 
    StringList& fragment_adducts, 
    String& can_cross_link); 

    /**
      @brief Get preset parameters for a given preset name (using default presets file)
      @param[in] p The preset name
      @param[out] nucleotides Output parameter for nucleotides
      @param[out] mapping Output parameter for mapping
      @param[out] modifications Output parameter for modifications
      @param[out] fragment_adducts Output parameter for fragment adducts
      @param[out] can_cross_link Output parameter for can_cross_link
    */
    OPENMS_DLLAPI void getPresets(const String& p, 
     StringList& nucleotides, 
     StringList& mapping, 
     StringList& modifications, 
     StringList& fragment_adducts, 
     String& can_cross_link);
  }

}
