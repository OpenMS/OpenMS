// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sachsenberg $
// $Authors: Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzymeDataProvider.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/config.h>

namespace OpenMS
{
  /**
    @ingroup Chemistry

    @brief Data provider for built-in protease definitions

    Provides the default set of proteases (Trypsin, Arg-C, Lys-C, etc.)
    as hardcoded definitions, eliminating the need for an XML file at runtime.
  */
  class OPENMS_DLLAPI BuiltInProteaseDataProvider : public DigestionEnzymeDataProvider<DigestionEnzymeProtein>
  {
  public:
    /// @brief Returns the hardcoded set of built-in protease definitions.
    std::vector<std::unique_ptr<DigestionEnzymeProtein>> loadEnzymes() override;
  };

} // namespace OpenMS
