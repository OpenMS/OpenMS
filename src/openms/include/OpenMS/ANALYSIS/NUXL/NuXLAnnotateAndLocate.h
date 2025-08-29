// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLParameterParsing.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLAnnotatedHit.h>

#include <set>
#include <map>
#include <vector>
#include <algorithm>

namespace OpenMS
{
class OPENMS_DLLAPI NuXLAnnotateAndLocate
{
  public:

  static void annotateAndLocate_(
    const PeakMap& exp, 
    std::vector<std::vector<NuXLAnnotatedHit>>& annotated_hits,
    const NuXLModificationMassesResult& mm,
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
    Size max_variable_mods_per_peptide, 
    double fragment_mass_tolerance, 
    bool fragment_mass_tolerance_unit_ppm, 
    const NuXLParameterParsing::PrecursorsToMS2Adducts & all_feasible_adducts);

};
} // namespace OpenMS


