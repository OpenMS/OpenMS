// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <string_view>

namespace OpenMS
{
  // add common enums here to avoid big includes of large classes and break circular dependencies

  /// Enum for different units which can be displayed on a plotting axis
  /// The order is arbitrary.
  enum class DIM_UNIT
  {
    RT = 0,   ///< RT in seconds
    MZ,       ///< m/z
    INT,      ///< intensity
    IM_MS,    ///< ion mobility milliseconds
    IM_VSSC,  ///< volt-second per square centimeter (i.e. 1/K_0)
    FAIMS_CV, ///< FAIMS compensation voltage
    SIZE_OF_DIM_UNITS
  };
  inline std::string_view DIM_NAMES[(int)DIM_UNIT::SIZE_OF_DIM_UNITS] = {"RT [s]", "m/z [Th]", "intensity", "IM [milliseconds]", "IM [vs / cm2]", "FAIMS CV"};
  inline std::string_view DIM_NAMES_SHORT[(int)DIM_UNIT::SIZE_OF_DIM_UNITS] = {"RT", "m/z", "int", "IM", "IM", "FCV"};


  enum class MZ_UNITS
  {
    DA = 0, ///< Dalton
    PPM,    ///< parts-per-million
    SIZE_OF_MZ_UNITS
  };
  inline std::string_view MZ_UNIT_NAMES[(int)MZ_UNITS::SIZE_OF_MZ_UNITS] = {"Da", "ppm"};

} // namespace OpenMS
