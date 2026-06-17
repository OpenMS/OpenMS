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
    IM_CCS,   ///< collisional cross section (square angstrom)
    SIZE_OF_DIM_UNITS
  };
  inline std::string_view DIM_NAMES[(int)DIM_UNIT::SIZE_OF_DIM_UNITS] = {"RT [s]", "m/z [Th]", "intensity", "IM [milliseconds]", "IM [vs / cm2]", "FAIMS CV", "CCS [angstrom^2]"};
  inline std::string_view DIM_NAMES_SHORT[(int)DIM_UNIT::SIZE_OF_DIM_UNITS] = {"RT", "m/z", "int", "IM", "IM", "FCV", "CCS"};


  enum class MZ_UNITS
  {
    DA = 0, ///< Dalton
    PPM,    ///< parts-per-million
    SIZE_OF_MZ_UNITS
  };
  inline std::string_view MZ_UNIT_NAMES[(int)MZ_UNITS::SIZE_OF_MZ_UNITS] = {"Da", "ppm"};

  /**
  @brief Enum for different hydrophobicity scales

  If a new scale is introduced, append it to the list below and assign it an enum value equal to the current maximum enum value plus one.

  Add the data for this scale here: @ref Residue::getHydrophobicity
  */
  enum class HydrophobicityScaleMethod 
    {
      KYTE_DOOLITTLE = 0,
      EISENBERG = 1,
      HOPP_WOODS = 2,
      BULL_BREESE = 3,
      BLACK_MOULD = 4,
      GUY = 5,
      EISENBERG_CONSENSUS = 6
    }; 

} // namespace OpenMS
