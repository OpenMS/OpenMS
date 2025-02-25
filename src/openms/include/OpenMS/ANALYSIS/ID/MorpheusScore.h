// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Macros.h>

namespace OpenMS
{

/**
 *  @brief An implementation of the Morpheus PSM scoring function
 *  Inspired by a C# implementation by C. Wenger released under MIT license
 */               
struct OPENMS_DLLAPI MorpheusScore
{

  /// score and subscores
  struct OPENMS_DLLAPI Result
  {
    Size matches = 0; ///< matched theoretical peaks
    Size n_peaks = 0; ///< number of theoretical peaks
    float score = 0; ///< Morpheus score (matched peaks + matched ion current / TIC)
    float MIC = 0; ///< ion current of matches (experimental peaks)
    float TIC = 0; ///< total ion current (experimental peak) 
    float err = 0; ///< average absolute mass error of matched fragments (in Da)
  };

  /// returns Morpheus Score, \#matched ions, \#total ions, \#matched intensities, \#total fragment intensities (TIC)
  static Result compute(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const PeakSpectrum& theo_spectrum);
};

}


