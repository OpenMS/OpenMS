// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <map>

namespace OpenMS
{

class MSExperiment;

/** @brief Preprocessing of nucleotide cross-link mass spectra.

  This class provides methods for preprocessing mass spectra used in NuXL analysis.
  It includes deisotoping, charge annotation, peak filtering, normalization, and
  interference peak removal.
  **/
class OPENMS_DLLAPI NuXLPreprocessSpectra
{
  public:

  /** @brief Preprocess spectra for NuXL analysis.

    This method performs comprehensive preprocessing of MS/MS spectra including:
    - Removal of zero intensities
    - Deisotoping and single charge conversion
    - Peak interference filtering based on precursor purity
    - Normalization
    - Window-based and top-N peak filtering
    - TIC calculation and annotation

    @param exp Input/output MS experiment containing spectra to preprocess
    @param single_charge_spectra Whether to convert all peaks to single charge
    @param annotate_charge Whether to annotate charge states in IntegerDataArray
    @param window_size Size of sliding window for WindowMower filtering
    @param peakcount Number of peaks to keep in WindowMower filtering
    @param purities Map of spectrum native IDs to precursor purity scores for interference filtering
    */
    static void preprocessSpectra(MSExperiment& exp, 
                                  bool single_charge_spectra, 
                                  bool annotate_charge,
                                  double window_size,
                                  Size peakcount,
                                  const std::map<String, PrecursorPurity::PurityScores>& purities);

  private:

  /** @brief Filter peaks that interfere with precursor isolation.

    Removes peaks that match interfering peaks identified during precursor isolation,
    based on the precursor purity analysis.

    @param spectra Input/output MS experiment containing spectra to filter
    @param purities Map of spectrum native IDs to precursor purity scores
    @param fragment_mass_tolerance Mass tolerance for peak matching (default: 20.0)
    @param fragment_mass_tolerance_unit_ppm Whether tolerance is in ppm (default: true)
    */
    static void filterPeakInterference_(MSExperiment& spectra, 
                                        const std::map<String, PrecursorPurity::PurityScores>& purities, 
                                        double fragment_mass_tolerance = 20.0, 
                                        bool fragment_mass_tolerance_unit_ppm = true);
};

}