// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <string>
#include <vector>

namespace OpenMS
{
  class MSSpectrum;

  /**
    @brief Finds precursors co-isolated with the selected one in the isolation window of an MS2 spectrum

    A chimeric MS2 spectrum holds the fragments of more than one precursor: those that fell into the isolation
    window next to the selected one. detect() looks for them in the survey (MS1) scan: each peak in the window, most
    intense first, is taken as a monoisotopic peak at each charge, and accepted if the isotope envelope found there
    fits an averagine envelope (Kullback-Leibler divergence of the observed from the averagine envelope at most
    @p max_kl). The selected precursor and its isotope peaks, and the isotope peaks of an envelope already accepted,
    are skipped. Since above about 1.8 kDa the M+1 peak is the most intense, an accepted envelope is moved down by
    one isotope step (up to three times) while there is a peak and the envelope from it fits averagine better.

    Ported from bigbio/andes (chimeric cascade, #118, #136).

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI CoIsolationDetector
  {
  public:
    /// a co-isolated precursor
    struct OPENMS_DLLAPI CoIsolated
    {
      double mono_mz = 0.0; ///< m/z of the monoisotopic peak
      int charge = 0;       ///< charge
      double kl = 0.0;      ///< divergence of its isotope envelope from averagine
    };

    /// meta value of an MS2 spectrum: the monoisotopic m/z of its co-isolated precursors (DoubleList, see annotate())
    static const std::string META_MZ;
    /// meta value of an MS2 spectrum: the charges of its co-isolated precursors (IntList, same order as META_MZ)
    static const std::string META_CHARGE;

    /**
      @brief Divergence of the isotope envelope at @p mono_mz from averagine

      The observed envelope is the most intense peak within @p tol_da of each of the first @p n_isotopes isotope
      positions (0 if none), scaled to sum 1. Returns sum_k o_k ln(o_k / a_k) over the isotopes with o_k > 0, with
      a the averagine envelope of the same mass (scaled to sum 1); 10 if no isotope position has a peak.

      @param ms1 Survey scan, sorted by m/z
    */
    static double isotopeKL(const MSSpectrum& ms1, double mono_mz, int charge, double tol_da, Size n_isotopes = 4);

    /**
      @brief Co-isolated precursors in the isolation window [@p window_lower, @p window_upper] of @p ms1

      @param ms1 Survey scan, sorted by m/z
      @param selected_mz m/z of the selected precursor (it and its isotope peaks are skipped)
      @param min_charge Lowest charge tried
      @param max_charge Highest charge tried
      @param tol_da m/z tolerance (Th) for the isotope peaks
      @param max_kl Largest accepted divergence from averagine
      @param max_n At most this many precursors (most intense monoisotopic peaks first); 0 finds none
    */
    static std::vector<CoIsolated> detect(const MSSpectrum& ms1, double window_lower, double window_upper,
                                          double selected_mz, int min_charge, int max_charge, double tol_da,
                                          double max_kl, Size max_n);

    /**
      @brief Annotate each MS2 spectrum of @p exp with the precursors co-isolated with its first precursor

      The survey scan is MSExperiment::getPrecursorSpectrum(); the isolation window is centred on the
      "isolation window target m/z" of the precursor (its m/z if absent) and spans its isolation window offsets.
      Spectra without a survey scan or an isolation window are left out. Writes the meta values META_MZ and
      META_CHARGE (also when nothing was found, as empty lists). Sorts the survey scans by m/z if they are not.

      @param tol_ppm m/z tolerance for the isotope peaks, in ppm of the selected precursor's m/z
      @return Number of co-isolated precursors found
    */
    static Size annotate(PeakMap& exp, int min_charge, int max_charge, double tol_ppm, double max_kl, Size max_n);
  };
} // namespace OpenMS
