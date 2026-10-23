// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/METADATA/PeptideHit.h>

#include <memory>
#include <string>
#include <vector>

namespace OpenMS
{
  class MSExperiment;
  class MSSpectrum;
  class PeptideIdentificationList;

/**
  @brief The parts of a QPX dataset that a search engine adds to the views of QPXFile

  A tool that writes a QPX dataset (e.g. ProSE with @c -out_qpx) takes the mandatory columns of the psm view from
  QPXFile and adds the optional ones of this class with ArrowIOHelpers::withColumns(). It writes the mz view of its
  identified spectra with mzView(). Columns are ArrowIOHelpers::ParquetColumn and tables are only forward-declared,
  so the caller needs no Arrow headers.

  All members are static.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI QPXDatasetExport
{
public:
  /// A column of a QPX view, without Arrow types
  using Column = ArrowIOHelpers::ParquetColumn;

  /**
    @brief The name of a score or meta value in the snake_case of the QPX @c additional_scores

    Percolator's scores get their names of bigbio/qpx#355 (e.g. "MS:1001492" -> "percolator_score"), "q-value" is
    "global_qvalue", others are converted ("ln(hyperscore)" -> "ln_hyperscore", "DeltaMass" -> "delta_mass",
    "COMET:xcorr" -> "comet_xcorr").

    @param[in] name The score type or meta value key
    @return the QPX name; "" for the posterior error probability, which has a column of its own
  */
  static std::string scoreName(const std::string& name);

  /**
    @brief The optional columns of the QPX psm view, one entry per row of the view (the best hit of each element of
    @p peptide_ids with hits; elements without hits have no row)

    The columns, in this order:
      - @c additional_scores (SCORE_LIST): the main score and every numeric meta value of the hit, named by
        scoreName() (each name once, the posterior error probability and non-score meta values such as
        @c target_decoy left out), with @c higher_better for the scores whose direction is known;
      - @c spectrum_reference (UTF8): the native ID; null if the element has no spectrum (@p has_spectrum) or no
        native ID;
      - @c isotope_error (INT16): the meta value @c isotope_error (or @c IsotopeError) as observed minus theoretical
        (also so in @c additional_scores); null without one;
      - @c protein_mappings (STRUCT_LIST of @c accession, @c start, @c end, @c aa_before, @c aa_after): one entry per
        peptide evidence, 1-based positions (null if unknown), "-" for a protein terminus, null for an unknown residue.

    @param[in] peptide_ids The identifications of one run, in the order of the rows of the psm view
    @param[in] has_spectrum Per element of @p peptide_ids whether the mz view has its spectrum; empty: all true
    @param[in] isotope_observed_minus_theoretical Whether the isotope error meta values are observed minus
               theoretical (false: theoretical minus observed, their sign is flipped)
    @return the four columns
    @throws Exception::IllegalArgument if @p has_spectrum is neither empty nor of the size of @p peptide_ids
  */
  static std::vector<Column> psmColumns(const PeptideIdentificationList& peptide_ids,
                                        const std::vector<bool>& has_spectrum,
                                        bool isotope_observed_minus_theoretical);

  /**
    @brief The peak filter of the mz view: the @p top most intense peaks per @p window Th (@p top 0: all peaks), and
    the annotated ones (keptPeaks())
  */
  struct PeakFilter
  {
    /// The width of a window (Th)
    double window = 100.0;
    /// The most intense peaks kept per window; 0: all peaks
    Size top = 100;
    /// The fragment annotation is of deisotoped spectra (e.g. ProSE's 'Search:fragment:deisotope' with a
    /// high-resolution fragment tolerance): the deisotoper moved each fragment of charge z (up to @p max_charge, or the
    /// precursor's charge with @p cap_at_precursor) to its 1+ m/z, which is no peak of the spectrum as the file has it.
    bool deisotoped = false;
    /// The highest fragment charge of the deisotoper
    int max_charge = 3;
    /// The deisotoper capped the fragment charge at the precursor's charge
    bool cap_at_precursor = true;
  };

  /**
    @brief The indices of the peaks of @p spectrum that @p filter keeps: the most intense ones per window and those of
    the fragment annotation @p annotations; all peaks if @p filter.top is 0

    The windows are jumping, as WindowMower's (a window starts at its first peak, the next at the first peak past it),
    but the last window keeps the full quota too. An annotated m/z keeps the peak nearest to it (within 0.01 Th). Of
    deisotoped spectra (@p filter.deisotoped) it also keeps, for each charge z from 2 to the deisotoper's maximum
    (@p filter.max_charge, or @p precursor_charge if lower and @p filter.cap_at_precursor), the peak nearest to
    (m/z + (z - 1) * proton) / z (within 0.01 Th): the monoisotopic peak of a fragment of charge z, which the
    deisotoper moved to the annotated 1+ m/z. Every annotated fragment thus keeps its peak in the spectrum as the file
    has it; a 1+ fragment may keep an unrelated peak at its 2+ or 3+ position too.

    @param[in] spectrum The spectrum, sorted by m/z
    @param[in] annotations The fragment annotation of a PSM of @p spectrum
    @param[in] precursor_charge The spectrum's recorded precursor charge (as in the search's preprocessing); 0 if
               unknown (then up to @p filter.max_charge)
    @param[in] filter The peak filter
    @return the indices of the kept peaks, ascending
  */
  static std::vector<Size> keptPeaks(const MSSpectrum& spectrum, const std::vector<PeptideHit::PeakAnnotation>& annotations,
                                     int precursor_charge, const PeakFilter& filter);

  /// The QPX mz view of the spectra of the PSMs of one run (mzView())
  struct MzView
  {
    /// The mz view (null if it could not be built; the error is logged)
    std::shared_ptr<arrow::Table> table;
    /// Per element of the peptide identifications whether the view has the row of its spectrum (false for an element
    /// without hits, and for all if @p table is null)
    std::vector<bool> has_spectrum;
    /// The PSMs without a spectrum of their native ID (and RT), whose native ID has no row either
    Size without_spectrum = 0;
    /// The PSMs of a spectrum that shares its native ID with another spectrum (at another RT), which has the row
    Size unresolved = 0;
    /// The peaks of the spectra of the rows
    Size peaks_read = 0;
    /// The peaks of the rows (kept by the peak filter)
    Size peaks_kept = 0;
  };

  /**
    @brief The QPX mz view of the MS2 spectra of the PSMs of one run, with the peaks @p filter keeps

    The spectrum of a PSM (the best hit of an element of @p peptide_ids with hits) is the MS2 spectrum of @p spectra of
    its native ID; of several spectra of its native ID the one at its RT (within 1e-6 s). The view has one row per
    native ID (@c id is its key), in the order of the first PSM of each: the row of the first spectrum of a native ID
    that a PSM has. A PSM whose spectrum is not found has none, nor has one whose spectrum shares its native ID with
    the spectrum of an earlier row (unresolved). The peaks of a row are those that keptPeaks() keeps for any of the
    PSMs of its spectrum (each keeps its own annotated peaks), with the precursor charge of the spectrum's first
    precursor (not the charge of a PSM: a spectrum of unknown charge was deisotoped up to @p filter.max_charge, whatever
    the charge of its PSMs). A spectrum not sorted by m/z is sorted (in a copy).

    The columns, in this order: @c id, @c run_file_name, @c scan (ArrowIOHelpers::qpxScanComponents()), @c ms_level,
    @c centroid, @c scan_start_time (minutes), @c inverse_ion_mobility (null unless in Vs/cm^2), @c ion_injection_time
    (null without MS:1000927), @c total_ion_current, @c precursors, @c mz, @c intensity, @c cv_params (null).

    @param[in] spectra The spectra of the run (the MS2 spectra of its PSMs suffice; other spectra are left out)
    @param[in] peptide_ids The identifications of the run, in the order of the rows of the psm view
    @param[in] run_file_name The run's @c run_file_name (e.g. ArrowIOHelpers::qpxRunFileName() of its file)
    @param[in] filter The peak filter
    @return the view, @c has_spectrum with one entry per element of @p peptide_ids (for psmColumns()), and its counts
  */
  static MzView mzView(const MSExperiment& spectra, const PeptideIdentificationList& peptide_ids,
                       const std::string& run_file_name, const PeakFilter& filter);
};

} // namespace OpenMS
