// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <vector>

namespace OpenMS
{

  /// Strict-weak ordering on @ref MSSpectrum by the first precursor's m/z (ascending). Convenient as a sort predicate when building a precursor-m/z-indexed spectral library.
  struct OPENMS_DLLAPI PrecursorMassComparator
  {
    bool operator() (const MSSpectrum& a, const MSSpectrum& b)
    {
      return a.getPrecursors()[0].getMZ() < b.getPrecursors()[0].getMZ();
    }

  } PrecursorMZLess;

  /**
    @brief One library-search match record: the observed spectrum, the matching DB spectrum, the score, and the metabolite metadata copied from the hit.

    Used by @ref MetaboliteSpectralMatching to accumulate all candidate matches before the
    top-N export step. The fields fall into three groups:
      - Observed-side: the index / native id of the experimental MS2 spectrum and its
        precursor mass / RT.
      - DB-side: the index of the matching DB spectrum, the DB precursor mass / charge,
        and the hyperscore (see @ref MetaboliteSpectralMatching::computeHyperScore).
      - Metabolite metadata copied from the DB spectrum: primary/secondary identifier,
        common name, sum formula, InChI, SMILES, and the precursor adduct.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI SpectralMatch
  {
  public:
    /// Default constructor (zero-initialised numeric fields, empty strings)
    SpectralMatch();

    /// Default destructor
    ~SpectralMatch();

    /// Copy constructor
    SpectralMatch(const SpectralMatch&);

    /// Assignment operator
    SpectralMatch& operator=(const SpectralMatch&);

    /// Precursor mass (Da) of the experimental MS2 spectrum being matched
    double getObservedPrecursorMass() const;
    /// Setter for the observed-side precursor mass (Da)
    void setObservedPrecursorMass(const double&);

    /// Retention time (seconds) at which the observed MS2 spectrum was acquired
    double getObservedPrecursorRT() const;
    /// Setter for the observed-side precursor RT (seconds)
    void setObservedPrecursorRT(const double&);

    /// Precursor mass (Da) recorded on the DB-side spectrum that matched
    double getFoundPrecursorMass() const;
    /// Setter for the DB-side precursor mass (Da)
    void setFoundPrecursorMass(const double&);

    /// Precursor charge recorded on the DB-side spectrum
    Int getFoundPrecursorCharge() const;
    /// Setter for the DB-side precursor charge
    void setFoundPrecursorCharge(const Int&);

    /// Hyperscore (higher is better; 0 if fewer than three peaks matched — see @ref MetaboliteSpectralMatching::computeHyperScore)
    double getMatchingScore() const;
    /// Setter for the hyperscore
    void setMatchingScore(const double&);

    /// Index of the observed MS2 spectrum in the input @ref PeakMap passed to @ref MetaboliteSpectralMatching::run
    Size getObservedSpectrumIndex() const;
    /// Setter for the observed-spectrum index
    void setObservedSpectrumIndex(const Size&);

    /// Index of the matching DB spectrum in the spectral-library @ref PeakMap
    Size getMatchingSpectrumIndex() const;
    /// Setter for the DB-spectrum index
    void setMatchingSpectrumIndex(const Size&);

    /// Native id (vendor-specific id) of the observed MS2 spectrum
    String getObservedSpectrumNativeID() const;
    /// Setter for the observed-spectrum native id
    void setObservedSpectrumNativeID(const String&);

    /// Primary identifier of the matched metabolite. The .cpp populates this via the first available of: @c GNPS_Spectrum_ID, @c Massbank_Accession_ID, the metabolite-name meta value, or — as a last fallback — the DB spectrum's native id.
    String getPrimaryIdentifier() const;
    /// Setter for the primary metabolite identifier
    void setPrimaryIdentifier(const String&);

    /// Secondary identifier of the matched metabolite. The .cpp populates this from the @c HMDB_ID meta value on the DB-side spectrum.
    String getSecondaryIdentifier() const;
    /// Setter for the secondary metabolite identifier
    void setSecondaryIdentifier(const String&);

    /// Common (trivial) name of the matched metabolite
    String getCommonName() const;
    /// Setter for the common name
    void setCommonName(const String&);

    /// Chemical sum formula of the matched metabolite (neutral)
    String getSumFormula() const;
    /// Setter for the sum formula
    void setSumFormula(const String&);

    /// InChI string of the matched metabolite
    String getInchiString() const;
    /// Setter for the InChI string
    void setInchiString(const String&);

    /// SMILES string of the matched metabolite
    String getSMILESString() const;
    /// Setter for the SMILES string
    void setSMILESString(const String&);

    /// Adduct annotation recorded on the DB-side precursor (e.g. "[M+H]+")
    String getPrecursorAdduct() const;
    /// Setter for the precursor adduct
    void setPrecursorAdduct(const String&);

    /// Observed collision cross section (CCS) of the experimental precursor in Angstrom^2; -1.0 if not available
    double getObservedCCS() const;
    /// Setter for the observed CCS (Angstrom^2); use -1.0 to indicate "not available"
    void setObservedCCS(const double&);

    /// CCS (Angstrom^2) of the matching DB-side spectrum; -1.0 if not available
    double getFoundCCS() const;
    /// Setter for the DB-side CCS (Angstrom^2); use -1.0 to indicate "not available"
    void setFoundCCS(const double&);


  private:
    double observed_precursor_mass_;   ///< Observed-side precursor mass (Da)
    double observed_precursor_rt_;     ///< Observed-side precursor RT (seconds)
    double found_precursor_mass_;      ///< DB-side precursor mass (Da)
    Int found_precursor_charge_;       ///< DB-side precursor charge
    double matching_score_;            ///< Hyperscore (see @ref MetaboliteSpectralMatching::computeHyperScore)
    Size observed_spectrum_idx_;       ///< Index of the observed spectrum in the input PeakMap
    Size matching_spectrum_idx_;       ///< Index of the DB spectrum in the spectral-library PeakMap
    String observed_spectrum_native_id_; ///< Native id (vendor-specific) of the observed spectrum

    // further meta information
    String primary_id_;     ///< Primary metabolite identifier
    String secondary_id_;   ///< Secondary metabolite identifier
    String common_name_;    ///< Trivial / common name
    String sum_formula_;    ///< Neutral chemical sum formula
    String inchi_string_;   ///< InChI string
    String smiles_string_;  ///< SMILES string
    String precursor_adduct_; ///< DB-side precursor adduct annotation

    double observed_ccs_;   ///< Observed-side collision cross section (Angstrom^2); -1.0 if not available
    double found_ccs_;      ///< DB-side collision cross section (Angstrom^2); -1.0 if not available

  };

  /// Strict-weak ordering on @ref SpectralMatch by hyperscore (descending — top-scoring match comes first when sorted)
  struct OPENMS_DLLAPI SpectralMatchScoreComparator
  {
    bool operator() (const SpectralMatch& a, const SpectralMatch& b)
    {
      return a.getMatchingScore() > b.getMatchingScore();
    }

  } SpectralMatchScoreGreater;

  /**
    @brief Metabolomics spectral-library search by hyperscore.

    Matches the MS2 spectra in an experimental @ref PeakMap against the MS2 spectra of a
    spectral library (the .cpp recognises GNPS- and MassBank-style metadata on the DB-side
    spectrum — see @ref SpectralMatch::getPrimaryIdentifier for the actual lookup chain)
    and emits the matching results as an @ref MzTab table.

    Pipeline:
      -# Sort the library by precursor m/z (using @ref PrecursorMassComparator) for fast lookup.
      -# Apply @ref WindowMower noise reduction to each experimental MS2 spectrum.
      -# For each experimental MS2 spectrum, select DB spectra whose precursor is within the
         configured precursor m/z tolerance and ion-mode filter (parameter
         @c "ionization_mode", values @c "positive" or @c "negative").
      -# Score each candidate via @ref computeHyperScore — the hyperscore formula popularised
         by X!Tandem, @c log(dot_product) + @c log(factorial(matched_ions)). Matches with
         fewer than three matched ions are scored 0.
      -# Sort matches per spectrum by score (using the @ref SpectralMatchScoreComparator
         instance @c SpectralMatchScoreGreater) and emit either the top three or only the
         single best match per spectrum according to parameter @c "report_mode" (only
         @c "top3" or @c "best" are valid values) to MzTab.

    Optionally, candidates can additionally be filtered by collision cross section (CCS / ion
    mobility): when parameter @c "ccs_error_percent" is greater than 0, a candidate is discarded
    if the relative CCS error @c |observed - library| / library * 100 exceeds the configured
    percentage. The experimental CCS is taken (in this order) from a CCS-unit drift time, a
    1/K0 (VSSC) drift time converted via @ref IMTypes::oneOverK0ToCCS, or a numeric @c "CCS"
    meta value; the library CCS is taken from the DB-side spectrum the same way. If either side
    has no CCS, the candidate is kept (CCS does not penalise data lacking ion mobility). The
    observed/library CCS and the relative CCS error are always exported to MzTab when available.

    Configuration is exposed through @ref DefaultParamHandler — see the defaults installed
    by the constructor for the supported keys (@c "prec_mass_error_value",
    @c "frag_mass_error_value", @c "mass_error_unit", @c "ionization_mode", @c "report_mode",
    @c "merge_spectra", @c "ccs_error_percent").

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI MetaboliteSpectralMatching :
  public DefaultParamHandler,
  public ProgressLogger
  {
  public:
    /// Default constructor; installs the search parameters (see class docs)
    MetaboliteSpectralMatching();

    /// Default destructor
    ~MetaboliteSpectralMatching() override;

    /**
      @brief Compute the X!Tandem-style hyperscore for one experimental vs. DB spectrum pair.

      Formula: @c log(dot_product) + @c log(factorial(matched_ions)). The dot product
      sums @c (db_intensity * exp_intensity) over peak pairs where the experimental peak
      is within @p fragment_mass_error of a DB peak (ppm or Th per @p fragment_mass_tolerance_unit_ppm).
      Matches that find fewer than three matched ions return 0; matched-ion counts above
      @c boost::math::max_factorial saturate at the maximum factorial. Negative hyperscores
      are clamped to 0.

      @param[in] fragment_mass_error           Fragment-ion tolerance (value).
      @param[in] fragment_mass_tolerance_unit_ppm If @c true, @p fragment_mass_error is in ppm; otherwise Th.
      @param[in] exp_spectrum                  Experimental MS2 spectrum (sorted by m/z).
      @param[in] db_spectrum                   DB MS2 spectrum (sorted by m/z).
      @param[in] mz_lower_bound                Minimum m/z to consider; peaks below this are ignored. Default 0.0.
      @return Hyperscore (>= 0; 0 if fewer than three matched peaks).
    */
    static double computeHyperScore(
      double fragment_mass_error,
      bool fragment_mass_tolerance_unit_ppm,
      const MSSpectrum& exp_spectrum,
      const MSSpectrum& db_spectrum,
      double mz_lower_bound = 0.0);

    /**
      @brief Same as the other overload, but also fills @p annotations with per-peak match metadata.

      Annotations are produced only when the DB spectrum carries a @ref MSSpectrum::StringDataArray
      (annotation string) and a @ref MSSpectrum::IntegerDataArray (charge) — otherwise this
      overload behaves like the no-annotation form.

      @param[in]  fragment_mass_error           Fragment-ion tolerance (value).
      @param[in]  fragment_mass_tolerance_unit_ppm If @c true, @p fragment_mass_error is in ppm; otherwise Th.
      @param[in]  exp_spectrum                  Experimental MS2 spectrum (sorted by m/z).
      @param[in]  db_spectrum                   DB MS2 spectrum (sorted by m/z).
      @param[out] annotations                   Per-matched-peak annotations (annotation string, charge, m/z, intensity); appended to.
      @param[in]  mz_lower_bound                Minimum m/z to consider; peaks below this are ignored.
      @return Hyperscore (>= 0; 0 if fewer than three matched peaks).
    */
    static double computeHyperScore(
      double fragment_mass_error,
      bool fragment_mass_tolerance_unit_ppm,
      const MSSpectrum& exp_spectrum,
      const MSSpectrum& db_spectrum,
      std::vector<PeptideHit::PeakAnnotation>& annotations,
      double mz_lower_bound = 0.0);

    /**
      @brief Run the metabolomics spectral-library search and populate the MzTab output.

      Implements the pipeline described in the class brief: sort the DB by precursor m/z,
      noise-reduce each experimental MS2 via @ref WindowMower, score candidate library
      spectra, keep the top-N per query, and export to MzTab.

      @param[in,out] msexp       Experimental MS2 spectra. Modified in place by the noise-reduction step.
      @param[in,out] spec_db     Spectral library; will be sorted by precursor m/z.
      @param[out]    mztab_out   MzTab result table populated with the top matches.
      @param[in]     out_spectra If non-empty, the (noise-reduced) experimental spectra are also written to this mzML path.
    */
    void run(PeakMap & msexp, PeakMap & spec_db, MzTab & mztab_out, String & out_spectra);

  protected:
    void updateMembers_() override;

    /// Shared implementation behind the two public computeHyperScore() overloads; uses a pointer for @p annotations because mutable references cannot carry a temporary default value.
    static double computeHyperScore_(
      double fragment_mass_error,
      bool fragment_mass_tolerance_unit_ppm,
      const MSSpectrum& exp_spectrum,
      const MSSpectrum& db_spectrum,
      std::vector<PeptideHit::PeakAnnotation>* annotations = 0,
      double mz_lower_bound = 0.0);

  private:
    /// Convert the accumulated @ref SpectralMatch records into an MzTab table.
    void exportMzTab_(const std::vector<SpectralMatch>&, MzTab&);

    double precursor_mz_error_;     ///< Precursor m/z tolerance (value); unit per @c mz_error_unit_
    double fragment_mz_error_;      ///< Fragment m/z tolerance (value); unit per @c mz_error_unit_
    String mz_error_unit_;          ///< "ppm" or "Da"
    String ion_mode_;               ///< Cached value of parameter @c "ionization_mode"; @c "positive" or @c "negative"; pre-filters DB spectra by precursor-charge sign

    String report_mode_;            ///< Cached value of parameter @c "report_mode"; only @c "top3" or @c "best" are valid — controls whether the top 3 or only the single best match per spectrum is exported to MzTab

    bool merge_spectra_;            ///< If true, merge same-precursor MS2 spectra before searching

    double ccs_error_percent_;      ///< Cached value of parameter @c "ccs_error_percent"; allowed CCS error in percent for filtering candidates (0 = filtering disabled)
  };

}

