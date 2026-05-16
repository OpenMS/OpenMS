// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka, Axel Walter $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <fstream>

namespace OpenMS
{

  /**
    @brief Writes a SIRIUS @c .ms file from an MSExperiment, optionally enriched with feature/adduct/formula annotations.

    Used by @ref TOPP_SiriusExport to translate centroided MS2 data (and, optionally, feature
    annotations from @ref FeatureFindingMetabo, @ref TOPP_MetaboliteAdductDecharger, and/or
    @ref TOPP_AccurateMassSearch) into the @c .ms compound format consumed by SIRIUS.

    The writer chooses one of three layouts based on the input @c FeatureMapping::FeatureToMs2Indices —
      - @b feature-driven (@c assignedMS2 non-empty): one compound block per feature, carrying
        its associated MS2 spectra plus adduct / formula / description metadata if those are
        present on the feature.
      - @b unassigned-MS2 (only when @c unassignedMS2 is non-empty and @p feature_only is @c false):
        an additional compound block per unassigned MS2 spectrum with @c UNKNOWN description /
        sumformula / adducts.
      - @b no-feature-information (both maps empty): every MS2 spectrum in @p spectra is
        emitted as its own UNKNOWN compound. This is the mzML-only fallback.
    The first two layouts can both fire in the same call (feature-driven followed by
    unassigned-MS2). For each compound emitted, a matching @ref CompoundInfo entry is
    appended to @p v_cmpinfo for downstream mzTab-M export.

    Constraints enforced during @ref store —
      - Spectra must be centroided. If the first spectrum is @c SpectrumSettings::PROFILE,
        @c OpenMS::Exception::IllegalArgument is thrown.
      - The input mzML must carry @c SourceFile annotation; an empty @ref MSExperiment::getSourceFiles
        also throws @c OpenMS::Exception::IllegalArgument.
      - SIRIUS supports only singly charged precursors: features (and spectra) with
        @c |charge| > 1 are dropped and counted; a per-call summary is written via
        @c OPENMS_LOG_WARN at the end of @ref store.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI SiriusMSFile
  {
public:

  /// Source-file accession metadata of the input mzML, captured for the mzTab-M @c MS_RUN section.
  class AccessionInfo
  {
  public:
    String sf_path; ///< sourcefile path for mztab-m
    String sf_type; ///< sourcefile type for mztab-m
    String sf_filename; ///< sourcefile name for mztab-m
    String sf_accession; ///< sourcefile accessions for mztab-m
    String native_id_accession; ///< nativeID accession for mztab-m
    String native_id_type; ///< nativeID type for mztab-m
  };

  /**
    @brief Per-compound metadata accumulated while writing the @c .ms file.

    One entry is appended for every compound block emitted by @ref store. The same data is
    later serialised via @ref saveFeatureCompoundInfoAsTSV and used by downstream tools to
    map SIRIUS results back to their originating spectra / features and to populate the
    mzTab-M @c SmallMolecule / @c SmallMoleculeFeature sections.
  */
  class CompoundInfo
  {
  public:
    String cmp; ///< query_id used compound in .ms file
    double pmass; ///< parent/precursor mass of the compound
    double pint_mono; ///< parent/precursor intensity of the compound
    double rt; ///< retention time of the compound
    double fmz; ///< annotated mass of a feature (if available)
    String fid; ///< annotated feature_id (if available)
    String formula; ///< sumformula of the compound
    int charge; ///< precursor/feature charge
    String ionization; ///< adduct information
    String des; ///< description/name of the compound
    String specref_format; ///< spectra ref format for mztab-m
    String source_file; ///< sourcefile for mztab-m
    String source_format; ///< format of the sourcefile for mztab-m
    std::vector<String> native_ids; ///< native ids of the associated spectra
    String native_ids_id; ///< concatenated list of the associated spectra
    std::vector<String> m_ids; ///< native ids and identifier for multiple possible identification via AMS ("|" separator)
    String m_ids_id; ///< concatenated list of native ids and identifier for multiple possible identification via AMS ("|" separator) used for mapping of compounds and the annotated spectrum.
    std::vector<String> scan_indices; ///< index of the associated spectra
    std::vector<String> specrefs; ///< spectra reference for mztab-m
    int file_index; ///< source file index >
  };

  /**
    @brief Write one mzML/featureXML pair to a SIRIUS @c .ms output stream.

    Selects the @c .ms layout (feature-driven, unassigned-MS2, or no-feature-information)
    according to the contents of @p feature_mapping; see the class documentation for the
    branching rules. For every compound block emitted, a @ref CompoundInfo record is
    appended to @p v_cmpinfo for downstream mzTab-M export. When adduct information is
    missing for a spectrum no adduct line is written — SIRIUS will then assume defaults.

    The accession-CV term for the @c sourcefile is resolved by loading the PSI-MS OBO
    (@c CV/psi-ms.obo via @c File::find) and locating the term whose name matches the
    @c SourceFile type. The output stream @p os is written to but not closed; ownership
    remains with the caller.

    @param[in]     spectra                          Peakmap from the input mzML. The first spectrum must be centroided.
    @param[in]     os                               Open output stream the @c .ms text is written to.
    @param[in]     feature_mapping                  Result of @ref FeatureMapping::assignMS2IndexToFeature; selects the output layout.
    @param[in]     feature_only                     If @c true, MS2 spectra that are not associated with any feature are dropped instead of being emitted as additional UNKNOWN compounds.
    @param[in]     isotope_pattern_iterations       Upper bound on the number of isotope-trace peaks searched per precursor when no feature mass-traces are available.
    @param[in]     no_masstrace_info_isotope_pattern If @c true, fall back to spectrum-based isotope-pattern extraction when feature mass-traces are absent.
    @param[in,out] v_cmpinfo                        Receives one @ref CompoundInfo per emitted compound (appended; existing entries preserved).
    @param[in]     file_index                       Numeric identifier mixed into compound IDs to disambiguate entries from different mzML files.

    @throws OpenMS::Exception::IllegalArgument if @p spectra contains profile data (centroiding is required).
    @throws OpenMS::Exception::IllegalArgument if @p spectra carries no @c SourceFile annotation.

    @note SIRIUS supports only singly charged precursors: features (and spectra) with
          @c |charge| @c > @c 1 are skipped. A summary of the skipped/assumed counts is
          emitted via @c OPENMS_LOG_WARN at the end of the call (even when all counts are
          zero).
  */
    static void store(const MSExperiment& spectra,
                      std::ofstream& os,
                      const FeatureMapping::FeatureToMs2Indices& feature_mapping,
                      const bool& feature_only,
                      const int& isotope_pattern_iterations,
                      const bool no_masstrace_info_isotope_pattern,
                      std::vector<SiriusMSFile::CompoundInfo>& v_cmpinfo,
                      const size_t& file_index);
    /**
      @brief Serialise the @ref CompoundInfo records collected during @ref store as a TSV file.

      The file is written with a fixed 16-column header line in the order:
      @c cmp, @c file_index, @c pmass, @c pint_mono, @c rt, @c fmz, @c fid, @c formula,
      @c charge, @c ionization, @c des, @c specref_format, @c source_file, @c source_format,
      @c native_ids_id, @c m_ids_id. The columns @c native_ids and @c m_ids (the raw vector
      versions) are not written — only their already-concatenated @c _id forms.

      @param[in] v_cmpinfo Compound records to write (typically the vector populated by @ref store).
      @param[in] filename  Destination path of the TSV file.

      @throws std::runtime_error if @p filename cannot be opened for writing.
    */
    static void saveFeatureCompoundInfoAsTSV(const std::vector<SiriusMSFile::CompoundInfo>& v_cmpinfo,
                                      const std::string& filename);

  protected:
    /**
    @brief Internal structure to write the .ms file (called in store function)

    @param[out] os: stream
    @param[in] spectra: spectra
    @param[in] ms2_spectra_index: vector of index ms2 spectra (in feature)
    @param[in] ainfo: accession information
    @param[in] adducts: vector of adducts
    @param[in] v_description: vector of descriptions
    @param[in] v_sumformula: vector of sumformulas
    @param[in] f_isotopes: isotope pattern of the feature
    @param[in] feature_charge: feature charge
    @param[in] feature_id: feature id
    @param[in] feature_rt: features retention time
    @param[in] feature_mz: feature mass to charge
    @param[out] writecompound: bool if new compound should be written in .ms file
    @param[in] no_masstrace_info_isotope_pattern: bool if isotope pattern should be extracted (if not in feature)
    @param[in] isotope_pattern_iterations: number of iterations (trying to find a C13 pattern)
    @param[in] count_skipped_spectra: count number of skipped spectra
    @param[in] count_assume_mono: count number of features where mono charge was assumed
    @param[in] count_no_ms1: count number of compounds without a valid ms1 spectrum
    @param[in] v_cmpinfo: vector of CompoundInfo
    @param[in] file_index: file index (to differentiate entries derived from different mzML files and resolve ambiguities)
    */

    static void writeMsFile_(std::ofstream& os,
                             const MSExperiment& spectra,
                             const std::vector<size_t>& ms2_spectra_index,
                             const SiriusMSFile::AccessionInfo& ainfo,
                             const StringList& adducts,
                             const std::vector<String>& v_description,
                             const std::vector<String>& v_sumformula,
                             const std::vector<std::pair<double,double>>& f_isotopes,
                             int& feature_charge,
                             uint64_t& feature_id,
                             const double& feature_rt,
                             const double& feature_mz,
                             bool& writecompound,
                             const bool& no_masstrace_info_isotope_pattern,
                             const int& isotope_pattern_iterations,
                             int& count_skipped_spectra,
                             int& count_assume_mono,
                             int& count_no_ms1,
                             std::vector<SiriusMSFile::CompoundInfo>& v_cmpinfo,
                             const size_t& file_index);
    /**
      @brief Return the index of the most-intense peak of @p spectrum whose m/z lies within a tolerance window around @p test_mz.

      @param[in] test_mz   Target m/z; the search window is built around this value.
      @param[in] spectrum  Spectrum to scan.
      @param[in] tolerance Half-width of the tolerance window (in Da or ppm, depending on @p ppm).
      @param[in] ppm       If @c true, @p tolerance is interpreted as ppm; otherwise as Da.

      @return Index of the most-intense peak inside the window, or @c -1 if the window
              contains no peak.
    */
    static Int getHighestIntensityPeakInMZRange_(double test_mz,
                                                 const MSSpectrum& spectrum,
                                                 double tolerance,
                                                 bool ppm);

    /**
      @brief Walk an isotope ladder from the precursor m/z, picking the most-intense peak at each C12-C13 step.

      Used by @ref store when no per-feature mass-trace information is available. The
      monoisotopic peak is located inside a fixed 10 ppm window around @p precursor_mz;
      subsequent traces are located at @c precursor_mz @c + @c k @c * @c C13C12_MASSDIFF_U / |charge|
      inside a 1 ppm window, stopping when @p iterations reaches zero or when no peak is
      found at the current step.

      @param[in]     precursor_mz       Precursor m/z to start the ladder from.
      @param[in]     precursor_spectrum Spectrum to search.
      @param[in,out] iterations         Maximum number of isotope-trace steps; decremented per successful step.
      @param[in]     charge             Precursor charge (used to scale the C12-C13 distance; @c 0 disables scaling).

      @return Sequence of picked isotope peaks, starting with the monoisotopic peak (empty if not even that is found).
    */
    static std::vector<Peak1D> extractPrecursorIsotopePattern_(const double& precursor_mz,
                                                               const MSSpectrum& precursor_spectrum,
                                                               int& iterations,
                                                               const int& charge);



  };

} // namespace OpenMS