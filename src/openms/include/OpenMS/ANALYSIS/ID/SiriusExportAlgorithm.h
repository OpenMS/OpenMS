// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka, Axel Walter $
// $Authors: Oliver Alka, Lukas Zimmermann $
// --------------------------------------------------------------------------

#pragma once 

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>

#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{
  /**
    @brief Algorithm class behind the @ref TOPP_SiriusExport tool: writes a SIRIUS @c .ms input file from MS data and (optionally) a featureXML.

    Drives the export pipeline that feeds SIRIUS / CSI:FingerID / AssayGeneratorMetabo:

      -# preprocessing() — optional: load a featureXML, filter features by
         number of mass traces, build a @ref KDTreeFeatureMaps over them, and
         assign every MS2 spectrum to its closest feature by precursor (RT, m/z).
      -# logFeatureSpectraNumber() — log how many features / additional MS2 spectra
         are going to be processed.
      -# run() — write the SIRIUS @c .ms file and the optional compound-info TSV.

    Configurable behaviour is exposed through @ref DefaultParamHandler defaults
    (see the constructor). The accessors below resolve those parameters by name
    for callers that prefer typed getters over Param lookups.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI SiriusExportAlgorithm : public DefaultParamHandler
    {
    public:
      /// Default constructor; installs the export parameters (see class docs)
      SiriusExportAlgorithm();

      /// @c true if only feature-mapped MS2 spectra should be written (parameter @c feature_only)
      bool isFeatureOnly() const { return param_.getValue("feature_only").toBool(); }
      /// Minimum number of mass traces a feature must have to be retained (parameter @c filter_by_num_masstraces)
      UInt getFilterByNumMassTraces() const { return param_.getValue("filter_by_num_masstraces"); }
      /// Half-width of the precursor m/z window used during MS2-to-feature mapping (unit per @ref precursorMzToleranceUnitIsPPM)
      double getPrecursorMzTolerance() const { return param_.getValue("precursor_mz_tolerance"); }
      /// Half-width of the precursor RT window (seconds) used during MS2-to-feature mapping
      double getPrecursorRtTolerance() const { return param_.getValue("precursor_rt_tolerance"); }
      /// @c true if @ref getPrecursorMzTolerance is interpreted as ppm; otherwise as Th
      bool precursorMzToleranceUnitIsPPM() const { return param_.getValue("precursor_mz_tolerance_unit") == "ppm"; }
      /// @c true to suppress mass-trace-derived isotope-pattern hints in the @c .ms output (parameter @c no_masstrace_info_isotope_pattern)
      bool isNoMasstraceInfoIsotopePattern() const { return param_.getValue("no_masstrace_info_isotope_pattern").toBool(); }
      /// Number of SIRIUS isotope-pattern iterations to record in the @c .ms output (parameter @c isotope_pattern_iterations)
      int getIsotopePatternIterations() const { return  param_.getValue("isotope_pattern_iterations"); }

      /**
        @brief Build the feature/MS2 mapping consumed by run() from a featureXML.

        If @p featureinfo is empty the call is a no-op (both output arguments stay
        untouched) — in that case run() will treat every MS2 spectrum as unassigned.

        If @p featureinfo is given but @c isFeatureOnly() is false, the mass-trace
        filter is silently clamped to 1 (with a warning) so that the adduct annotation
        on multi-trace features is preserved for MS2 spectra. Features with fewer than
        @c getFilterByNumMassTraces() mass traces are then dropped, the survivors are
        wrapped in a @ref KDTreeFeatureMaps, and MS2 spectra are matched to features
        via @ref FeatureMapping::assignMS2IndexToFeature using the configured
        precursor tolerances.

        @param[in]  featureinfo           Path to featureXML; if empty, the function is a no-op.
        @param[in]  spectra               Spectra container (MS1 + MS2); only MS2 are mapped.
        @param[out] feature_mapping_info  Receives the loaded feature maps and the kd-tree built from them.
        @param[out] feature_ms2_indices   Receives the MS2-to-feature assignment / unassigned-MS2 list.
        @throws Exception::FileEmpty If @p featureinfo is non-empty but the file does not exist or is empty.
      */
      void preprocessing(const String& featureinfo,
                               const MSExperiment& spectra,
                               FeatureMapping::FeatureMappingInfo& feature_mapping_info,
                               FeatureMapping::FeatureToMs2Indices& feature_ms2_indices) const;

      /**
        @brief Emit @c OPENMS_LOG_INFO messages summarising how many features / MS2 spectra will be processed.

        Three modes depending on whether a featureXML is in play:
          - @c isFeatureOnly() and @p featureinfo non-empty: logs assigned-feature count only.
          - @p featureinfo non-empty (but not feature-only): logs both assigned-feature count
            and the number of additional unassigned MS2 spectra.
          - @p featureinfo empty: logs the raw MS2 spectrum count in @p spectra.

        @param[in] featureinfo         Path to featureXML (controls which log branch is taken).
        @param[in] feature_ms2_indices Output of @ref preprocessing — feature/MS2 mapping.
        @param[in] spectra             Spectra container used when no featureXML is given.
      */
      void logFeatureSpectraNumber(const String& featureinfo,
                                   const FeatureMapping::FeatureToMs2Indices& feature_ms2_indices,
                                   const MSExperiment& spectra) const;

      /**
        @brief Export the SIRIUS @c .ms input file (and optional compound-info TSV).

        Iterates over @p mzML_files (and, if non-empty, the same-index entry of
        @p featureXML_files for feature-level annotation), runs @ref preprocessing
        internally, and writes one SIRIUS @c .ms file at @p out_ms. If
        @p out_compoundinfo is non-empty, also writes a per-compound TSV with
        provenance/metadata.

        @param[in] mzML_files       Input mzML paths to export.
        @param[in] featureXML_files Optional per-mzML featureXML paths (parallel to @p mzML_files, or empty).
        @param[in] out_ms           Destination SIRIUS @c .ms file (overwritten).
        @param[in] out_compoundinfo Destination compound-info TSV (overwritten); skipped if empty.
        @throws Exception::UnableToCreateFile If @p out_ms cannot be opened for writing.
      */
      void run(const StringList& mzML_files,
               const StringList& featureXML_files,
               const String& out_ms,
               const String& out_compoundinfo) const;

    };
} // namespace OpenMS
