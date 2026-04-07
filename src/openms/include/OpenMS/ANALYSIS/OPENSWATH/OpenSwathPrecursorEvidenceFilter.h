// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Prefilter transition-library precursors by quick raw-data evidence.

    This module scans MS1 and/or SWATH MS2 spectra directly through the OpenSWATH
    spectrum access layer. It keeps a target precursor when the run contains
    precursor m/z evidence in MS1 or enough distinct library fragment hits in MS2.
    The output experiment is target-only and can be used as a data-supported
    auto-iRT sampling pool or as the target selection step for a caller that
    adds matching decoys before downstream scoring.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI OpenSwathPrecursorEvidenceFilter :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Compact evidence summary for one target precursor candidate.
    struct PrecursorEvidence
    {
      /// Compound/peptide identifier in the transition experiment.
      String compound_id;
      /// Peptide sequence, if present in the transition experiment.
      String sequence;
      /// Precursor m/z used for matching.
      double precursor_mz{0.0};
      /// Library precursor ion mobility, or -1 if unavailable.
      double precursor_im{-1.0};
      /// Whether MS1 precursor evidence was observed.
      bool supported_ms1{false};
      /// Whether MS2 fragment evidence was observed.
      bool supported_ms2{false};
      /// Number of MS1 peak matches.
      Size ms1_hit_count{0};
      /// Maximum matched MS1 peak intensity.
      double ms1_max_intensity{0.0};
      /// Sum of matched MS1 peak intensities.
      double ms1_sum_intensity{0.0};
      /// RT of the strongest matched MS1 peak, or -1 if unsupported.
      double ms1_best_rt{-1.0};
      /// Total number of distinct MS2 fragment matches across supporting spectra.
      Size ms2_hit_count{0};
      /// Best number of distinct MS2 fragment hits in a single spectrum.
      Size ms2_best_fragment_hits{0};
      /// Maximum matched MS2 peak intensity.
      double ms2_max_intensity{0.0};
      /// Sum of matched MS2 peak intensities from supporting spectra.
      double ms2_sum_intensity{0.0};
      /// RT of the strongest supporting MS2 spectrum, or -1 if unsupported.
      double ms2_best_rt{-1.0};
    };

    /// Result of filtering one run against one transition experiment.
    struct Result
    {
      /// Target-only filtered transition experiment for auto-iRT sampling.
      OpenSwath::LightTargetedExperiment filtered_targets;
      /// Per-candidate evidence, in the same order as the internal target candidates.
      std::vector<PrecursorEvidence> evidence;
      /// Total number of non-decoy precursor candidates considered.
      Size total_target_precursors{0};
      /// Number of candidates supported by the configured evidence source(s).
      Size supported_precursors{0};
      /// Number of candidates with MS1 evidence.
      Size ms1_supported{0};
      /// Number of candidates with MS2 evidence.
      Size ms2_supported{0};
      /// Number of candidates with both MS1 and MS2 evidence.
      Size hybrid_supported{0};
      /// Human-readable summary for logging.
      String summary;
      /// Scale factor applied to precursor ion mobility values for matching and filtered output.
      double precursor_im_scale{1.0};
    };

    /** @name Constructors and Destructors
    */
    //@{
    /// Default constructor
    OpenSwathPrecursorEvidenceFilter();

    /// Destructor
    ~OpenSwathPrecursorEvidenceFilter() override = default;
    //@}

    /// Synchronize members with the parameter object.
    void updateMembers_() override;

    /**
      @brief Filter target precursors by MS1/MS2 evidence observed in the current run.

      @param[in] swath_maps Raw SWATH maps for the current run
      @param[in] transition_exp Input transition experiment
      @param[in] ms1_params MS1 extraction parameters used as precursor m/z/IM tolerances
      @param[in] ms2_params MS2/iRT extraction parameters used as fragment m/z/IM tolerances
      @param[in] pasef Whether the run uses PASEF windows with ion mobility boundaries
      @param[in] threads Number of OpenMP worker threads to use (1 = serial, -1 = OpenMP maximum)
      @return Filtered target experiment and compact evidence summaries

      @throws Exception::IllegalArgument if the configured evidence source cannot be used or too few precursors remain
    */
    Result filter(const std::vector<OpenSwath::SwathMap>& swath_maps,
                  const OpenSwath::LightTargetedExperiment& transition_exp,
                  const ChromExtractParams& ms1_params,
                  const ChromExtractParams& ms2_params,
                  bool pasef,
                  int threads = 1) const;

private:
    bool enabled_{false};
    String evidence_sources_{"hybrid"};
    Size ms1_top_peaks_per_spectrum_{1000};
    Size ms2_top_peaks_per_spectrum_{1000};
    Size ms2_top_transitions_per_precursor_{3};
    Size ms2_min_fragment_hits_{2};
    Size min_supported_precursors_{3};
    bool peak_picking_enabled_{false};
    bool peak_picking_use_gauss_{true};
  };
}
