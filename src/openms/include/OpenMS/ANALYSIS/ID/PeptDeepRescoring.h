// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <string>
#include <vector>

namespace OpenMS
{
  class PeptideIdentificationList;
  class ProteinIdentification;
  class PeptDeepMS2Inference;
  class PeptDeepRTInference;

  /**
    @brief Annotates PSMs with rescoring features derived from PeptDeep predictions.

    Compares every PSM against an AlphaPeptDeep prediction of the same peptide and
    stores the agreement as per-hit features, which are then registered in the search
    parameters' @p extra_features list so that PercolatorAdapter picks them up:

    - @p ms2_cosine, @p ms2_spectral_angle, @p ms2_pearson, @p ms2_frac_pred_found --
      predicted versus observed fragment intensities
    - @p rt_abs_error -- observed retention time versus a per-run calibration of the
      predicted retention time

    Observed intensities are taken from the PSM's own peak annotations, so no second
    pass over the spectra is needed; the search engine must have been run with
    fragment annotation enabled (ProSE's @p annotate:PSM default of @p ALL does this).
    Predicted ions that were not matched contribute zero on the observed side, which is
    what makes a wrong peptide score badly rather than simply score on fewer ions.

    Two per-run quantities are calibrated from the data rather than being asked of the
    user, because both are properties of the run and neither has a safe global default:

    - <b>Normalised collision energy.</b> The model's NCE scale does not coincide with
      the number the instrument reports, so the mzML value is used only as the centre of
      a small search grid. Each candidate is scored on a sample of confident PSMs and the
      one with the highest median cosine wins. With @p nce set to a positive value the
      scan is skipped and that value is used directly.
    - <b>Retention time.</b> Predicted RT is in the model's own iRT-like units, so it is
      mapped onto the run's observed gradient by fitting a TransformationModel on
      confident PSMs before the residual is taken.

    Confident PSMs are selected by search score rather than by spectral similarity: using
    similarity here would tie the RT feature to the MS2 features and let one feature's
    errors propagate into the other.

    @htmlinclude OpenMS_PeptDeepRescoring.parameters

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI PeptDeepRescoring :
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:
    /// Constructor. Sets the default parameters; the model paths start empty and must be set.
    PeptDeepRescoring();

    /// Destructor.
    ~PeptDeepRescoring() override;

    /// Names of the features this class adds to each PSM, in a stable order.
    static StringList featureNames();

    /**
      @brief Annotates @p peptide_ids in place and registers the features for rescoring.

      @param[in] spectra Spectra the identifications were made on; supplies retention
             time and, for the NCE scan, the instrument's reported collision energy.
      @param[in,out] protein_ids Search metadata; the feature names are appended to the
             first run's @p extra_features.
      @param[in,out] peptide_ids PSMs to annotate. Every hit receives every feature, so
             that no run ends up with a feature present on only some of its PSMs.

      @throws Exception::FileNotFound if a configured model file does not exist.
      @throws Exception::MissingInformation if no PSM carries peak annotations, since
              every MS2 feature would otherwise be silently zero.
    */
    void annotate(const PeakMap& spectra,
                  std::vector<ProteinIdentification>& protein_ids,
                  PeptideIdentificationList& peptide_ids);

    /// Normalised collision energy used for the run annotated last (after any scan).
    /// -1 if annotate() has not run, or returned before selecting one.
    double getUsedNCE() const;

    /// Median absolute RT residual over the calibration set of the run annotated last.
    /// -1 if annotate() has not run, or could not calibrate.
    double getRTCalibrationError() const;

  protected:
    void updateMembers_() override;

    /// Row of the flat PSM work list built by annotate().
    struct Row
    {
      Size pi;      ///< index into the PeptideIdentificationList
      Size hit;     ///< index into that identification's hits
      Size len;     ///< peptide length
      Size key;     ///< index into the de-duplicated (sequence, charge) list
      double score; ///< search score, used to pick the calibration set
    };

    /// Annotates the PSMs of a single identification run. Each run gets its own
    /// collision energy and its own retention-time calibration, since neither
    /// transfers across runs with different gradients or acquisition settings.
    void annotateRun_(const PeakMap& spectra,
                      PeptideIdentificationList& peptide_ids,
                      const std::vector<Row>& rows,
                      const std::vector<Size>& run_rows,
                      const std::vector<std::string>& uniq_seq,
                      const std::vector<float>& uniq_charge,
                      bool higher_score_better,
                      PeptDeepMS2Inference& ms2_inference,
                      PeptDeepRTInference& rt_inference);

    std::string ms2_model_;
    std::string rt_model_;
    double nce_{-1.0};
    double nce_grid_halfwidth_{6.0};
    double nce_grid_step_{2.0};
    double nce_fallback_{30.0};
    Size nce_sample_size_{2000};
    std::string instrument_;
    double calibration_quantile_{0.5};
    std::string rt_model_type_{"b_spline"};
    Size rt_num_nodes_{5};
    int ms2_min_ordinal_{1};
    double ms2_strong_fraction_{0.05};
    Size batch_size_{500};
    Int threads_{4};

    double used_nce_{-1.0};
    double rt_calibration_error_{-1.0};
  };

} // namespace OpenMS
