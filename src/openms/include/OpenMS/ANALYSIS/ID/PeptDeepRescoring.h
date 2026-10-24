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

#include <map>
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
    errors propagate into the other. Only the best-scoring hit of each spectrum, and only
    if it is not a decoy, is a candidate, so the calibration does not depend on how many
    hits per spectrum the search reports: with ten, lower-ranked candidates and decoys
    would otherwise make up most of the better-scoring half.

    Parameters are not pulled into the generated documentation with
    @c \@htmlinclude: that requires the class to be registered in
    DefaultParamHandlerDocumenter, which cannot instantiate a type that only
    exists in ONNX-enabled builds. ProSEAlgorithm omits it for the same reason.

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

    /// Name of the run meta value ("peptdeep_nce", a double) in which annotate() stores, on each annotated
    /// ProteinIdentification, the normalised collision energy its predictions used. It is written to idXML with the
    /// run, so consumers of stored results can rely on the name.
    static constexpr const char* NCE_META_VALUE = "peptdeep_nce";

    /**
      @brief Annotates @p peptide_ids in place and registers the features for rescoring.

      @param[in] spectra Spectra the identifications were made on; supplies retention
             time and, for the NCE scan, the instrument's reported collision energy.
      @param[in,out] protein_ids Search metadata; the feature names are appended to the
             @p extra_features of each run that has PSMs, and the NCE its predictions used is
             stored on the run as the double meta value NCE_META_VALUE (see below).
      @param[in,out] peptide_ids PSMs to annotate. Every hit receives every feature, so
             that no run ends up with a feature present on only some of its PSMs.

      Multiple runs in one call are fine: PSMs are partitioned by
      PeptideIdentification::getIdentifier() and each run gets its own NCE scan and its
      own RT calibration.

      Each annotated ProteinIdentification gets the meta value NCE_META_VALUE
      ("peptdeep_nce", a double):
      the NCE its run's predictions used -- the one selected by the scan, or @p nce when
      that is set. It lets a consumer of the result predict spectra at the same energy the
      annotation used. A ProteinIdentification without PSMs is not annotated and gets none.
      If no run identifier matches any ProteinIdentification, the features go on the first
      one as before, and so does the value, but only when all runs used the same NCE.

      @note That partition is only as good as the run identifiers, so call this
            <em>before</em> any step that collapses several runs into one -- ProSE's
            @p out_merged and IDMergerAlgorithm, which remap every identifier onto a
            single merged run, or @ref TOPP_IDMerger with @p merge_proteins_add_PSMs
            (its default appends runs and is therefore safe). Both calibrated quantities
            genuinely differ between runs: six files of one HLA ligand atlas dataset
            (PXD019643) calibrated to median RT residuals between 98 s and 625 s. On
            already-collapsed input this silently fits one NCE and one RT curve across
            all of them, and nothing here can detect it after the fact -- the merged
            input is indistinguishable from a genuine single run.

      @throws Exception::FileNotFound if a configured model file does not exist.
      @throws Exception::MissingInformation if no PSM carries peak annotations, since
              every MS2 feature would otherwise be silently zero.
    */
    void annotate(const PeakMap& spectra,
                  std::vector<ProteinIdentification>& protein_ids,
                  PeptideIdentificationList& peptide_ids);

    /**
      @brief Selects the normalised collision energy of each run of @p peptide_ids as annotate() does, without
      predicting the spectra of all PSMs or adding features.

      For a caller that predicts spectra itself at the energy annotate() would use (e.g. ProSE's -qpx_extended for
      runs that were not rescored with PeptDeep): the same selection on the same sample of confident PSMs (the scan of
      the grid, or @p nce if that is set), and the NCE stored on the runs as NCE_META_VALUE, exactly as annotate()
      stores it. Needs only @p ms2_model. Logs no "[PeptDeepRescoring]" lines of its own besides warnings about the
      input (shared with annotate()); its progress is that of the scan ("Selecting collision energy...").

      @param[in] spectra As for annotate()
      @param[in,out] protein_ids The runs; each run with PSMs gets NCE_META_VALUE
      @param[in] peptide_ids The PSMs, as for annotate(); not changed

      @throws Exception::FileNotFound if @p ms2_model does not exist.
      @throws Exception::MissingInformation if @p ms2_model is not set, or no PSM carries peak annotations.
    */
    void selectNCE(const PeakMap& spectra,
                   std::vector<ProteinIdentification>& protein_ids,
                   const PeptideIdentificationList& peptide_ids);

    /// Normalised collision energy used for the run annotated last (after any scan).
    /// -1 if annotate() has not run, or returned before selecting one.
    double getUsedNCE() const;

    /// Median absolute RT residual over the calibration set of the run annotated last.
    /// -1 if annotate() has not run, or could not calibrate.
    double getRTCalibrationError() const;

    /// Number of PSMs (peptide hits) seen by the last annotate() call.
    Size getNumPSMs() const;

    /// Number of PSMs of the last annotate() call that received a model prediction
    /// and carry peak annotations, i.e. whose MS2 features are computed rather than
    /// zero-filled. 0 if annotate() returned early or predicted nothing.
    Size getNumPredictedPSMs() const;

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
      bool calibration_candidate; ///< best hit of its spectrum and not a decoy
    };

    /// The PSMs of annotate() and selectNCE(): a row per hit with a sequence, the distinct (sequence, charge) pairs to
    /// predict, and the rows of each run (by run identifier, in a stable order)
    struct PSMRows
    {
      std::vector<Row> rows;
      std::vector<std::string> uniq_seq;
      std::vector<float> uniq_charge;
      std::map<std::string, std::vector<Size>> rows_by_run;
    };

    /// The PSMs of @p peptide_ids. Throws Exception::MissingInformation if no PSM carries peak annotations, and warns
    /// if only some runs do.
    static PSMRows gatherPSMs_(const PeptideIdentificationList& peptide_ids);

    /// The collision energy of each MS2 spectrum of @p spectra by native ID (native IDs recorded with differing
    /// energies are left out, with a warning)
    static std::map<std::string, double> collisionEnergies_(const PeakMap& spectra);

    /// The confident rows of a run (@p run_rows of @p rows): its calibration candidates (all rows if it has none),
    /// best search score first, without the calibration_quantile held out
    std::vector<Size> confidentRows_(const std::vector<Row>& rows, const std::vector<Size>& run_rows, bool higher_better) const;

    /// The NCE of a run: 'nce' if set, else the grid scan on a sample of its @p confident rows. With the scan, @p centre
    /// is the centre of its grid and @p median_cosine the median cosine of the selected NCE (else both unchanged).
    double selectRunNCE_(const std::map<std::string, double>& collision_energies,
                         const PeptideIdentificationList& peptide_ids,
                         const std::vector<Row>& rows,
                         const std::vector<Size>& run_rows,
                         const std::vector<Size>& confident,
                         const std::vector<std::string>& uniq_seq,
                         const std::vector<float>& uniq_charge,
                         PeptDeepMS2Inference& ms2_inference,
                         double& centre,
                         double& median_cosine);

    /// Stores the NCE of each run of @p nce_by_run on the run of @p protein_ids with its identifier (NCE_META_VALUE);
    /// if none has one of these identifiers, on the first run, provided that all NCEs agree
    static void storeNCE_(std::vector<ProteinIdentification>& protein_ids, const std::map<std::string, double>& nce_by_run);

    /// Annotates the PSMs of a single identification run. Each run gets its own
    /// collision energy and its own retention-time calibration, since neither
    /// transfers across runs with different gradients or acquisition settings.
    void annotateRun_(const std::string& run_id,
                      const std::map<std::string, double>& collision_energies,
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
    Size batch_size_{64};
    Int threads_{4};

    double used_nce_{-1.0};
    double rt_calibration_error_{-1.0};
    Size num_psms_{0};
    Size num_predicted_psms_{0};
  };

} // namespace OpenMS
