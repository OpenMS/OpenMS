// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepCCSInference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepRTInference.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/config.h>

#include <memory>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Input precursor for native PeptDeep spectral-library prediction.

    The caller owns precursor identity. In particular, @c id must be unique across the
    complete library, not only within one prediction batch. This keeps IDs stable when
    the same predictor is later used for bounded/on-demand library generation.
  */
  struct OPENMS_DLLAPI PeptDeepLibraryPrecursor
  {
    AASequence peptide;
    std::string id;
    Int charge = 0;
    float nce = 30.0f;
    Int64 instrument_index = 0;
    std::vector<std::string> protein_refs;
  };

  /**
    @brief Convert peptide precursors into an OpenSWATH assay batch using native PeptDeep ONNX models.

    Predicts retention time and MS2 fragment intensities for every precursor. CCS prediction is
    optional; when enabled, predicted CCS values are converted to reduced inverse mobility (1/K0)
    through IMTypes and stored as the compound drift time and transition precursor IM.

    MS2 output is mapped using the AlphaPeptDeep channel order
    @c [b_z1,b_z2,y_z1,y_z2,b_modloss_z1,b_modloss_z2,y_modloss_z1,y_modloss_z2].
    Only the four canonical non-loss channels are emitted in this initial implementation.

    This class deliberately does not perform digestion, assay filtering, decoy generation, or
    OpenSWATH ID normalization. Those remain responsibilities of the surrounding library workflow.

    Available in builds with @c WITH_ONNX.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI PeptDeepLibraryPredictor
  {
  public:
    struct OPENMS_DLLAPI Config
    {
      std::string rt_model_path;
      std::string ccs_model_path;
      std::string ms2_model_path;
      int intra_op_threads = 4;
      Size batch_size = 500;
      bool predict_ccs = true;
    };

    /// Construct with the installed OpenMS PeptDeep models and default inference settings.
    PeptDeepLibraryPredictor();

    /// Construct with explicit inference/model configuration.
    explicit PeptDeepLibraryPredictor(const Config& config);

    /**
      @brief Predict one batch of precursor assays.

      @param precursors Precursors with caller-assigned unique IDs.
      @return A LightTargetedExperiment containing one compound per precursor and four b/y
              fragment channels per cleavage position.
      @throws Exception::IllegalArgument for empty input, empty/duplicate IDs, non-positive
              precursor charges, peptide sequences shorter than two residues, or empty protein IDs.
      @throws Exception::InvalidValue if an inference output has an unexpected shape.
    */
    OpenSwath::LightTargetedExperiment predict(const std::vector<PeptDeepLibraryPrecursor>& precursors);

  private:
    PeptDeepRTInference rt_predictor_;
    std::unique_ptr<PeptDeepCCSInference> ccs_predictor_;
    PeptDeepMS2Inference ms2_predictor_;
  };
} // namespace OpenMS
