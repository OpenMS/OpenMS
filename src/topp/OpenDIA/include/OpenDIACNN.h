// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief Numeric acquisition geometry for fixed [6,32,8] trace windows.
  */
  struct CNNGeometry
  {
    double mz_tolerance = 0.02; ///< Fragment m/z tolerance in Da.
    double rt_half_width = 15.0; ///< RT half-width in seconds.
    double im_half_width = 0.03; ///< Ion mobility half-width in Vs/cm^2.
  };

  /**
    @brief Trace metadata, explicit supervision and peptide-family provenance.

    Only library intensities and the mobility availability mask enter the network.
    Fragment channels must have descending library intensity. A zero fragment m/z
    denotes padding and requires zero library and raw intensities in that channel.
  */
  struct CNNSample
  {
    std::string sample_id;
    std::string group_id;
    int label = -1; ///< Unknown (-1), negative (0) or positive (1).
    bool is_decoy = false;
    bool has_im = false;
    std::array<double, 6> library_intensities{};
    std::array<double, 6> fragment_mz{};
  };

  /**
    @brief Non-owning view of contiguous raw float32 [samples,6,32,8] windows.

    The caller must retain the raw storage for the duration of training or scoring.
    Preprocessing allocates only the requested tensor batch, not the full dataset.
  */
  struct CNNData
  {
    const float* intensities = nullptr;
    std::size_t intensity_count = 0;
    std::vector<CNNSample> samples;
  };

  /**
    @brief Deterministic CPU training options for the fixed 6141-parameter CNN.
  */
  struct CNNOptions
  {
    int seed = 42;
    int threads = 5;
    int epochs = 40; ///< Valid range: 1..100.
    int patience = 6; ///< Valid range: 1..20.
    int batch_size = 32;
  };

  /**
    @brief Raw model logits, family partitions and diagnostic training metrics.

    Logits are not identification probabilities and make no FDR guarantee.
    A heldout_auc of -1 means that the heldout rows lack one of the two classes.
  */
  struct CNNResult
  {
    std::vector<float> logits;
    std::vector<std::string> splits;
    int epochs_run = 0;
    double validation_loss = 0;
    double heldout_auc = -1;
  };

  /**
    @brief Train the raw-trace CNN and save its best validation state.

    Families are partitioned 60/20/20 using FNV1a64 over decimal seed, NUL and
    group ID. Unknown labels are excluded from loss but retained in provenance.
    The checkpoint records the split algorithm, preprocessing, geometry, training
    options and LibTorch version. Loading requires the same LibTorch version.

    @param[in] data Raw windows and metadata, retained by the caller.
    @param[in] geometry Positive finite acquisition sampling geometry.
    @param[in] options CPU thread count, random seed and training limits.
    @param[in] checkpoint Destination native LibTorch archive path.
    @return Logits and splits in input order, with best validation diagnostics.
    @throws std::invalid_argument If data, geometry, options or labeled partitions are invalid.
  */
  CNNResult trainTraceCNN(const CNNData& data, const CNNGeometry& geometry,
                          const CNNOptions& options, const std::string& checkpoint);

  /**
    @brief Score raw windows using a compatible native CNN checkpoint.

    Known families retain train/validation/heldout provenance. New families are
    marked external. Training-family overlap is refused unless explicitly allowed.

    @param[in] data Raw windows and metadata, retained by the caller.
    @param[in] geometry Acquisition geometry, which must match the checkpoint.
    @param[in] checkpoint Source native LibTorch archive path.
    @param[in] threads Positive number of CPU execution threads.
    @param[in] allow_training_overlap Allow diagnostic scoring of training families.
    @return Logits and provenance splits in input order.
    @throws std::invalid_argument If inputs, checkpoint metadata or overlap policy are invalid.
  */
  CNNResult scoreTraceCNN(const CNNData& data, const CNNGeometry& geometry,
                          const std::string& checkpoint, int threads = 5,
                          bool allow_training_overlap = false);
}
