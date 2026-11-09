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
#include <functional>
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
    /// TOPAZ only: whether the MS1 traces of this sample were observed (any MS1 spectrum in its RT window).
    bool has_ms1 = false;
    /// TOPAZ multiple-instance training only: candidates sharing a bag ID (precursor) form one bag with one label.
    std::string bag_id;
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
    /// TOPAZ only, optional: non-owning float32 [samples,3,32] MS1 traces (M+0, M+1, M+2) over the same RT bins.
    const float* ms1 = nullptr;
    std::size_t ms1_count = 0;
    /// TOPAZ only, optional: non-owning float32 [samples, feature_dim] heuristic features, standardised by the caller.
    /// They enter the TOPAZ head next to the trace embeddings (as in redeem's CandidateScorer).
    const float* features = nullptr;
    std::size_t feature_dim = 0;
    /// TOPAZ augmentation only, optional: writes all feature_dim standardised columns for an augmented [6,32,8] window
    /// built from sample `sample`. Only the columns marked in `recomputable` are kept; the others are set to 0.
    std::function<void(const float* window, std::size_t sample, float* features)> recompute_features;
    std::vector<bool> recomputable; ///< Size feature_dim when recompute_features is set.
  };

  /**
    @brief Network architecture.

    GRID2D is the original 6141-parameter 2D CNN over the [6,32,8] RT x IM grid.
    TOPAZ follows redeem's TOPAZ scorer: 1D convolutional branches over the fragment XIC
    (IM summed), the MS1 isotope traces and the fragment mobilogram (core RT bins summed),
    each with an explicit co-elution head, fused by an MLP. The MS1 and mobilogram branch
    outputs are multiplied by has_ms1 and has_im, so a run without them adds exactly zero.
  */
  enum class CNNArchitecture { GRID2D, TOPAZ };

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
    int outer_folds = 0; ///< Zero for the original split, otherwise 2..10 outer family folds.
    int heldout_fold = -1; ///< Outer fold excluded from both fitting and early stopping.
    CNNArchitecture architecture = CNNArchitecture::GRID2D;
    /// TOPAZ only: masked-max multiple-instance training over bag_id bags instead of per-window labels.
    bool mil = false;
    int mil_bags_per_epoch = 40000; ///< MIL: training bags drawn per epoch; also the fixed validation subsample size.
    int topaz_width = 32; ///< TOPAZ only: base width w (conv w-2w-2w, embeddings 2w, MLP 4w-2w); 32 is the TOPAZ layout.
    double weight_decay = 0.0; ///< TOPAZ only: decoupled (AdamW) weight decay; 0 is plain Adam.
    /// TOPAZ only: redeem-topaz training defaults. Per-sample block max normalisation of each trace block
    /// (normalize_max), cosine warm-up LR (5%, floor 0.05), gradient clipping at norm 5, and auxiliary losses
    /// 0.1 x pairwise positive-negative softplus and, with MIL, 0.1 x in-bag ranking (margin 1). Saved in the checkpoint.
    bool topaz_faithful = false;
    /// TOPAZ seed training only (not MIL): each epoch, every training window is replaced with this probability by an
    /// augmented copy s x window + partner, where the partner is a random decoy window of a training family (real
    /// background and interference) and s is log-uniform in [augment_scale_min, 1]. MS1 traces are mixed alike. Heuristic
    /// features are recomputed where recomputable and set to 0 otherwise; with the same probability, unaugmented windows
    /// get the same feature dropout, so missing features do not mark augmented windows. Both classes are treated alike,
    /// so the augmentation carries no label; it simulates weaker, interfered signal. 0 disables (the default).
    double augment_probability = 0.0;
    double augment_scale_min = 0.1;
  };

  /**
    @brief Deterministic outer fold shared by all windows of a peptide/decoy family.

    @param[in] group_id Nonempty peptide-family identifier.
    @param[in] seed Nonnegative partition seed.
    @param[in] folds Number of outer folds, from 2 to 10.
    @return Zero-based outer fold, independent of labels and raw intensities.
    @throws std::invalid_argument If an identifier, seed or fold count is invalid.
  */
  int traceCNNFold(const std::string& group_id, int seed, int folds);

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
    std::size_t augmented_first_epoch = 0; ///< TOPAZ augmentation: windows augmented in the first epoch.
  };

  /**
    @brief Train the raw-trace CNN and save its best validation state.

    Families are partitioned 60/20/20 using FNV1a64 over decimal seed, NUL and
    group ID. Unknown labels are excluded from loss but retained in provenance.
    The checkpoint records the split algorithm, preprocessing, geometry, training
    options and LibTorch version. Loading requires the same LibTorch version.
    With outer_folds enabled, the heldout family fold is excluded from fitting
    and early stopping. Remaining families use the original 60/40 train/validation
    partition. The checkpoint preserves this outer partition contract.
    With the TOPAZ architecture and outer folds, only held-out rows are scored and
    the other logits are NaN. MIL training (options.mil) ignores labels and fits
    one bag per nonempty bag_id (target 1, decoy 0), scored by its best candidate.

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
