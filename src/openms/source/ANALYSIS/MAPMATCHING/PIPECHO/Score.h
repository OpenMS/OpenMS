// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

namespace OpenMS::PipEcho
{

/******************************************************************************/
/**
 * How the retention-time feature is fed to the transfer-FDR SVM (and whether it
 * also enters the bootstrap geometric mean). Only meaningful on the local
 * adaptive RT path; the global/baseline path always uses Raw.
 *  - Raw:        the raw |Δrt| (rt_diff_error) is the SVM predictor; RT excluded
 *                from the geometric mean. This is the byte-identical baseline.
 *  - Svm:        a CALIBRATED [0,1] RT-agreement score (rt_score) replaces the raw
 *                rt_diff_error as the SVM predictor; RT still excluded from the mean.
 *  - SvmAndMbr:  as Svm, and the calibrated rt_score is additionally multiplied
 *                into the bootstrap geometric mean (mbr_score).
 */
enum class RtScoreMode
{
  Raw,
  Svm,
  SvmAndMbr
};

/******************************************************************************/
/**
 * Peak score for an acceptor.
 */
struct Score
{
  /// Intensity score calculated from a distribution of intensities.
  double intensity;

  /// Difference between donor and acceptor retention times.
  double rt_diff_error;

  /// Difference (in PPM) between experimental and theoretical mass.
  double mass_error;

  /// Ion-mobility agreement score in [0, 1] (1 = identical mobility), used
  /// only when the data carries ion mobility. A negative value is a sentinel
  /// meaning "ion mobility not applicable" (no IM in the data, or this
  /// donor/acceptor lacks an IM annotation); such a score is excluded from
  /// both the geometric-mean MBR score and the SVM feature set.
  double im_diff_score;

  /// Isotope-distribution agreement score in [0, 1]: the Pearson correlation
  /// between the acceptor's observed isotope envelope and the donor peptide's
  /// theoretical envelope, mapped from [-1, 1] to [0, 1] (1 = identical shape).
  /// 0.5 ("uncorrelated") is used when the envelope is unavailable. Currently an
  /// SVM-only feature (not part of the geometric-mean MBR score).
  double isotope_score;

  /// Calibrated retention-time agreement score in [0, 1] (1 = the acceptor sits
  /// exactly at the locally predicted RT), computed via the FlashLFQ-style
  /// two-tailed CDF against a data-driven RT prediction-error distribution. A
  /// negative value is a sentinel meaning "not calibrated" (the global/baseline
  /// RT path, where the raw rt_diff_error is used instead). Set only on the local
  /// adaptive RT path; see RtScoreMode.
  double rt_score = -1.0;

  /// Single score used for bootstrapping.
  double mbr_score;
};

} // namespace OpenMS::PipEcho
