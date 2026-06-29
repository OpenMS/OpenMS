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

  /// Single score used for bootstrapping.
  double mbr_score;
};

} // namespace OpenMS::PipEcho
