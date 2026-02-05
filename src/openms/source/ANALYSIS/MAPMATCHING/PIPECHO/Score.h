// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

namespace OpenMS {
namespace PipEcho {

  /****************************************************************************/
  /**
   * Peak score for an acceptor.
   */
  struct Score {
    /// Intensity score calculated from a distribution of intensities.
    double intensity;

    /// Difference between donor and acceptor retention times.
    double rt_diff_error;

    /// Difference (in PPM) between experimental and theoretical mass.
    double mass_error;

    /// Single score used for bootstrapping.
    double mbr_score;
  };

}}
