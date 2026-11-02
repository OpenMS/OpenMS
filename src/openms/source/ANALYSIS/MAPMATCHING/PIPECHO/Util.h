// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/METADATA/PeptideHit.h>

#include <optional>
#include <vector>

namespace OpenMS::PipEcho::Util
{

/******************************************************************************/
/**
 * Return the first peptide hit from a feature.
 *
 * This code is here to isolate checking hits because the process
 * might change in a future version of OpenMS.
 */
std::optional<PeptideHit> feature_hit(const Feature&);

/******************************************************************************/
/**
 * Return `true` if the peptide hit referenced by the given feature is
 * a decoy.
 */
bool feature_is_decoy(const Feature&);

/******************************************************************************/
/**
 * Compute the mass error in PPM for the given feature.
 */
std::optional<double> feature_mass_error(const Feature&);

/******************************************************************************/
/**
 * Return the feature's observed isotope envelope (per-isotope intensities), or
 * nullopt if absent. Reads the "pipecho_obs_envelope" meta value snapshotted by
 * ProteomicsLFQ from the FeatureFinderIdentification subordinate intensities
 * before they are stripped.
 */
std::optional<std::vector<double>> feature_obs_envelope(const Feature&);

/******************************************************************************/
/**
 * Return the feature's ion-mobility value (e.g. 1/K0), or nullopt if the
 * feature carries no ion-mobility annotation.
 *
 * Reads the "IM_median" meta value written by FeatureFinderIdentification
 * (and kept through ProteomicsLFQ), falling back to the generic
 * Constants::UserParam::IM key used by other detectors (e.g. Biosaur2).
 */
std::optional<double> feature_ion_mobility(const Feature&);

/******************************************************************************/
/**
 * Return the feature's ion-mobility spread, i.e. "IM_max" - "IM_min", or
 * nullopt if those annotations are absent or non-positive. Used to derive a
 * data-driven ion-mobility tolerance.
 */
std::optional<double> feature_im_width(const Feature&);

} // namespace OpenMS::PipEcho::Util
