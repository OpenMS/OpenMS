// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <string>

namespace OpenMS
{
/**
  @brief Pseudo features of MetaProSIP: identifications and spectra that feature finding gave no feature.

  Besides the features of the feature map, MetaProSIP evaluates the unassigned identifications (@c -use_unassigned_ids)
  and the MS2 spectra without identification (@c -use_averagine_ids) at their own positions. Their pseudo features
  carry the given type as meta value @c feature_type.
*/
class MetaProSIPPseudoFeatures
{
public:
  /**
    @brief Add a pseudo feature for every unassigned identification of @p features whose first match has a charge

    The pseudo feature lies at the RT of the identification and at the m/z of the peptide of its first match, and links
    the matches of the identification.

    @return The number of pseudo features added
  */
  static Size addUnassignedIdentifications(FeatureMap& features, const std::string& feature_type);

  /**
    @brief Add a pseudo feature for every MS2 spectrum of @p ms2 with a precursor that no identification of @p features
    was made from

    The identifications with matches, linked by a feature or unassigned, claim a spectrum each: the first one from
    their RT on (minus 1e-5 s). A pseudo feature lies at the RT of its spectrum and at the m/z of the spectrum's
    precursor, with its charge, and links nothing: MetaProSIP models it as an averagine peptide.

    @param[in,out] features Features with their identifications
    @param[in] ms2 MS2 spectra, sorted by RT
    @param[in] feature_type The meta value @c feature_type of the pseudo features
    @return The number of pseudo features added
  */
  static Size addUnidentifiedSpectra(FeatureMap& features, const PeakMap& ms2, const std::string& feature_type);
};
} // namespace OpenMS
