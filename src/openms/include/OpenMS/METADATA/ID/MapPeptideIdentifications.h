// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <map>
#include <vector>

namespace OpenMS
{
  /**
    @brief The identifications of a feature or consensus map as peptide and protein identifications, apart from the
    map: for formats that hold them (e.g. featureXML, consensusXML) and code that reads them

    The peptide identifications of a feature are those at its path: its index, then the indices of its subordinates
    (a consensus feature has a path of one index). Paths are ordered like the features: a feature before its
    subordinates. See IdentificationDataConverter::exportMapIdentifications() and importMapIdentifications().

    @ingroup Metadata
  */
  struct OPENMS_DLLAPI MapPeptideIdentifications
  {
    std::vector<ProteinIdentification> proteins;
    std::map<std::vector<Size>, PeptideIdentificationList> features;
    PeptideIdentificationList unassigned;

    /// The peptide identifications of the feature at @p path (none if it has none)
    const PeptideIdentificationList& at(const std::vector<Size>& path) const;

    bool empty() const
    {
      return proteins.empty() && features.empty() && unassigned.empty();
    }

    /// All peptide identifications: those of the features (in the order of their paths), then the unassigned ones
    std::vector<const PeptideIdentification*> peptides() const;

    /// Remove the peptide identifications of the feature at @p path and of its subordinates
    void eraseFeature(const std::vector<Size>& path);
  };
} // namespace OpenMS
