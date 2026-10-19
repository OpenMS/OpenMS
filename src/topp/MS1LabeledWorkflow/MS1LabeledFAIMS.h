// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <map>
#include <optional>

namespace OpenMS
{
class IDMapper;
class FeatureGroupingAlgorithmQT;

/** @brief Keep FAIMS acquisitions separate while retaining physical run/channel columns. */
class MS1LabeledFAIMS
{
public:
  /// A missing CV denotes conventional, non-FAIMS data; NaN is never used as a map key.
  using CV = std::optional<double>;

  /// Recover identification CVs from spectrum references before any RT-based reference repair.
  static void annotateCompensationVoltages(const MSExperiment& spectra, PeptideIdentificationList& ids);

  /**
    @brief Map IDs onto multiplets through the elution regions of their channels, only at the same CV.

    @p channels holds one feature per (multiplet, channel) of @p map, with the unique id of the channel's
    feature handle and its convex hulls (see the workflow). The IDs are matched against these hulls
    (IDMapper's feature map mode, i.e. against where the channel elutes rather than its apex) and every
    match is attached to the channel's multiplet, with the channel as @c map_index. An ID matching
    several channels of one multiplet is attached once. IDs without match become unassigned. @p map
    receives @p proteins and the data processing IDMapper records.
  */
  static void annotate(IDMapper& mapper, ConsensusMap& map, const FeatureMap& channels, const PeptideIdentificationList& ids,
                       const std::vector<ProteinIdentification>& proteins);

  /// Link only within a CV, preserving the input maps' channel subelements.
  static void group(FeatureGroupingAlgorithmQT& linker, const std::vector<ConsensusMap>& maps, ConsensusMap& result);

private:
  static CV getCV_(const MetaInfoInterface& value);
  static std::map<CV, std::vector<ConsensusMap>> split_(const std::vector<ConsensusMap>& maps);
  static void append_(ConsensusMap& result, ConsensusMap&& part);
};
} // namespace OpenMS
