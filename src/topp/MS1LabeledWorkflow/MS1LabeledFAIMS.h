// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>
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

  /// Map IDs only onto multiplets acquired at the same CV.
  static void annotate(IDMapper& mapper, ConsensusMap& map, const PeptideIdentificationList& ids, const std::vector<ProteinIdentification>& proteins);

  /// Link only within a CV, preserving the input maps' channel subelements.
  static void group(FeatureGroupingAlgorithmQT& linker, const std::vector<ConsensusMap>& maps, ConsensusMap& result);

private:
  static CV getCV_(const MetaInfoInterface& value);
  static std::map<CV, std::vector<ConsensusMap>> split_(const std::vector<ConsensusMap>& maps);
  static void append_(ConsensusMap& result, ConsensusMap&& part);
};
} // namespace OpenMS
