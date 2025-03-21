// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------


#pragma once

#include <vector>
#include <OpenMS/config.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

namespace OpenMS
{
  class TransformationDescription;
  class ConsensusMap;
  class PeptideIdentification;
  class ConsensusFeature;

  /**
   * @brief This class collects functions for applying retention time transformations to data structures.
   */
  class OPENMS_DLLAPI MapAlignmentTransformer
  {

  public:
    /// Applies the given transformation to a peak map
    static void transformRetentionTimes(PeakMap& msexp,
                                        const TransformationDescription& trafo,
                                        bool store_original_rt = false);

    /// Applies the given transformation to a feature map
    static void transformRetentionTimes(FeatureMap& fmap,
                                        const TransformationDescription& trafo,
                                        bool store_original_rt = false);

    /// Applies the given transformation to a consensus map
    static void transformRetentionTimes(ConsensusMap& cmap,
                                        const TransformationDescription& trafo,
                                        bool store_original_rt = false);

    /// Applies the given transformation to peptide identifications
    static void transformRetentionTimes(
      std::vector<PeptideIdentification>& pep_ids,
      const TransformationDescription& trafo, bool store_original_rt = false);

    /// Applies the given transformation to input items in IdentificationData
    static void transformRetentionTimes(IdentificationData& id_data,
                                        const TransformationDescription& trafo,
                                        bool store_original_rt = false);

  private:
    /// Applies a transformation to a feature
    static void applyToFeature_(Feature& feature,
                                const TransformationDescription& trafo,
                                bool store_original_rt = false);

    /// Applies a transformation to a basic feature
    static void applyToBaseFeature_(BaseFeature& feature,
                                    const TransformationDescription& trafo,
                                    bool store_original_rt = false);

    /// Applies a transformation to a consensus feature
    static void applyToConsensusFeature_(
      ConsensusFeature& feature, const TransformationDescription& trafo,
      bool store_original_rt = false);

    /**
       @brief Stores the original RT in a meta value

       The original RT is written to a meta value "original_RT", but only if that meta value doesn't already exist.

       @returns True if the meta value was written, false if it already exists.
    */
    static bool storeOriginalRT_(MetaInfoInterface& meta_info,
                                 double original_rt);

  };
} // namespace OpenMS


