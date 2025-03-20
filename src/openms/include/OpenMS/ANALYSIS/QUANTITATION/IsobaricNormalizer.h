// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>

#include <map>

namespace OpenMS
{
  class IsobaricQuantitationMethod;
  class ConsensusMap;

  /**
    @brief Performs median normalization on the extracted ratios of isobaric labeling experiment.
  */
  class OPENMS_DLLAPI IsobaricNormalizer
  {
public:
    /// Default c'tor
    explicit IsobaricNormalizer(const IsobaricQuantitationMethod* const quant_method);

    /// Copy c'tor
    IsobaricNormalizer(const IsobaricNormalizer& other);

    /// Assignment operator
    IsobaricNormalizer& operator=(const IsobaricNormalizer& rhs);

    /**
      @brief Normalizes the intensity ratios in the given input map (using median).

      @param consensus_map The map to normalize.
    */
    void normalize(ConsensusMap& consensus_map);

private:
    /// The selected quantitation method that will be used for the normalization.
    const IsobaricQuantitationMethod* quant_meth_;

    /// The name of the reference channel as given in the IsobaricChannelInformation.
    String reference_channel_name_;

    /**
      @brief Given a ConsensusFeature the method will return an iterator pointing to the consensus element representing the reference channel.

      @param cf The ConsensusFeature for which the reference element should be found.
      @param consensus_map The ConsensusMap in which the reference element should be found.
      @return An iterator pointing to the consensus element of the reference channel. ConsensusFeature::end() if the reference channel is not contained.
    */
    ConsensusFeature::HandleSetType::iterator findReferenceChannel_(ConsensusFeature& cf, const ConsensusMap& consensus_map) const;

    /**
      @brief Constructs a mapping from file description to the index in the corresponding ratio/intensity vectors.

      @param consensus_map The consensus map for which the mapping should be build.
    */
    void buildVectorIndex_(const ConsensusMap& consensus_map);

    /**
      @brief Collects ratios and intensities for a given ConsensusFeature.

      @param cf The consensus feature to evaluate.
      @param ref_intensity The intensity of the reference channel.
    */
    void collectRatios_(const ConsensusFeature& cf,
                        const Peak2D::IntensityType& ref_intensity);

    /**
      @brief Computes the normalization factors from the given peptide ratios.

      @param normalization_factors The normalization factors to compute.
    */
    void computeNormalizationFactors_(std::vector<Peak2D::IntensityType>& normalization_factors);

    /// The mapping between map indices and the corresponding indices in the peptide ratio/intensity vectors.
    std::map<Size, Size> map_to_vec_index_;

    /// The index of the reference channel in the peptide ratio/intensity vectors.
    Size ref_map_id_;

    /// Collection containing the collected peptide ratios for the individual channels.
    std::vector<std::vector<Peak2D::IntensityType> > peptide_ratios_;

    /// Collection containing the collected peptide intensities for the individual channels.
    std::vector<std::vector<Peak2D::IntensityType> > peptide_intensities_;

  };
} // namespace

