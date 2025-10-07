// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>

namespace OpenMS
{
  class IsobaricQuantitationMethod;
  class ConsensusMap;

  /**
    @brief Given the extracted channel intensities the IsobaricQuantifier corrects and normalizes
           the intensities for further processing.

    @htmlinclude OpenMS_IsobaricQuantifier.parameters
  */
  class OPENMS_DLLAPI IsobaricQuantifier :
    public DefaultParamHandler
  {
public:
    /**
      @brief Constructor given an IsobaricQuantitationMethod (e.g., iTRAQ 4 plex).

      @param quant_method The quantification method used for the data set to analyze.
    */
    explicit IsobaricQuantifier(const IsobaricQuantitationMethod* const quant_method);

    /// Copy c'tor
    IsobaricQuantifier(const IsobaricQuantifier& other);

    /// Assignment operator
    IsobaricQuantifier& operator=(const IsobaricQuantifier& rhs);

    /**
      @brief Using the raw isobaric intensities we apply isotope correction, normalization (using median).

      @param consensus_map_in Raw isobaric channel intensities from channel extraction.
      @param consensus_map_out Corrected and normalized isobaric channel ratios for peptides.

      @throws Exception::FailedAPICall is least-squares fit fails
      @throws Exception::InvalidParameter if parameter is invalid (e.g. reference_channel)
    */
    void quantify(const ConsensusMap& consensus_map_in, ConsensusMap& consensus_map_out);

protected:
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

    /// implemented for DefaultParamHandler
    void updateMembers_() override;

private:
    /// Stats of current quantitation run.
    IsobaricQuantifierStatistics stats_;

    /// The quantification method used for the dataset to be analyzed.
    const IsobaricQuantitationMethod* quant_method_;

    /// Is true if isotope correction is enabled, false otherwise.
    bool isotope_correction_enabled_;

    /// Is true if normalization is enabled, false otherwise.
    bool normalization_enabled_;

    /// Computes labeling statistics (efficiency, number of empty scans,...)
    void computeLabelingStatistics_(ConsensusMap& consensus_map_out);
  };
} // namespace

