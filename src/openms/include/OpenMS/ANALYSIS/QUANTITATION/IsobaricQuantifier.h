// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_ISOBARICQUANTIFIER_H
#define OPENMS_ANALYSIS_QUANTITATION_ISOBARICQUANTIFIER_H

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

#endif // OPENMS_ANALYSIS_QUANTITATION_ISOBARICQUANTIFIER_H
