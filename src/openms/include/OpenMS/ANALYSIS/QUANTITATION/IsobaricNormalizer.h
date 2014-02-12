// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_ISOBARICNORMALIZER_H
#define OPENMS_ANALYSIS_QUANTITATION_ISOBARICNORMALIZER_H

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{
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
    Int reference_channel_name_;

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
    Map<Size, Size> map_to_vec_index_;

    /// The index of the reference channel in the peptide ratio/intensity vectors.
    Size ref_map_id_;

    /// Collection containing the collected peptide ratios for the individual channels.
    std::vector<std::vector<Peak2D::IntensityType> > peptide_ratios_;

    /// Collection containing the collected peptide intensities for the individual channels.
    std::vector<std::vector<Peak2D::IntensityType> > peptide_intensities_;

  };
} // namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_ISOBARICNORMALIZER_H
