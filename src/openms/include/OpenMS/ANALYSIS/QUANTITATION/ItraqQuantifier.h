// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_ITRAQQUANTIFIER_H
#define OPENMS_ANALYSIS_QUANTITATION_ITRAQQUANTIFIER_H

#include <vector>

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  class ConsensusMap;

  /**
    @brief Does post-processing on raw iTRAQ channel quantitation

    Using the raw consensus map from ItraqChannelExtractor, a non-negative isotope correction, normalization (using median)
    and [optionally] protein inference is computed.

    @htmlinclude OpenMS_ItraqQuantifier.parameters
  */
  class OPENMS_DLLAPI ItraqQuantifier :
    public DefaultParamHandler,
    public ItraqConstants
  {

public:

    typedef ItraqConstants::ChannelInfo ChannelInfo;
    typedef ItraqConstants::ChannelMapType ChannelMapType;
    typedef ItraqConstants::IsotopeMatrices IsotopeMatrices;

    /// Constructor (assuming 4-plex experiment)
    ItraqQuantifier();

    /// Constructor with iTRAQ-type (either ItraqConstants::FOURPLEX or ItraqConstants::EIGHTPLEX)
    explicit ItraqQuantifier(Int itraq_type);

    /// Constructor with iTRAQ-type (either ItraqConstants::FOURPLEX or ItraqConstants::EIGHTPLEX) and Param
    ItraqQuantifier(Int itraq_type, const Param & param);

    /// copy constructor
    ItraqQuantifier(const ItraqQuantifier & cp);

    /// assignment operator
    ItraqQuantifier & operator=(const ItraqQuantifier & rhs);

    /**
      @brief using the raw iTRAQ intensities we apply isotope correction, normalization (using median)

      @param consensus_map_in Raw iTRAQ intensities from previous step
      @param consensus_map_out Post-processed iTRAQ ratios for peptides

      @throws Exception::FailedAPICall is least-squares fit fails
      @throws Exception::InvalidParameter if parameter is invalid (e.g. reference_channel)
    */
    void run(const ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out);

    /**
      @brief Statistics for quantitation performance and comparison of NNLS vs. naive method (aka matrix inversion)
    */
    struct ItraqQuantifierStats
    {
      ItraqQuantifierStats() :
        channel_count(0),
        iso_number_ms2_negative(0),
        iso_number_reporter_negative(0),
        iso_number_reporter_different(0),
        iso_solution_different_intensity(0),
        iso_total_intensity_negative(0),
        number_ms2_total(0),
        number_ms2_empty(0),
        empty_channels()
      {
      }

      Size channel_count;  //< 4plex, 6plex, or 8 plex?!
      Size iso_number_ms2_negative; //< number of MS2 spectra where one or more channels had negative solution
      Size iso_number_reporter_negative;  //< number of channels where naive solution was negative
      Size iso_number_reporter_different; //< number of channels >0 where naive solution was different; happens when naive solution is negative in other channels
      double iso_solution_different_intensity; //< absolute intensity difference between both solutions (for channels > 0)
      double iso_total_intensity_negative; //< only for spectra where naive solution is negative
      Size number_ms2_total; //< total number of MS2 spectra
      Size number_ms2_empty; //< number of empty MS2 (no reporters at all)
      std::map<Size, Size> empty_channels; //< Channel_ID -> Missing; indicating the number of empty channels from all MS2 scans, i.e., numbers are between number_ms2_empty and number_ms2_total
    };

    ItraqQuantifierStats getStats() const;

protected:

    void setDefaultParams_();

    void updateMembers_();

private:

    /// initialize
    void initIsotopeCorrections_();

    void reconstructChannelInfo_(const ConsensusMap & consensus_map);

    /** 
      @brief Check if the given channel_frequency matrix is an identity matrix
     
      @param The matrix to check.
    */
    bool isIdentityCorrectionMatrix_(const Matrix<double>& channel_frequency) const;
    
    /// either ItraqConstants::FOURPLEX or ItraqConstants::EIGHTPLEX
    Int itraq_type_;

    /// map the channel-name (e.g., 114) onto its channel_info
    /// the channel-description is also the id-string in the mapList section of the ConsensusMap
    ChannelMapType channel_map_;

    /// Matrices with isotope correction values (one for each plex-type)
    IsotopeMatrices isotope_corrections_;

    /// stats for isotope correction
    ItraqQuantifierStats stats_;

  };   // !class

  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const ItraqQuantifier::ItraqQuantifierStats & stats);

} // !namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_ITRAQQUANTIFIER_H
