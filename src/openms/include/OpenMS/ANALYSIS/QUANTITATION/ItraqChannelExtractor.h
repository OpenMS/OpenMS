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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_ITRAQCHANNELEXTRACTOR_H
#define OPENMS_ANALYSIS_QUANTITATION_ITRAQCHANNELEXTRACTOR_H

#include <vector>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  
  class ConsensusMap;

  /**
    @brief [experimental class] extracts the iTRAQ channels from tandem MS data and stores intensity values in a consensus map

    [experimental class]
    This class supports 4 and 8 channel iTRAQ and 6 channel TMT and will optionally do peak picking
    before the quantitation step. Quantitation is done by adding all signals within a small delta
    around the expected m/z of each channel. When all channels are found to be empty, the
    ConsensusFeature is not created. No post-processing is done here. Use ItraqQuantifier for that!

    @htmlinclude OpenMS_ItraqChannelExtractor.parameters
  */
  class OPENMS_DLLAPI ItraqChannelExtractor :
    public DefaultParamHandler,
    public ItraqConstants
  {

public:

    typedef ItraqConstants::ChannelMapType ChannelMapType;

    /// Constructor (assuming 4plex)
    ItraqChannelExtractor();

    /// Constructor with iTRAQ type (from enum ItraqConstants::ITRAQ_TYPES)
    explicit ItraqChannelExtractor(Int itraq_type);

    /// Constructor with iTRAQ type (from enum ItraqConstants::ITRAQ_TYPES) and param
    ItraqChannelExtractor(Int itraq_type, const Param & param);

    /// copy constructor
    ItraqChannelExtractor(const ItraqChannelExtractor & cp);

    /// assignment operator
    ItraqChannelExtractor & operator=(const ItraqChannelExtractor & rhs);


    /// @brief extracts the iTRAQ channels from the tandem MS data and stores intensity values in a consensus map
    ///
    /// @param ms_exp_data Raw data to read
    /// @param consensus_map
    void run(const MSExperiment<Peak1D> & ms_exp_data, ConsensusMap & consensus_map);

protected:

    void setDefaultParams_();

    /// implemented for DefaultParamHandler
    void updateMembers_();

private:

    /// initialize
    void init_();

    /// set to either ItraqConstants::FOURPLEX, ItraqConstants::EIGHTPLEX, or ItraqConstants::TMT_SIXPLEX
    Int itraq_type_;

    /// map the channel-name (e.g. 114) onto its description and the centroid mass
    /// the channel-name is also the id-string in the mapList section of the ConsensusMap
    ChannelMapType channel_map_;

  };   // !class

} // !namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_ITRAQCHANNELEXTRACTOR_H
