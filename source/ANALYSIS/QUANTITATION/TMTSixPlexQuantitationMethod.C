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
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>

namespace OpenMS
{
  const String TMTSixPlexQuantitationMethod::name_ = "tmt6plex";

  TMTSixPlexQuantitationMethod::TMTSixPlexQuantitationMethod()
  {
    setName("TMTSixPlexQuantitationMethod");

    // create the channel map
    channels_.push_back(IsobaricChannelInformation(126, 0, "", 126.127725));
    channels_.push_back(IsobaricChannelInformation(127, 1, "", 127.124760));
    channels_.push_back(IsobaricChannelInformation(128, 2, "", 128.134433));
    channels_.push_back(IsobaricChannelInformation(129, 3, "", 129.131468));
    channels_.push_back(IsobaricChannelInformation(130, 4, "", 130.141141));
    channels_.push_back(IsobaricChannelInformation(131, 5, "", 131.138176));

    // we assume 126 to be the reference
    reference_channel_ = 0;

    setDefaultParams_();
  }

  TMTSixPlexQuantitationMethod::~TMTSixPlexQuantitationMethod()
  {
  }

  void TMTSixPlexQuantitationMethod::setDefaultParams_()
  {
    defaults_.setValue("channel_126_description", "", "Description for the content of the 126 channel.");
    defaults_.setValue("channel_127_description", "", "Description for the content of the 127 channel.");
    defaults_.setValue("channel_128_description", "", "Description for the content of the 128 channel.");
    defaults_.setValue("channel_129_description", "", "Description for the content of the 129 channel.");
    defaults_.setValue("channel_130_description", "", "Description for the content of the 130 channel.");
    defaults_.setValue("channel_131_description", "", "Description for the content of the 131 channel.");
    defaults_.setValue("reference_channel", 126, "Number of the reference channel (126-131).");
    defaults_.setMinInt("reference_channel", 126);
    defaults_.setMaxInt("reference_channel", 131);

    //    {0.0, 1.0, 5.9, 0.2},   //114
    //    {0.0, 2.0, 5.6, 0.1},
    //    {0.0, 3.0, 4.5, 0.1},
    //    {0.1, 4.0, 3.5, 0.1}    //117
    defaults_.setValue("correction_matrix", StringList::create("0.0/0.0/0.0/0.0,"
                                                               "0.0/0.0/0.0/0.0,"
                                                               "0.0/0.0/0.0/0.0,"
                                                               "0.0/0.0/0.0/0.0,"
                                                               "0.0/0.0/0.0/0.0,"
                                                               "0.0/0.0/0.0/0.0"),
                       "Override default values (see Documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '0/0.3/4/0' , '0.1/0.3/3/0.2'");

    defaultsToParam_();
  }

  void TMTSixPlexQuantitationMethod::updateMembers_()
  {
    channels_[0].description = param_.getValue("channel_126_description");
    channels_[1].description = param_.getValue("channel_127_description");
    channels_[2].description = param_.getValue("channel_128_description");
    channels_[3].description = param_.getValue("channel_129_description");
    channels_[4].description = param_.getValue("channel_130_description");
    channels_[5].description = param_.getValue("channel_131_description");

    // compute the index of the reference channel
    reference_channel_ = ((Int) param_.getValue("reference_channel")) - 126;
  }

  TMTSixPlexQuantitationMethod::TMTSixPlexQuantitationMethod(const TMTSixPlexQuantitationMethod& other)
  {
    channels_.clear();
    channels_.insert(channels_.begin(), other.channels_.begin(), other.channels_.end());

    reference_channel_ = other.reference_channel_;
  }

  TMTSixPlexQuantitationMethod& TMTSixPlexQuantitationMethod::operator=(const TMTSixPlexQuantitationMethod& rhs)
  {
    if (this == &rhs)
      return *this;

    channels_.clear();
    channels_.insert(channels_.begin(), rhs.channels_.begin(), rhs.channels_.end());

    reference_channel_ = rhs.reference_channel_;

    return *this;
  }

  const String& TMTSixPlexQuantitationMethod::getName() const
  {
    return TMTSixPlexQuantitationMethod::name_;
  }

  const IsobaricQuantitationMethod::IsobaricChannelList& TMTSixPlexQuantitationMethod::getChannelInformation() const
  {
    return channels_;
  }

  Size TMTSixPlexQuantitationMethod::getNumberOfChannels() const
  {
    return 6;
  }

  Matrix<double> TMTSixPlexQuantitationMethod::getIsotopeCorrectionMatrix() const
  {
    StringList iso_correction = getParameters().getValue("correction_matrix");
    return stringListToIsotopCorrectionMatrix_(iso_correction);
  }

  Size TMTSixPlexQuantitationMethod::getReferenceChannel() const
  {
    return reference_channel_;
  }

} // namespace
