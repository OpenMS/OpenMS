// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>

#include <OpenMS/DATASTRUCTURES/Matrix.h>

namespace OpenMS
{
  const String ItraqFourPlexQuantitationMethod::name_ = "itraq4plex";

  ItraqFourPlexQuantitationMethod::ItraqFourPlexQuantitationMethod()
  {
    setName("ItraqFourPlexQuantitationMethod");

    // create the channel map
    channels_.push_back(IsobaricChannelInformation("114", 0, "", 114.1112, -1, -1, 1, 2));
    channels_.push_back(IsobaricChannelInformation("115", 1, "", 115.1082, -1, 0, 2, 3));
    channels_.push_back(IsobaricChannelInformation("116", 2, "", 116.1116, 0, 1, 3, -1));
    channels_.push_back(IsobaricChannelInformation("117", 3, "", 117.1149, 1, 2, -1, -1));

    // we assume 114 to be the reference
    reference_channel_ = 0;

    setDefaultParams_();
  }

  ItraqFourPlexQuantitationMethod::~ItraqFourPlexQuantitationMethod()
  {
  }

  void ItraqFourPlexQuantitationMethod::setDefaultParams_()
  {
    defaults_.setValue("channel_114_description", "", "Description for the content of the 114 channel.");
    defaults_.setValue("channel_115_description", "", "Description for the content of the 115 channel.");
    defaults_.setValue("channel_116_description", "", "Description for the content of the 116 channel.");
    defaults_.setValue("channel_117_description", "", "Description for the content of the 117 channel.");
    defaults_.setValue("reference_channel", 114, "Number of the reference channel (114-117).");
    defaults_.setMinInt("reference_channel", 114);
    defaults_.setMaxInt("reference_channel", 117);

//    {0.0, 1.0, 5.9, 0.2},   //114
//    {0.0, 2.0, 5.6, 0.1},
//    {0.0, 3.0, 4.5, 0.1},
//    {0.1, 4.0, 3.5, 0.1}    //117
    defaults_.setValue("correction_matrix", ListUtils::create<String>("0.0/1.0/5.9/0.2,"
                                                               "0.0/2.0/5.6/0.1,"
                                                               "0.0/3.0/4.5/0.1,"
                                                               "0.1/4.0/3.5/0.1"),
                       "Correction matrix for isotope distributions (see documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da>; e.g. '0/0.3/4/0', '0.1/0.3/3/0.2'");

    defaultsToParam_();
  }

  void ItraqFourPlexQuantitationMethod::updateMembers_()
  {
    channels_[0].description = param_.getValue("channel_114_description");
    channels_[1].description = param_.getValue("channel_115_description");
    channels_[2].description = param_.getValue("channel_116_description");
    channels_[3].description = param_.getValue("channel_117_description");

    // compute the index of the reference channel
    reference_channel_ = ((Int) param_.getValue("reference_channel")) - 114;
  }

  ItraqFourPlexQuantitationMethod::ItraqFourPlexQuantitationMethod(const ItraqFourPlexQuantitationMethod& other):
  IsobaricQuantitationMethod(other)
  {
    channels_.clear();
    channels_.insert(channels_.begin(), other.channels_.begin(), other.channels_.end());

    reference_channel_ = other.reference_channel_;
  }

  ItraqFourPlexQuantitationMethod& ItraqFourPlexQuantitationMethod::operator=(const ItraqFourPlexQuantitationMethod& rhs)
  {
    if (this == &rhs)
      return *this;

    channels_.clear();
    channels_.insert(channels_.begin(), rhs.channels_.begin(), rhs.channels_.end());

    reference_channel_ = rhs.reference_channel_;

    return *this;
  }

  const String& ItraqFourPlexQuantitationMethod::getMethodName() const
  {
    return ItraqFourPlexQuantitationMethod::name_;
  }

  const IsobaricQuantitationMethod::IsobaricChannelList& ItraqFourPlexQuantitationMethod::getChannelInformation() const
  {
    return channels_;
  }

  Size ItraqFourPlexQuantitationMethod::getNumberOfChannels() const
  {
    return 4;
  }

  Matrix<double> ItraqFourPlexQuantitationMethod::getIsotopeCorrectionMatrix() const
  {
    StringList iso_correction = getParameters().getValue("correction_matrix");
    return stringListToIsotopCorrectionMatrix_(iso_correction);
  }

  Size ItraqFourPlexQuantitationMethod::getReferenceChannel() const
  {
    return reference_channel_;
  }

} // namespace
