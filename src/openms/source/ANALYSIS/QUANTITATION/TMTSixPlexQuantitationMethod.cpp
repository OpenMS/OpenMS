// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

namespace OpenMS
{
  const String TMTSixPlexQuantitationMethod::name_ = "tmt6plex";

  TMTSixPlexQuantitationMethod::TMTSixPlexQuantitationMethod()
  {
    setName("TMTSixPlexQuantitationMethod");

    // create the channel map
    channels_.push_back(IsobaricChannelInformation("126", 0, "", 126.127725, {-1, -1, 1, 2}));
    channels_.push_back(IsobaricChannelInformation("127", 1, "", 127.124760, {-1, 0, 2, 3}));
    channels_.push_back(IsobaricChannelInformation("128", 2, "", 128.134433, {0, 1, 3, 4}));
    channels_.push_back(IsobaricChannelInformation("129", 3, "", 129.131468, {1, 2, 4, 5}));
    channels_.push_back(IsobaricChannelInformation("130", 4, "", 130.141141, {2, 3, 5, -1}));
    channels_.push_back(IsobaricChannelInformation("131", 5, "", 131.138176, {3, 4, -1, -1}));

    // we assume 126 to be the reference
    reference_channel_ = 0;

    setDefaultParams_();
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

    // default: Product Number: 90061 Lot Number: ZE386964
    defaults_.setValue("correction_matrix", std::vector<std::string>{
      "0.0/0.0/8.6/0.3",
      "0.0/0.1/7.8/0.1",
      "0.0/1.5/6.2/0.2",
      "0.0/1.5/5.7/0.1",
      "0.0/3.1/3.6/0.0",
      "0.1/2.9/3.8/0.0"
      },
      "Correction matrix for isotope distributions (see documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da>; e.g. '0/0.3/4/0', '0.1/0.3/3/0.2'");

    defaultsToParam_();
  }

  void TMTSixPlexQuantitationMethod::updateMembers_()
  {
    channels_[0].description = param_.getValue("channel_126_description").toString();
    channels_[1].description = param_.getValue("channel_127_description").toString();
    channels_[2].description = param_.getValue("channel_128_description").toString();
    channels_[3].description = param_.getValue("channel_129_description").toString();
    channels_[4].description = param_.getValue("channel_130_description").toString();
    channels_[5].description = param_.getValue("channel_131_description").toString();

    // compute the index of the reference channel
    reference_channel_ = ((Int) param_.getValue("reference_channel")) - 126;
  }

  TMTSixPlexQuantitationMethod::TMTSixPlexQuantitationMethod(const TMTSixPlexQuantitationMethod& other):
  IsobaricQuantitationMethod(other)
  {
    channels_.clear();
    channels_.insert(channels_.begin(), other.channels_.begin(), other.channels_.end());

    reference_channel_ = other.reference_channel_;
  }

  TMTSixPlexQuantitationMethod& TMTSixPlexQuantitationMethod::operator=(const TMTSixPlexQuantitationMethod& rhs) = default;

  const String& TMTSixPlexQuantitationMethod::getMethodName() const
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
    StringList iso_correction = ListUtils::toStringList<std::string>(getParameters().getValue("correction_matrix"));
    return stringListToIsotopeCorrectionMatrix_(iso_correction);
  }

  Size TMTSixPlexQuantitationMethod::getReferenceChannel() const
  {
    return reference_channel_;
  }

} // namespace
