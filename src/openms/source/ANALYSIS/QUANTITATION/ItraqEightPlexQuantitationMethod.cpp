// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  const String ItraqEightPlexQuantitationMethod::name_ = "itraq8plex";

  ItraqEightPlexQuantitationMethod::ItraqEightPlexQuantitationMethod()
  {
    setName("ItraqFourPlexQuantitationMethod");

    // create the channel map
    channels_.push_back(IsobaricChannelInformation("113", 0, "", 113.1078, {-1, -1, 1, 2}));
    channels_.push_back(IsobaricChannelInformation("114", 1, "", 114.1112, {-1, 0, 2, 3}));
    channels_.push_back(IsobaricChannelInformation("115", 2, "", 115.1082, {0, 1, 3, 4}));
    channels_.push_back(IsobaricChannelInformation("116", 3, "", 116.1116, {1, 2, 4, 5}));
    channels_.push_back(IsobaricChannelInformation("117", 4, "", 117.1149, {2, 3, 5, 6}));
    channels_.push_back(IsobaricChannelInformation("118", 5, "", 118.1120, {3, 4, 6, 7}));
    channels_.push_back(IsobaricChannelInformation("119", 6, "", 119.1153, {4, 5, -1, 7}));
    channels_.push_back(IsobaricChannelInformation("121", 7, "", 121.1220, {6, -1, -1, -1}));

    // we assume 114 to be the reference
    reference_channel_ = 0;

    setDefaultParams_();
  }

  ItraqEightPlexQuantitationMethod::~ItraqEightPlexQuantitationMethod() = default;

  ItraqEightPlexQuantitationMethod::ItraqEightPlexQuantitationMethod(const ItraqEightPlexQuantitationMethod& other):
  IsobaricQuantitationMethod(other)
  {
    channels_.clear();
    channels_.insert(channels_.begin(), other.channels_.begin(), other.channels_.end());

    reference_channel_ = other.reference_channel_;
  }

  ItraqEightPlexQuantitationMethod& ItraqEightPlexQuantitationMethod::operator=(const ItraqEightPlexQuantitationMethod& rhs)
  {
    if (this == &rhs)
      return *this;

    channels_.clear();
    channels_.insert(channels_.begin(), rhs.channels_.begin(), rhs.channels_.end());

    reference_channel_ = rhs.reference_channel_;

    return *this;
  }

  void ItraqEightPlexQuantitationMethod::setDefaultParams_()
  {
    defaults_.setValue("channel_113_description", "", "Description for the content of the 113 channel.");
    defaults_.setValue("channel_114_description", "", "Description for the content of the 114 channel.");
    defaults_.setValue("channel_115_description", "", "Description for the content of the 115 channel.");
    defaults_.setValue("channel_116_description", "", "Description for the content of the 116 channel.");
    defaults_.setValue("channel_117_description", "", "Description for the content of the 117 channel.");
    defaults_.setValue("channel_118_description", "", "Description for the content of the 118 channel.");
    defaults_.setValue("channel_119_description", "", "Description for the content of the 119 channel.");
    defaults_.setValue("channel_121_description", "", "Description for the content of the 121 channel.");
    defaults_.setValue("reference_channel", 113, "Number of the reference channel (113-121). Please note that 120 is not valid.");
    defaults_.setMinInt("reference_channel", 113);
    defaults_.setMaxInt("reference_channel", 121);

//    {0.00, 0.00, 6.89, 0.22},   //113
//    {0.00, 0.94, 5.90, 0.16},
//    {0.00, 1.88, 4.90, 0.10},
//    {0.00, 2.82, 3.90, 0.07},
//    {0.06, 3.77, 2.99, 0.00},
//    {0.09, 4.71, 1.88, 0.00},
//    {0.14, 5.66, 0.87, 0.00},
//    {0.27, 7.44, 0.18, 0.00}    //121
    defaults_.setValue("correction_matrix", std::vector<std::string>{"0.00/0.00/6.89/0.22", //113
                                                                               "0.00/0.94/5.90/0.16",
                                                                               "0.00/1.88/4.90/0.10",
                                                                               "0.00/2.82/3.90/0.07",
                                                                               "0.06/3.77/2.99/0.00",
                                                                               "0.09/4.71/1.88/0.00",
                                                                               "0.14/5.66/0.87/0.00",
                                                                               "0.27/7.44/0.18/0.00"}, //121
                       "Correction matrix for isotope distributions (see documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da>; e.g. '0/0.3/4/0', '0.1/0.3/3/0.2'");

    defaultsToParam_();
  }

  void ItraqEightPlexQuantitationMethod::updateMembers_()
  {
    channels_[0].description = param_.getValue("channel_113_description").toString();
    channels_[1].description = param_.getValue("channel_114_description").toString();
    channels_[2].description = param_.getValue("channel_115_description").toString();
    channels_[3].description = param_.getValue("channel_116_description").toString();
    channels_[4].description = param_.getValue("channel_117_description").toString();
    channels_[5].description = param_.getValue("channel_118_description").toString();
    channels_[6].description = param_.getValue("channel_119_description").toString();
    channels_[7].description = param_.getValue("channel_121_description").toString();

    // compute the index of the reference channel
    Int ref_ch = param_.getValue("reference_channel");
    if (ref_ch == 121)
    {
      reference_channel_ = 7;
    }
    else if (ref_ch == 120)
    {
      OPENMS_LOG_WARN << "Invalid channel selection." << std::endl;
    }
    else
    {
      reference_channel_ = ref_ch - 113;
    }
  }

  const String& ItraqEightPlexQuantitationMethod::getMethodName() const
  {
    return ItraqEightPlexQuantitationMethod::name_;
  }

  const IsobaricQuantitationMethod::IsobaricChannelList& ItraqEightPlexQuantitationMethod::getChannelInformation() const
  {
    return channels_;
  }

  Size ItraqEightPlexQuantitationMethod::getNumberOfChannels() const
  {
    return 8;
  }

  Matrix<double> ItraqEightPlexQuantitationMethod::getIsotopeCorrectionMatrix() const
  {
    StringList iso_correction = ListUtils::toStringList<std::string>(getParameters().getValue("correction_matrix"));
    return stringListToIsotopeCorrectionMatrix_(iso_correction);
  }

  Size ItraqEightPlexQuantitationMethod::getReferenceChannel() const
  {
    return reference_channel_;
  }

} // namespace
