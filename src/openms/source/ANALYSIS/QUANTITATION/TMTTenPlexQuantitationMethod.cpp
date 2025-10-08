// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

namespace OpenMS
{
const String TMTTenPlexQuantitationMethod::name_ = "tmt10plex";
const std::vector<std::string> TMTTenPlexQuantitationMethod::channel_names_ = {"126","127N","127C","128N","128C","129N","129C","130N","130C","131"};

TMTTenPlexQuantitationMethod::TMTTenPlexQuantitationMethod()
{
    setName("TMTTenPlexQuantitationMethod");

    //    // mass map outline - for further details please see #2427
    //    "126", 126.127726, x, x, 127C, 128C
    //    "127N", 127.124761, x, x, 128N, 129N
    //    "127C", 127.131081, x, 126, 128C, 129C
    //    "128N", 128.128116, x, 127N, 129N, 130N
    //    "128C", 128.134436, 126, 127C, 129C, 130C
    //    "129N", 129.131471, 127N, 128N, 130N, 131
    //    "129C", 129.137790, 127C, 128C, 130C, x
    //    "130N", 130.134825, 128N, 129N, 131, x
    //    "130C", 130.141145, 128C, 129C, x, x
    //    "131", 131.138180, 129N, 130N, x, x

    // create the channel map                                               //-2  -1  +1  +2
    channels_.push_back(IsobaricChannelInformation("126",  0, "", 126.127726, {-1, -1,  2,  4}));
    channels_.push_back(IsobaricChannelInformation("127N", 1, "", 127.124761, {-1, -1,  3,  5}));
    channels_.push_back(IsobaricChannelInformation("127C", 2, "", 127.131081, {-1,  0,  4,  6}));
    channels_.push_back(IsobaricChannelInformation("128N", 3, "", 128.128116, {-1,  1,  5,  7}));
    channels_.push_back(IsobaricChannelInformation("128C", 4, "", 128.134436, {0,  2,  6,  8}));
    channels_.push_back(IsobaricChannelInformation("129N", 5, "", 129.131471, {1,  3,  7,  9}));
    channels_.push_back(IsobaricChannelInformation("129C", 6, "", 129.137790, {2,  4,  8, -1}));
    channels_.push_back(IsobaricChannelInformation("130N", 7, "", 130.134825, {3,  5,  9, -1}));
    channels_.push_back(IsobaricChannelInformation("130C", 8, "", 130.141145, {4,  6, -1, -1}));
    channels_.push_back(IsobaricChannelInformation("131",  9, "", 131.138180, {5,  7, -1, -1}));

    // we assume 126 to be the reference
    reference_channel_ = 0;

    setDefaultParams_();
}

void TMTTenPlexQuantitationMethod::setDefaultParams_()
{
    defaults_.setValue("channel_126_description", "", "Description for the content of the 126 channel.");
    defaults_.setValue("channel_127N_description", "", "Description for the content of the 127N channel.");
    defaults_.setValue("channel_127C_description", "", "Description for the content of the 127C channel.");
    defaults_.setValue("channel_128N_description", "", "Description for the content of the 128N channel.");
    defaults_.setValue("channel_128C_description", "", "Description for the content of the 128C channel.");
    defaults_.setValue("channel_129N_description", "", "Description for the content of the 129N channel.");
    defaults_.setValue("channel_129C_description", "", "Description for the content of the 129C channel.");
    defaults_.setValue("channel_130N_description", "", "Description for the content of the 130N channel.");
    defaults_.setValue("channel_130C_description", "", "Description for the content of the 130C channel.");
    defaults_.setValue("channel_131_description", "", "Description for the content of the 131 channel.");

    defaults_.setValue("reference_channel", "126", "The reference channel (126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131).");
    defaults_.setValidStrings("reference_channel", TMTTenPlexQuantitationMethod::channel_names_);

    defaults_.setValue("correction_matrix", std::vector<std::string>{"0.0/0.0/5.09/0.0",
                                                                      "0.0/0.25/5.27/0.0",
                                                                      "0.0/0.37/5.36/0.15",
                                                                      "0.0/0.65/4.17/0.1",
                                                                      "0.08/0.49/3.06/0.0",
                                                                      "0.01/0.71/3.07/0.0",
                                                                      "0.0/1.32/2.62/0.0",
                                                                      "0.02/1.28/2.75/2.53",
                                                                      "0.03/2.08/2.23/0.0",
                                                                      "0.08/1.99/1.65/0.0"},
                       "Correction matrix for isotope distributions (see documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da>; e.g. '0/0.3/4/0', '0.1/0.3/3/0.2'");

    defaultsToParam_();
}

void TMTTenPlexQuantitationMethod::updateMembers_()
{
    channels_[0].description = param_.getValue("channel_126_description").toString();
    channels_[1].description = param_.getValue("channel_127N_description").toString();
    channels_[2].description = param_.getValue("channel_127C_description").toString();
    channels_[3].description = param_.getValue("channel_128N_description").toString();
    channels_[4].description = param_.getValue("channel_128C_description").toString();
    channels_[5].description = param_.getValue("channel_129N_description").toString();
    channels_[6].description = param_.getValue("channel_129C_description").toString();
    channels_[7].description = param_.getValue("channel_130N_description").toString();
    channels_[8].description = param_.getValue("channel_130C_description").toString();
    channels_[9].description = param_.getValue("channel_131_description").toString();

    // compute the index of the reference channel
    std::vector<std::string>::const_iterator t_it = std::find(TMTTenPlexQuantitationMethod::channel_names_.begin(),
                                                         TMTTenPlexQuantitationMethod::channel_names_.end(),
                                                         param_.getValue("reference_channel"));

    reference_channel_ = t_it - TMTTenPlexQuantitationMethod::channel_names_.begin();
}

TMTTenPlexQuantitationMethod::TMTTenPlexQuantitationMethod(const TMTTenPlexQuantitationMethod& other):
IsobaricQuantitationMethod(other)
{
    channels_.clear();
    channels_.insert(channels_.begin(), other.channels_.begin(), other.channels_.end());

    reference_channel_ = other.reference_channel_;
}

TMTTenPlexQuantitationMethod& TMTTenPlexQuantitationMethod::operator=(const TMTTenPlexQuantitationMethod& rhs)
= default;

const String& TMTTenPlexQuantitationMethod::getMethodName() const
{
    return TMTTenPlexQuantitationMethod::name_;
}

const IsobaricQuantitationMethod::IsobaricChannelList& TMTTenPlexQuantitationMethod::getChannelInformation() const
{
    return channels_;
}

Size TMTTenPlexQuantitationMethod::getNumberOfChannels() const
{
    return 10;
}

Matrix<double> TMTTenPlexQuantitationMethod::getIsotopeCorrectionMatrix() const
{
    StringList iso_correction = ListUtils::toStringList<std::string>(getParameters().getValue("correction_matrix"));
    return stringListToIsotopeCorrectionMatrix_(iso_correction);
}

Size TMTTenPlexQuantitationMethod::getReferenceChannel() const
{
    return reference_channel_;
}

} // namespace
