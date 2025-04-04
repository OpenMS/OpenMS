// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// // --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/TMTThirtyTwoPlexQuantitationMethod.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>

#include <algorithm>

namespace OpenMS{
    const String TMTThirtyTwoPlexQuantitationMethod::name_ = "tmt32plex";
    const std::vector<std::string> TMTThirtyTwoPlexQuantitationMethod::channel_names_ = {"127D", "128ND", "128CD", "129ND", "129CD", "130ND", "130CD", "131ND", "131CD", "132ND", "132CD", "133ND", "133CD", "134ND", "134CD", "135ND"}

    TMTThirtyTwoPlexQuantitationMethod::TMTThirtyTwoPlexQuantitationMethod()
    {
        setName("TMTThirtyTwoPlexQuantitationMethod");

        // create the channel map
        // todo VERIFY the 
       channels_.push_back(IsobaricChannelInformation("127D", 0, "", 127.124761, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("128ND", 1, "", 128.128116, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("128CD", 2, "", 128.134436, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("129ND", 3, "", 129.131471, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("129CD", 4, "", 129.137790, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("130ND", 5, "", 130.134825, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("130CD", 6, "", 130.141145, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("131ND", 7, "", 131.138180, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("131CD", 8, "", 131.144500, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("132ND", 9, "", 132.141535, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("132CD", 10, "", 132.147855, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("133ND", 11, "", 133.144890, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("133CD", 12, "", 133.151210, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("134ND", 13, "", 134.148245, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("134CD", 14, "", 134.154565, {-1, -1, -1, -1, -1, -1, -1, -1}));
       channels_.push_back(IsobaricChannelInformation("135ND", 15, "", 135.151600, {-1, -1, -1, -1, -1, -1, -1, -1}));

        // we assume 127D to be the reference
        reference_channel_ = 0;

        setDefaultParams_();
    }

    void TMTThirtyTwoPlexQuantitationMethod::setDefaultParams_()
    {
        defaults_.setValue("channel_127D_description", "", "Description for the content of the 127D channel.");
        defaults_.setValue("channel_128ND_description", "", "Description for the content of the 128ND channel.");
        defaults_.setValue("channel_128CD_description", "", "Description for the content of the 128CD channel.");
        defaults_.setValue("channel_129ND_description", "", "Description for the content of the 129ND channel.");
        defaults_.setValue("channel_129CD_description", "", "Description for the content of the 129CD channel.");
        defaults_.setValue("channel_130ND_description", "", "Description for the content of the 130ND channel.");
        defaults_.setValue("channel_130CD_description", "", "Description for the content of the 130CD channel.");
        defaults_.setValue("channel_131ND_description", "", "Description for the content of the 131ND channel.");
        defaults_.setValue("channel_131CD_description", "", "Description for the content of the 131CD channel.");
        defaults_.setValue("channel_132ND_description", "", "Description for the content of the 132ND channel.");
        defaults_.setValue("channel_132CD_description", "", "Description for the content of the 132CD channel.");
        defaults_.setValue("channel_133ND_description", "", "Description for the content of the 133ND channel.");
        defaults_.setValue("channel_133CD_description", "", "Description for the content of the 133CD channel.");
        defaults_.setValue("channel_134ND_description", "", "Description for the content of the 134ND channel.");
        defaults_.setValue("channel_134CD_description", "", "Description for the content of the 134CD channel.");
        defaults_.setValue("channel_135ND_description", "", "Description for the content of the 135ND channel.");

        default_.setValue("reference_channel", "127D","The reference channel (127D, 128ND, 128CD, 129ND, 129CD, 130ND, 130CD, 131ND, 131CD, 132ND, 132CD, 133ND, 133CD, 134ND, 134CD, 135ND).");
        defaults_.setValidStrings("reference_channel", TMTThirtyTwoPlexQuantitationMethod::channel_names_);

        // TODO :verify these below values
        defaults_.setValue("correction_matrix", std::vector<std::string>{
                                            "NA/NA/NA/NA            /0.82/NA/NA/       /0.30/8.71/0.33    /0.00/0.00/0.26/0.00",
                                            "NA/NA/0.00/NA          /1.13/NA/0.59/     /NA/8.61/0.24      /NA/NA/0.28/0.00",
                                            "0.00/NA/NA/NA          /1.47/0.70/NA      /0.30/7.55/0.27    /0.00/0.00/0.21/0.00",
                                            "0.00/NA/0.00/0.00      /1.97/0.91/0.61    /NA/7.58/0.35      /NA/NA/0.19/0.00",
                                            "0.00/0.00/NA/NA        /1.39/1.37/NA      /0.31/6.02/0.50    /0.00/0.00/0.08/0.00",
                                            "0.00/0.08/0.00/0.00    /1.77/1.50/0.59    /NA/6.10/0.54      /NA/NA/0.09/0.00",
                                            "0.00/0.10/NA/NA        /1.44/2.18/NA      /0.33/5.49/0.40    /0.00/0.00/0.06/0.00",
                                            "0.01/0.19/0.00/0.00    /1.47/2.23/0.61    /NA/5.48/0.41      /NA/NA/0.05/0.00",
                                            "0.02/0.02/NA/NA        /2.52/2.61/NA      /0.98/4.14/1.15    /0.01/0.00/0.02/0.01",
                                            "0.02/0.02/0.00/0.00    /3.52/1.70/0.60    /NA/4.19/0.85      /NA/NA/0.01/0.00",
                                            "0.00/0.03/NA/NA        /1.00/3.35/NA      /1.00/3.07/0.94    /0.00/0.00/0.00/0.00",
                                            "0.00/0.01/0.00/0.00    /0.89/2.48/0.63    /NA/3.14/0.51      /NA/NA/0.00/0.00",
                                            "0.01/0.01/NA/NA        /1.86/3.73/NA      /0.31/1.63/0.79    /0.00/0.00/0.00/0.00",
                                            "0.02/0.10/0.00/0.00    /1.66/3.76/0.22    /NA /1.73/0.31     /NA/NA/0.00/0.00",
                                            "0.03/0.23/NA/NA        /1.53/4.94/NA      /0.31/0.95/0.31    /0.00/0.00/NA/0.00",
                                            "0.04/0.23/0.00/0.00    /1.54/4.85/0.21    /NA/1.03/0.32      /NA/NA/NA/0.00",
                                            },
                                    "Correction matrix for isotope distributions in percent from the Thermo data sheet (see documentation);"
                                    "Please provide 16 entries (rows), separated by comma, where each entry contains 14 values in the following format: <-C13-H2>/<-2C13>/<-N15-H2>/<-C13-N15>/<-H2>/<-C13>/<-N15>/<+N15>/<+C13>/<+H2>/<+N15+C13>/<+N15+H2>/<+2C13>/<+C13+H2> e.g. one row may look like this: 'NA/NA/NA/NA   /0.82/NA/NA/  /0.30/8.71/0.33  /0.00/0.00/0.26/0.00'. You may use whitespaces at your leisure to ease reading.");
        defaultsToParam_();
    }

void TMTThirtyTwoPlexQuantitationMethod::updateMembers_(){
    channels_[0].description    =  param_.getValue("channel_127D_description").toString();
    channels_[1].description    =  param_.getValue("channel_128ND_description").toString();
    channels_[2].description    =  param_.getValue("channel_128CD_description").toString();
    channels_[3].description    =  param_.getValue("channel_129ND_description").toString();
    channels_[4].description    =  param_.getValue("channel_129CD_description").toString();
    channels_[5].description    =  param_.getValue("channel_130ND_description").toString();
    channels_[6].description    =  param_.getValue("channel_130CD_description").toString();
    channels_[7].description    =  param_.getValue("channel_131ND_description").toString();
    channels_[8].description    =  param_.getValue("channel_131CD_description").toString();
    channels_[9].description    =  param_.getValue("channel_132ND_description").toString();
    channels_[10].description   =  param_.getValue("channel_132CD_description").toString();
    channels_[11].description   =  param_.getValue("channel_133ND_description").toString();
    channels_[12].description   =  param_.getValue("channel_133CD_description").toString();
    channels_[13].description   =  param_.getValue("channel_134ND_description").toString();
    channels_[14].description   =  param_.getValue("channel_134CD_description").toString();
    channels_[15].description   =  param_.getValue("channel_135ND_description").toString();

    std::vector<std::string>::const_iterator t_it = std::find(TMTThirtyTwoPlexQuantitationMethod::channel_names_.begin(),
                                                            TMTThirtyTwoPlexQuantitationMethod::channel_names_.end(),
                                                            param_.getValue("reference_channel"));

    reference_channel_ = t_it - TMTThirtyTwoPlexQuantitationMethod::channel_names_.begin();
}

TMTThirtyTwoPlexQuantitationMethod::TMTThirtyTwoPlexQuantitationMethod(const TMTThirtyTwoPlexQuantitationMethod& other) :IsobaricQuantitationMethod(other){
    channels_.clear();
    channels_.insert(channels_.begin(), other.channels_.begin(), other.channels_.end());

    reference_channel_ = other.reference_channel_;
}

TMTThirtyTwoPlexQuantitationMethod& TMTThirtyTwoPlexQuantitationMethod::operator=(const TMTThirtyTwoPlexQuantitationMethod& rhs)
= default;

const String& TMTThirtyTwoPlexQuantitationMethod::getMethodName() const
{
    return TMTThirtyTwoPlexQuantitationMethod::name_;
}

const IsobaricQuantitationMethod::IsobaricChannelList& TMTThirtyTwoPlexQuantitationMethod::getChannelInformation() const
{
    return channels_;
}

Size TMTThirtyTwoPlexQuantitationMethod::getNumberOfChannels() const
{
    return channels_.size();
}

Matrix<double> TMTThirtyTwoPlexQuantitationMethod::getIsotopeCorrectionMatrix() const
{
    StringList iso_correction = ListUtils::toStringList<std::string>(getParameters().getValue("correction_matrix"));
    return stringListToIsotopeCorrectionMatrix_(iso_correction);
}

Size TMTThirtyTwoPlexQuantitationMethod::getReferenceChannel() const
{
    return reference_channel_;
}

} // namespace OpenMS