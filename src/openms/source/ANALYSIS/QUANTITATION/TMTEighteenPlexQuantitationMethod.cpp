// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: Stephan Aiche, Radu Suciu $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/TMTEighteenPlexQuantitationMethod.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>


#include <algorithm>


namespace OpenMS
{
const String TMTEighteenPlexQuantitationMethod::name_ = "tmt18plex";
const std::vector<std::string> TMTEighteenPlexQuantitationMethod::channel_names_ = {"126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C","132N","132C","133N","133C","134N","134C","135N"};

TMTEighteenPlexQuantitationMethod::TMTEighteenPlexQuantitationMethod()
{
    setName("TMTEighteenPlexQuantitationMethod");

    //    // mass map outline - for further details please see #2427 (was adapted for tmt18plex)
    // create the channel map                                                //-2  -1  +1  +2
    channels_.push_back(IsobaricChannelInformation("126",   0, "", 126.127726, -1, -1,  2,  4));
    channels_.push_back(IsobaricChannelInformation("127N",  1, "", 127.124761, -1, -1,  3,  5));
    channels_.push_back(IsobaricChannelInformation("127C",  2, "", 127.131081, -1,  0,  4,  6));
    channels_.push_back(IsobaricChannelInformation("128N",  3, "", 128.128116, -1,  1,  5,  7));
    channels_.push_back(IsobaricChannelInformation("128C",  4, "", 128.134436,  -1,  2,  6,  8));
    channels_.push_back(IsobaricChannelInformation("129N",  5, "", 129.131471,  -1,  3,  7,  9)));
    channels_.push_back(IsobaricChannelInformation("129C",  6, "", 129.137790,  2,  4,  8, 10));
    channels_.push_back(IsobaricChannelInformation("130N",  7, "", 130.134825,  3,  5,  9, 11));
    channels_.push_back(IsobaricChannelInformation("130C",  8, "", 130.141145,  4,  6, 10, 12));
    channels_.push_back(IsobaricChannelInformation("131N",  9, "", 131.138180,  5,  7, 11, 13));
    channels_.push_back(IsobaricChannelInformation("131C", 10, "", 131.144500,  6,  8, 12, 14));
    channels_.push_back(IsobaricChannelInformation("132N", 11, "", 132.141535,  7,  9, 13, 15));
    channels_.push_back(IsobaricChannelInformation("132C", 12, "", 132.147855,  8,  10, 14, 16));
    channels_.push_back(IsobaricChannelInformation("133N", 13, "", 133.144890,  9,  11, 15, 17));
    channels_.push_back(IsobaricChannelInformation("133C", 14, "", 133.151210,  10,  12, 16, -1));
    channels_.push_back(IsobaricChannelInformation("134N", 15, "", 134.148245,  11,  13, 17, -1));
    channels_.push_back(IsobaricChannelInformation("134C", 16, "", 134.154565,  12,  14, -1, -1));
    channels_.push_back(IsobaricChannelInformation("135N", 17, "", 135.151600,  13,  15, -1, -1));

    // we assume 126 to be the reference
    reference_channel_ = 0;

    setDefaultParams_();
}

void TMTEighteenPlexQuantitationMethod::setDefaultParams_()
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
    defaults_.setValue("channel_131N_description", "", "Description for the content of the 131N channel.");
    defaults_.setValue("channel_131C_description", "", "Description for the content of the 131C channel.");
    defaults_.setValue("channel_132N_description", "", "Description for the content of the 132N channel.");
    defaults_.setValue("channel_132C_description", "", "Description for the content of the 132C channel.");
    defaults_.setValue("channel_133N_description", "", "Description for the content of the 133N channel.");
    defaults_.setValue("channel_133C_description", "", "Description for the content of the 133C channel.");
    defaults_.setValue("channel_134N_description", "", "Description for the content of the 134N channel.");
    defaults_.setValue("channel_134C_description", "", "Description for the content of the 134C channel.");
    defaults_.setValue("channel_135N_description", "", "Description for the content of the 135N channel.");

    defaults_.setValue("reference_channel", "126", "The reference channel (126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N, 134C, 135N).");
    defaults_.setValidStrings("reference_channel", TMTEighteenPlexQuantitationMethod::channel_names_);

    // TODO: verify these
    defaults_.setValue("correction_matrix", std::vector<std::string>{"0.0/0.0/9.40/0.34",
                                                                      "0.0/0.78/9.41/0.33",
                                                                      "0.0/0.93/8.98/0.28",
                                                                      "0.0/1.47/8.13/0.26",
                                                                      "0.0/1.47/7.25/0.15",
                                                                      "0.0/2.74/6.86/0.15",
                                                                      "0.13/2.59/6.39/0.10",
                                                                      "0.13/2.68/5.58/0.10",
                                                                      "0.04/3.10/5.24/0.08",
                                                                      "0.03/3.41/4.57/0.12",
                                                                      "0.08/3.90/4.04/0.04",
                                                                      "0.16/4.30/1.80/0.0",
                                                                      "0.11/4.55/2.29/0.0",
                                                                      "0.08/3.87/3.49/0.03",
                                                                      "0.22/4.96/1.37/0.0",
                                                                      "0.33/6.11/1.14/0.0",
                                                                      "0.14/5.81/0.31/0.0",
                                                                      "0.21/5.42/0.00/0.0"},
                       "Correction matrix for isotope distributions in percent (see documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da>; e.g. '0/0.3/4/0', '0.1/0.3/3/0.2'");

    defaultsToParam_();
}

void TMTEighteenPlexQuantitationMethod::updateMembers_()
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
    channels_[9].description = param_.getValue("channel_131N_description").toString();
    channels_[10].description = param_.getValue("channel_131C_description").toString();
    channels_[11].description = param_.getValue("channel_132N_description").toString();
    channels_[12].description = param_.getValue("channel_132C_description").toString();
    channels_[13].description = param_.getValue("channel_133N_description").toString();
    channels_[14].description = param_.getValue("channel_133C_description").toString();
    channels_[15].description = param_.getValue("channel_134N_description").toString();
    channels_[16].description = param_.getValue("channel_134C_description").toString();
    channels_[17].description = param_.getValue("channel_135N_description").toString();

    // compute the index of the reference channel
    std::vector<std::string>::const_iterator t_it = std::find(TMTEighteenPlexQuantitationMethod::channel_names_.begin(),
                                                         TMTEighteenPlexQuantitationMethod::channel_names_.end(),
                                                         param_.getValue("reference_channel"));

    reference_channel_ = t_it - TMTEighteenPlexQuantitationMethod::channel_names_.begin();
}

TMTEighteenPlexQuantitationMethod::TMTEighteenPlexQuantitationMethod(const TMTEighteenPlexQuantitationMethod& other):
IsobaricQuantitationMethod(other)
{
    channels_.clear();
    channels_.insert(channels_.begin(), other.channels_.begin(), other.channels_.end());

    reference_channel_ = other.reference_channel_;
}

TMTEighteenPlexQuantitationMethod& TMTEighteenPlexQuantitationMethod::operator=(const TMTEighteenPlexQuantitationMethod& rhs)
{
    if (this == &rhs)
    {
        return *this;
    }
    channels_.clear();
    channels_.insert(channels_.begin(), rhs.channels_.begin(), rhs.channels_.end());

    reference_channel_ = rhs.reference_channel_;

    return *this;
}

const String& TMTEighteenPlexQuantitationMethod::getMethodName() const
{
    return TMTEighteenPlexQuantitationMethod::name_;
}

const IsobaricQuantitationMethod::IsobaricChannelList& TMTEighteenPlexQuantitationMethod::getChannelInformation() const
{
    return channels_;
}

Size TMTEighteenPlexQuantitationMethod::getNumberOfChannels() const
{
    return 18;
}

Matrix<double> TMTEighteenPlexQuantitationMethod::getIsotopeCorrectionMatrix() const
{
    StringList iso_correction = ListUtils::toStringList<std::string>(getParameters().getValue("correction_matrix"));
    return stringListToIsotopCorrectionMatrix_(iso_correction);
}

Size TMTEighteenPlexQuantitationMethod::getReferenceChannel() const
{
    return reference_channel_;
}

} // namespace
