// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/ANALYSIS/QUANTITATION/TMTElevenPlexQuantitationMethod.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>

#include <algorithm>

namespace OpenMS
{
const String TMTElevenPlexQuantitationMethod::name_ = "tmt11plex";
const std::vector<String> TMTElevenPlexQuantitationMethod::channel_names_ = {"126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C"};

TMTElevenPlexQuantitationMethod::TMTElevenPlexQuantitationMethod()
{
    setName("TMTElevenPlexQuantitationMethod");

    //    // mass map outline - for further details please see #2427 (was adapted for tmt11plex)
    //    "126", 126.127726, x, x, 127C, 128C
    //    "127N", 127.124761, x, x, 128N, 129N
    //    "127C", 127.131081, x, 126, 128C, 129C
    //    "128N", 128.128116, x, 127N, 129N, 130N
    //    "128C", 128.134436, 126, 127C, 129C, 130C
    //    "129N", 129.131471, 127N, 128N, 130N, 131N
    //    "129C", 129.137790, 127C, 128C, 130C, 131C
    //    "130N", 130.134825, 128N, 129N, 131N, x
    //    "130C", 130.141145, 128C, 129C, 131C, x
    //    "131N", 131.138180, 129N, 130N, x, x
    //    "131C", 131.144500, 129C, 130C, x, x

    // create the channel map                                                //-2  -1  +1  +2
    channels_.push_back(IsobaricChannelInformation("126",   0, "", 126.127726, -1, -1,  2,  4));
    channels_.push_back(IsobaricChannelInformation("127N",  1, "", 127.124761, -1, -1,  3,  5));
    channels_.push_back(IsobaricChannelInformation("127C",  2, "", 127.131081, -1,  0,  4,  6));
    channels_.push_back(IsobaricChannelInformation("128N",  3, "", 128.128116, -1,  1,  5,  7));
    channels_.push_back(IsobaricChannelInformation("128C",  4, "", 128.134436,  0,  2,  6,  8));
    channels_.push_back(IsobaricChannelInformation("129N",  5, "", 129.131471,  1,  3,  7,  9));
    channels_.push_back(IsobaricChannelInformation("129C",  6, "", 129.137790,  2,  4,  8, 10));
    channels_.push_back(IsobaricChannelInformation("130N",  7, "", 130.134825,  3,  5,  9, -1));
    channels_.push_back(IsobaricChannelInformation("130C",  8, "", 130.141145,  4,  6, 10, -1));
    channels_.push_back(IsobaricChannelInformation("131N",  9, "", 131.138180,  5,  7, -1, -1));
    channels_.push_back(IsobaricChannelInformation("131C", 10, "", 131.144500,  6,  8, -1, -1));


    // Original 10plex channel
    // channels_.push_back(IsobaricChannelInformation("131", 9, "", 131.138180, 5, 7, -1, -1));

    // we assume 126 to be the reference
    reference_channel_ = 0;

    setDefaultParams_();
}

void TMTElevenPlexQuantitationMethod::setDefaultParams_()
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

    defaults_.setValue("reference_channel", "126", "The reference channel (126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C).");
    defaults_.setValidStrings("reference_channel", TMTElevenPlexQuantitationMethod::channel_names_);

    defaults_.setValue("correction_matrix", ListUtils::create<String>("0.0/0.0/0.0/0.0,"
                                                                      "0.0/0.0/0.0/0.0,"
                                                                      "0.0/0.0/0.0/0.0,"
                                                                      "0.0/0.0/0.0/0.0,"
                                                                      "0.0/0.0/0.0/0.0,"
                                                                      "0.0/0.0/0.0/0.0,"
                                                                      "0.0/0.0/0.0/0.0,"
                                                                      "0.0/0.0/0.0/0.0,"
                                                                      "0.0/0.0/0.0/0.0,"
                                                                      "0.0/0.0/0.0/0.0,"
                                                                      "0.0/0.0/0.0/0.0"),
                       "Correction matrix for isotope distributions (see documentation); use the following format: <-2Da>/<-1Da>/<+1Da>/<+2Da>; e.g. '0/0.3/4/0', '0.1/0.3/3/0.2'");

    defaultsToParam_();
}

void TMTElevenPlexQuantitationMethod::updateMembers_()
{
    channels_[0].description = param_.getValue("channel_126_description");
    channels_[1].description = param_.getValue("channel_127N_description");
    channels_[2].description = param_.getValue("channel_127C_description");
    channels_[3].description = param_.getValue("channel_128N_description");
    channels_[4].description = param_.getValue("channel_128C_description");
    channels_[5].description = param_.getValue("channel_129N_description");
    channels_[6].description = param_.getValue("channel_129C_description");
    channels_[7].description = param_.getValue("channel_130N_description");
    channels_[8].description = param_.getValue("channel_130C_description");
    channels_[9].description = param_.getValue("channel_131N_description");
    channels_[10].description = param_.getValue("channel_131C_description");

    // compute the index of the reference channel
    std::vector<String>::const_iterator t_it = std::find(TMTElevenPlexQuantitationMethod::channel_names_.begin(),
                                                         TMTElevenPlexQuantitationMethod::channel_names_.end(),
                                                         (String) param_.getValue("reference_channel"));

    reference_channel_ = t_it - TMTElevenPlexQuantitationMethod::channel_names_.begin();
}

TMTElevenPlexQuantitationMethod::TMTElevenPlexQuantitationMethod(const TMTElevenPlexQuantitationMethod& other):
IsobaricQuantitationMethod(other)
{
    channels_.clear();
    channels_.insert(channels_.begin(), other.channels_.begin(), other.channels_.end());

    reference_channel_ = other.reference_channel_;
}

TMTElevenPlexQuantitationMethod& TMTElevenPlexQuantitationMethod::operator=(const TMTElevenPlexQuantitationMethod& rhs)
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

const String& TMTElevenPlexQuantitationMethod::getMethodName() const
{
    return TMTElevenPlexQuantitationMethod::name_;
}

const IsobaricQuantitationMethod::IsobaricChannelList& TMTElevenPlexQuantitationMethod::getChannelInformation() const
{
    return channels_;
}

Size TMTElevenPlexQuantitationMethod::getNumberOfChannels() const
{
    return 11;
}

Matrix<double> TMTElevenPlexQuantitationMethod::getIsotopeCorrectionMatrix() const
{
    StringList iso_correction = getParameters().getValue("correction_matrix");
    return stringListToIsotopCorrectionMatrix_(iso_correction);
}

Size TMTElevenPlexQuantitationMethod::getReferenceChannel() const
{
    return reference_channel_;
}

} // namespace
