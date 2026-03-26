// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Julianus Pfeuffer $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/TMTThirtyFivePlexQuantitationMethod.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <array>
#include <algorithm>

namespace OpenMS{
    const String TMTThirtyFivePlexQuantitationMethod::name_ = "tmt35plex";
    /*
        For 35plex experiments use TMTpro 16plex Deuterated Label Reagent Set (Cat. No. A40000817), 
        TMTpro-135CD Label Reagent (Cat. No. A40000818), and TMTpro 18plex Label Reagent Set (Cat. No. A52045).
    */
    const std::vector<std::string> TMTThirtyFivePlexQuantitationMethod::channel_names_ = {"126",
                                                                                        "127N","127C","127D",
                                                                                        "128N","128C","128ND","128CD",
                                                                                        "129N","129C","129ND","129CD",
                                                                                        "130N","130C","130ND","130CD",
                                                                                        "131N","131C","131ND","131CD",
                                                                                        "132N","132C","132ND","132CD",
                                                                                        "133N","133C","133ND","133CD",
                                                                                        "134N","134C","134ND","134CD",
                                                                                        "135N","135ND","135CD"};

    // Channel adjacency topology: maps each channel to its neighbors in mass space
    // for isotope correction. Each row has 14 entries matching the 14-column correction
    // matrix format. Values are channel indices (-1 = no neighbor at that offset).
    // With the default all-NA correction_matrix, no correction is applied regardless
    // of these values. Users supplying calibrated correction percentages need this
    // topology to route corrections to the correct target channels.
    static const std::array<std::array<int, 14>, 35> interaction_vector = {{
                                        {{ -1, -1, -1, -1, -1, -1, -1,  1,  2,  3,  4,  6,  5,  7 }},
                                        {{ -1, -1, -1, -1, -1, -1,  0, -1,  4,  6, -1, -1,  8, 10 }},
                                        {{ -1, -1, -1, -1, -1,  0, -1,  4,  5,  7,  8, 10,  9, 11 }},
                                        {{ -1, -1, -1, -1,  0, -1, -1,  6,  7, -1, 10,  9, 11, -1 }},
                                        {{ -1, -1, -1,  0, -1,  1,  2, -1,  8, 10, -1, -1, 12, 14 }},
                                        {{ -1,  0, -1, -1, -1,  2, -1,  8,  9, 11, 12, 14, 13, 15 }},
                                        {{ -1, -1,  0, -1,  1, -1,  3, -1, 10,  9, -1, 12, 14, 13 }},
                                        {{  0, -1, -1, -1,  2,  3, -1, 10, 11, -1, 14, 13, 15, -1 }},
                                        {{ -1,  1, -1,  2, -1,  4,  5, -1, 12, 14, -1, -1, 16, 18 }},
                                        {{ -1,  2,  3, -1,  6,  5, -1, 12, 13, 15, 16, 18, 17, 19 }},
                                        {{  1, -1,  2,  3,  4,  6,  7, -1, 14, 13, -1, 16, 18, 17 }},
                                        {{  2,  3, -1, -1,  5,  7, -1, 14, 15, -1, 18, 17, 19, -1 }},
                                        {{ -1,  4,  6,  5, -1,  8,  9, -1, 16, 18, -1, -1, 20, 22 }},
                                        {{  6,  5,  7, -1, 10,  9, -1, 16, 17, 19, 20, 22, 21, 23 }},
                                        {{  4,  6,  5,  7,  8, 10, 11, -1, 18, 17, -1, 20, 22, 21 }},
                                        {{  5,  7, -1, -1,  9, 11, -1, 18, 19, -1, 22, 21, 23, -1 }},
                                        {{ -1,  8, 10,  9, -1, 12, 13, -1, 20, 22, -1, -1, 24, 26 }},
                                        {{ 10,  9, 11, -1, 14, 13, -1, 20, 21, 23, 24, 26, 25, 27 }},
                                        {{  8, 10,  9, 11, 12, 14, 15, -1, 22, 21, -1, 24, 26, 25 }},
                                        {{  9, 11, -1, -1, 13, 15, -1, 22, 23, -1, 26, 25, 27, -1 }},
                                        {{ -1, 12, 14, 13, -1, 16, 17, -1, 24, 26, -1, -1, 28, 30 }},
                                        {{ 14, 13, 15, -1, 18, 17, -1, 24, 25, 27, 28, 30, 29, 31 }},
                                        {{ 12, 14, 13, 15, 16, 18, 19, -1, 26, 25, -1, 28, 30, 29 }},
                                        {{ 13, 15, -1, -1, 17, 19, -1, 26, 27, -1, 30, 29, 31, -1 }},
                                        {{ -1, 16, 18, 17, -1, 20, 21, -1, 28, 30, -1, -1, 32, 33 }},
                                        {{ 18, 17, 19, -1, 22, 21, -1, 28, 29, 31, 32, 33, -1, 34 }},
                                        {{ 16, 18, 17, 19, 20, 22, 23, -1, 30, 29, -1, 32, 33, -1 }},
                                        {{ 17, 19, -1, -1, 21, 23, -1, 30, 31, -1, 33, -1, 34, -1 }},
                                        {{ -1, 20, 22, 21, -1, 24, 25, -1, 32, 33, -1, -1, -1, -1 }},
                                        {{ 22, 21, 23, -1, 26, 25, -1, 32, -1, 34, -1, -1, -1, -1 }},
                                        {{ 20, 22, 21, 23, 24, 26, 27, -1, 33, -1, -1, -1, -1, -1 }},
                                        {{ 21, 23, -1, -1, 25, 27, -1, 33, 34, -1, -1, -1, -1, -1 }},
                                        {{ -1, 24, 26, 25, -1, 28, 29, -1, -1, -1, -1, -1, -1, -1 }},
                                        {{ 24, 26, 25, 27, 28, 30, 31, -1, -1, -1, -1, -1, -1, -1 }},
                                        {{ 25, 27, -1, -1, 29, 31, -1, -1, -1, -1, -1, -1, -1, -1 }}
                                    }};
    static const std::array<double, 35> o_mass35 = {
        126.127726, 127.124761, 127.131081, 127.134003, 128.128116,
        128.134436, 128.131038, 128.137358, 129.131471, 129.137790,
        129.134393, 129.140713, 130.134825, 130.141145, 130.137748,
        130.144068, 131.138180, 131.144500, 131.141103, 131.147423,
        132.141535, 132.147855, 132.144458, 132.150778, 133.144890,
        133.151210, 133.147813, 133.154133, 134.148245, 134.154566,
        134.151171, 134.157491, 135.151601,135.154526 ,135.160846
    };

    TMTThirtyFivePlexQuantitationMethod::TMTThirtyFivePlexQuantitationMethod(){
        setName("TMTThirtyFivePlexQuantitationMethod");

        // Reporter ion masses from Thermo TMTpro documentation:
        // https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0018773_TMTproMassTagLabelingReagentsandKits_UG.pdf
        auto make_channel = [&](const std::string& name, int idx, double mass, int vec_idx){
            return IsobaricChannelInformation(name, idx, "", mass, std::vector<int>(interaction_vector[vec_idx].begin(), interaction_vector[vec_idx].end()));
        };
        for(size_t i =0; i<o_mass35.size(); i++){
            channels_.push_back(make_channel(channel_names_[i],i,o_mass35[i],i));
        };                                                                        
        // we assume 126 to be the reference
        reference_channel_ = 0;
        setDefaultParams_();
    }

    void TMTThirtyFivePlexQuantitationMethod::setDefaultParams_()
    {
        defaults_.setValue("channel_126_description", "", "Description for the content of the 126 channel.");
        defaults_.setValue("channel_127N_description", "", "Description for the content of the 127N channel.");
        defaults_.setValue("channel_127C_description", "", "Description for the content of the 127C channel.");
        defaults_.setValue("channel_127D_description", "", "Description for the content of the 127D channel.");
        defaults_.setValue("channel_128N_description", "", "Description for the content of the 128N channel.");
        defaults_.setValue("channel_128C_description", "", "Description for the content of the 128C channel.");
        defaults_.setValue("channel_128ND_description", "", "Description for the content of the 128ND channel.");
        defaults_.setValue("channel_128CD_description", "", "Description for the content of the 128CD channel.");
        defaults_.setValue("channel_129N_description", "", "Description for the content of the 129N channel.");
        defaults_.setValue("channel_129C_description", "", "Description for the content of the 129C channel.");
        defaults_.setValue("channel_129ND_description", "", "Description for the content of the 129ND channel.");
        defaults_.setValue("channel_129CD_description", "", "Description for the content of the 129CD channel.");
        defaults_.setValue("channel_130N_description", "", "Description for the content of the 130N channel.");
        defaults_.setValue("channel_130C_description", "", "Description for the content of the 130C channel.");
        defaults_.setValue("channel_130ND_description", "", "Description for the content of the 130ND channel.");
        defaults_.setValue("channel_130CD_description", "", "Description for the content of the 130CD channel.");
        defaults_.setValue("channel_131N_description", "", "Description for the content of the 131N channel.");
        defaults_.setValue("channel_131C_description", "", "Description for the content of the 131C channel.");
        defaults_.setValue("channel_131ND_description", "", "Description for the content of the 131ND channel.");
        defaults_.setValue("channel_131CD_description", "", "Description for the content of the 131CD channel.");
        defaults_.setValue("channel_132N_description", "", "Description for the content of the 132N channel.");
        defaults_.setValue("channel_132C_description", "", "Description for the content of the 132C channel.");
        defaults_.setValue("channel_132ND_description", "", "Description for the content of the 132ND channel.");
        defaults_.setValue("channel_132CD_description", "", "Description for the content of the 132CD channel.");
        defaults_.setValue("channel_133N_description", "", "Description for the content of the 133N channel.");
        defaults_.setValue("channel_133C_description", "", "Description for the content of the 133C channel.");
        defaults_.setValue("channel_133ND_description", "", "Description for the content of the 133ND channel.");
        defaults_.setValue("channel_133CD_description", "", "Description for the content of the 133CD channel.");
        defaults_.setValue("channel_134N_description", "", "Description for the content of the 134N channel.");
        defaults_.setValue("channel_134C_description", "", "Description for the content of the 134C channel.");
        defaults_.setValue("channel_134ND_description", "", "Description for the content of the 134ND channel.");
        defaults_.setValue("channel_134CD_description", "", "Description for the content of the 134CD channel.");
        defaults_.setValue("channel_135N_description", "", "Description for the content of the 135N channel.");
        defaults_.setValue("channel_135ND_description", "", "Description for the content of the 135ND channel.");
        defaults_.setValue("channel_135CD_description", "", "Description for the content of the 135CD channel.");

        defaults_.setValue("reference_channel", "126","The reference channel (126, 127N, 127C, 127D, 128N, 128C, 128ND, 128CD, 129N, 129C, 129ND, 129CD, 130N, 130C, 130ND, 130CD, 131N, 131C, 131ND, 131CD, 132N, 132C, 132ND, 132CD, 133N, 133C, 133ND, 133CD, 134N, 134C, 134ND, 134CD, 135N, 135ND, 135CD).");
        defaults_.setValidStrings("reference_channel", TMTThirtyFivePlexQuantitationMethod::channel_names_);

        // Identity matrix (no isotope correction). Calibrated correction values are
        // not yet available for TMT 35-plex. The correction_matrix parameter is ignored;
        // getIsotopeCorrectionMatrix() always returns an identity matrix regardless of user input.
        //
        // Each row has 14 NA entries → all off-diagonal contributions are 0, diagonal is 1.0.
        const std::string identity_row = "NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA";
        defaults_.setValue("correction_matrix",
                           std::vector<std::string>(35, identity_row),
                           "Ignored for TMT 35-plex (identity matrix is always used). "
                           "Isotope correction is not yet supported for this method.");
        defaultsToParam_();
    }

void TMTThirtyFivePlexQuantitationMethod::updateMembers_(){
    channels_[0].description = param_.getValue("channel_126_description").toString();
    channels_[1].description = param_.getValue("channel_127N_description").toString();
    channels_[2].description = param_.getValue("channel_127C_description").toString();
    channels_[3].description =  param_.getValue("channel_127D_description").toString();
    channels_[4].description = param_.getValue("channel_128N_description").toString();
    channels_[5].description = param_.getValue("channel_128C_description").toString();
    channels_[6].description =  param_.getValue("channel_128ND_description").toString();
    channels_[7].description =  param_.getValue("channel_128CD_description").toString();
    channels_[8].description = param_.getValue("channel_129N_description").toString();
    channels_[9].description = param_.getValue("channel_129C_description").toString();
    channels_[10].description =  param_.getValue("channel_129ND_description").toString();
    channels_[11].description =  param_.getValue("channel_129CD_description").toString();
    channels_[12].description = param_.getValue("channel_130N_description").toString();
    channels_[13].description = param_.getValue("channel_130C_description").toString();
    channels_[14].description =  param_.getValue("channel_130ND_description").toString();
    channels_[15].description =  param_.getValue("channel_130CD_description").toString();
    channels_[16].description = param_.getValue("channel_131N_description").toString();
    channels_[17].description = param_.getValue("channel_131C_description").toString();
    channels_[18].description  =  param_.getValue("channel_131ND_description").toString();
    channels_[19].description  =  param_.getValue("channel_131CD_description").toString();
    channels_[20].description = param_.getValue("channel_132N_description").toString();
    channels_[21].description = param_.getValue("channel_132C_description").toString();
    channels_[22].description  =  param_.getValue("channel_132ND_description").toString();
    channels_[23].description =  param_.getValue("channel_132CD_description").toString();
    channels_[24].description = param_.getValue("channel_133N_description").toString();
    channels_[25].description = param_.getValue("channel_133C_description").toString();
    channels_[26].description =  param_.getValue("channel_133ND_description").toString();
    channels_[27].description =  param_.getValue("channel_133CD_description").toString();
    channels_[28].description = param_.getValue("channel_134N_description").toString();
    channels_[29].description = param_.getValue("channel_134C_description").toString();
    channels_[30].description =  param_.getValue("channel_134ND_description").toString();
    channels_[31].description =  param_.getValue("channel_134CD_description").toString();
    channels_[32].description = param_.getValue("channel_135N_description").toString();
    channels_[33].description =  param_.getValue("channel_135ND_description").toString();
    channels_[34].description =  param_.getValue("channel_135CD_description").toString();

    std::vector<std::string>::const_iterator t_it = std::find(TMTThirtyFivePlexQuantitationMethod::channel_names_.begin(),
                                                            TMTThirtyFivePlexQuantitationMethod::channel_names_.end(),
                                                            param_.getValue("reference_channel"));

    reference_channel_ = t_it - TMTThirtyFivePlexQuantitationMethod::channel_names_.begin();
}

TMTThirtyFivePlexQuantitationMethod::TMTThirtyFivePlexQuantitationMethod(const TMTThirtyFivePlexQuantitationMethod&) = default;


TMTThirtyFivePlexQuantitationMethod& TMTThirtyFivePlexQuantitationMethod::operator=(const TMTThirtyFivePlexQuantitationMethod& rhs)
= default;

const String& TMTThirtyFivePlexQuantitationMethod::getMethodName() const
{
    return TMTThirtyFivePlexQuantitationMethod::name_;
}

const IsobaricQuantitationMethod::IsobaricChannelList& TMTThirtyFivePlexQuantitationMethod::getChannelInformation() const
{
    return channels_;
}

Size TMTThirtyFivePlexQuantitationMethod::getNumberOfChannels() const
{
    return channels_.size();
}

Matrix<double> TMTThirtyFivePlexQuantitationMethod::getIsotopeCorrectionMatrix() const
{
    // Always return identity matrix — isotope correction is not yet validated for TMT 35-plex.
    // User-supplied correction_matrix parameter is ignored.
    const std::string identity_row = "NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA";
    StringList iso_correction(35, identity_row);
    return stringListToIsotopeCorrectionMatrix_(iso_correction);
}

Size TMTThirtyFivePlexQuantitationMethod::getReferenceChannel() const
{
    return reference_channel_;
}

} // namespace OpenMS
