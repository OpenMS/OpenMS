// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// // --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $, Julianus Pfeuffer $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/TMTThirtyTwoPlexQuantitationMethod.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>

#include <algorithm>

namespace OpenMS{
    const String TMTThirtyTwoPlexQuantitationMethod::name_ = "tmt32plex";
    /*
        For 32plex experiments, use TMTpro 32plex Label Reagent Matched Set (Cat. No. A40000839) 
        or TMTpro 16plex Deuterated Label Reagent Set (Cat. No. A40000817) and TMTpro 16plex Label Reagent Set (Cat. No. A44520)
    */
    const std::vector<std::string> TMTThirtyTwoPlexQuantitationMethod::channel_names_ = {"126",
                                                                                        "127N","127C","127D",
                                                                                        "128N","128C","128ND","128CD",
                                                                                        "129N","129C","129ND","129CD",
                                                                                        "130N","130C","130ND","130CD",
                                                                                        "131N","131C","131ND","131CD",
                                                                                        "132N","132C","132ND","132CD",
                                                                                        "133N","133C","133ND","133CD",
                                                                                        "134N","134ND", "134CD",
                                                                                        "135ND"};

    static const int interaction_vector[32][14] = {  { -1, -1, -1, -1, -1, -1, -1,  1,  2,  3,  4,  6,  5,  7 },
                                        { -1, -1, -1, -1, -1, -1,  0, -1,  4,  6, -1, -1,  8, 10 },
                                        { -1, -1, -1, -1, -1,  0, -1,  4,  5,  7,  8, 10,  9, 11 },
                                        { -1, -1, -1, -1,  0, -1, -1,  6,  7, -1, 10,  9, 11, -1 },
                                        { -1, -1, -1,  0, -1,  1,  2, -1,  8, 10, -1, -1, 12, 14 },
                                        { -1,  0, -1, -1, -1,  2, -1,  8,  9, 11, 12, 14, 13, 15 },
                                        { -1, -1,  0, -1,  1, -1,  3, -1, 10,  9, -1, 12, 14, 13 },
                                        {  0, -1, -1, -1,  2,  3, -1, 10, 11, -1, 14, 13, 15, -1 },
                                        { -1,  1, -1,  2, -1,  4,  5, -1, 12, 14, -1, -1, 16, 18 },
                                        { -1,  2,  3, -1,  6,  5, -1, 12, 13, 15, 16, 18, 17, 19 },
                                        {  1, -1,  2,  3,  4,  6,  7, -1, 14, 13, -1, 16, 18, 17 },
                                        {  2,  3, -1, -1,  5,  7, -1, 14, 15, -1, 18, 17, 19, -1 },
                                        { -1,  4,  6,  5, -1,  8,  9, -1, 16, 18, -1, -1, 20, 22 },
                                        {  6,  5,  7, -1, 10,  9, -1, 16, 17, 19, 20, 22, 21, 23 },
                                        {  4,  6,  5,  7,  8, 10, 11, -1, 18, 17, -1, 20, 22, 21 },
                                        {  5,  7, -1, -1,  9, 11, -1, 18, 19, -1, 22, 21, 23, -1 },
                                        { -1,  8, 10,  9, -1, 12, 13, -1, 20, 22, -1, -1, 24, 26 },
                                        { 10,  9, 11, -1, 14, 13, -1, 20, 21, 23, 24, 26, 25, 27 },
                                        {  8, 10,  9, 11, 12, 14, 15, -1, 22, 21, -1, 24, 26, 25 },
                                        {  9, 11, -1, -1, 13, 15, -1, 22, 23, -1, 26, 25, 27, -1 },
                                        { -1, 12, 14, 13, -1, 16, 17, -1, 24, 26, -1, -1, 28, 29 },
                                        { 14, 13, 15, -1, 18, 17, -1, 24, 25, 27, 28, 29, -1, 30 },
                                        { 12, 14, 13, 15, 16, 18, 19, -1, 26, 25, -1, 28, 29, -1 },
                                        { 13, 15, -1, -1, 17, 19, -1, 26, 27, -1, 29, -1, 30, -1 },
                                        { -1, 16, 18, 17, -1, 20, 21, -1, 28, 29, -1, -1, -1, 31 },
                                        { 18, 17, 19, -1, 22, 21, -1, 28, -1, 30, -1, 31, -1, -1 },
                                        { 16, 18, 17, 19, 20, 22, 23, -1, 29, -1, -1, -1, 31, -1 },
                                        { 17, 19, -1, -1, 21, 23, -1, 29, 30, -1, 31, -1, -1, -1 },
                                        { -1, 20, 22, 21, -1, 24, 25, -1, -1, 31, -1, -1, -1, -1 },
                                        { 20, 22, 21, 23, 24, 26, 27, -1, 31, -1, -1, -1, -1, -1 },
                                        { 21, 23, -1, -1, 25, 27, -1, 31, -1, -1, -1, -1, -1, -1 },
                                        { 24, 26, 25, 27, 28, 29, 30, -1, -1, -1, -1, -1, -1, -1 }
    };
    double o_mass32[32] = {
        126.127726, 127.124761, 127.131081, 127.134003, 128.128116, 128.134436,
        128.131038, 128.137358, 129.131471, 129.137790, 129.134393, 129.140713,
        130.134825, 130.141145, 130.137748, 130.144068, 131.138180, 131.144500,
        131.141103, 131.147423, 132.141535, 132.147855, 132.144458, 132.150778,
        133.144890, 133.151210, 133.147813, 133.154133, 134.148245, 134.151171,
        134.157491, 135.154526
    };
    // auto make_channel = [&](const char* name, int idx, double mass, int vec_idx){
    //     return IsobaricChannelInformation(name, idx, "", mass, std::vector<int>(interaction_vector[vec_idx], interaction_vector[vec_idx]+14));
    // }
    TMTThirtyTwoPlexQuantitationMethod::TMTThirtyTwoPlexQuantitationMethod(){
        setName("TMTThirtyTwoPlexQuantitationMethod");

        // create the channel map
        // the following values for the reportor ion mass is taken from
        //  https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0018773_TMTproMassTagLabelingReagentsandKits_UG.pdf  
                                                                                // here the affected channels int the 
                                                                                // isobaric channel information as to updated
                                                                                // for the 16plex deuterated data.
                                                                                // 
        auto make_channel = [&](const std::string& name, int idx, double mass, int vec_idx){
            return IsobaricChannelInformation(name, idx, "", mass, std::vector<int>(interaction_vector[vec_idx], interaction_vector[vec_idx]+14));
        };

        for(int i =0; i < 32; i++){
            channels_.push_back(make_channel(channel_names_[i], i, o_mass32[i], i));
        };
        // we assume 126 to be the reference
        reference_channel_ = 0;

        setDefaultParams_();
    }

    void TMTThirtyTwoPlexQuantitationMethod::setDefaultParams_()
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
        defaults_.setValue("channel_134ND_description", "", "Description for the content of the 134ND channel.");
        defaults_.setValue("channel_134CD_description", "", "Description for the content of the 134CD channel.");
        defaults_.setValue("channel_135ND_description", "", "Description for the content of the 135ND channel.");

        defaults_.setValue("reference_channel", "126","The reference channel (126, 127N, 127C, 127D, 128N, 128C, 128ND, 128CD, 129N, 129C, 129ND, 129CD, 130N, 130C, 130ND, 130CD, 131N, 131C, 131ND, 131CD, 132N, 132C, 132ND, 132CD, 133N, 133C, 133ND, 133CD, 134N, 134ND, 134CD, 135ND).");
        defaults_.setValidStrings("reference_channel", TMTThirtyTwoPlexQuantitationMethod::channel_names_);

        // TODO :verify these below values
        defaults_.setValue("correction_matrix", std::vector<std::string>{
                                            "NA/NA/NA/NA            /NA/NA/NA          /0.31/9.09/NA      /0.02/NA/0.32/NA",//126 
                                            "NA/NA/NA/NA            /NA/NA/0.78        /NA/9.41/NA        /NA/NA/0.33/NA",//127N
                                            "NA/NA/NA/NA            /NA/0.93/NA        /0.35/8.63/NA      /0.01/NA/0.27/NA",//127C
                                            "NA/NA/NA/NA            /0.82/NA/NA        /0.30/8.71/0.33    /0.00/0.00/0.26/0.00",//127D
                                            "NA/NA/NA/0.00          /NA/0.82/0.65      /NA/8.13/NA        /NA/NA/0.26/NA",//128N
                                            "NA/0.00/NA/NA          /NA/1.47/NA        /0.34/6.91/NA      /0.00/NA/0.15/NA",//128C
                                            "NA/NA/0.00/NA          /1.13/NA/0.59     /NA/8.61/0.24      /NA/NA/0.28/0.00",//128ND
                                            "0.00/NA/NA/NA          /1.47/0.70/NA      /0.30/7.55/0.27    /0.00/0.00/0.21/0.00",//128CD
                                            "NA/0.00/NA/0.00        /NA/1.46/1.28      /NA/6.86/NA        /NA/NA/0.15/NA",//129N
                                            "NA/0.13/NA/NA          /NA/2.59/NA        /0.32/6.07/NA      /0.1/NA/0.09/NA",//129C
                                            "0.00/NA/0.00/0.00      /1.97/0.91/0.61    /NA/7.58/0.35      /NA/NA/0.19/0.00",//129ND
                                            "0.00/0.00/NA/NA        /1.39/1.37/NA      /0.31/6.02/0.50    /0.00/0.00/0.08/0.00",//129CD
                                            "NA/0.13/NA/0.00        /NA/2.41/0.27      /NA/5.58/NA        /NA/NA/0.10/NA",//130N
                                            "NA/0.04/NA/NA          /NA/3.10/NA        /0.42/4.82/NA      /0.02/NA/0.06/NA",//130C
                                            "0.00/0.08/0.00/0.00    /1.77/1.50/0.59    /NA/6.10/0.54      /NA/NA/0.09/0.00",//130ND
                                            "0.00/0.10/NA/NA        /1.44/2.18/NA      /0.33/5.49/0.40    /0.00/0.00/0.06/0.00",//130CD
                                            "NA/0.03/NA/0.00        /NA/2.78/0.63      /NA/4.57/NA        /NA/NA/0.12/NA",//131N
                                            "NA/0.08/NA/NA          /NA/3.90/NA        /0.47/3.57/NA      /0.00/NA/0.04/NA",//131C
                                            "0.01/0.19/0.00/0.00    /1.47/2.23/0.61    /NA/5.48/0.41      /NA/NA/0.05/0.00",//131ND
                                            "0.02/0.02/NA/NA        /2.52/2.61/NA      /0.98/4.14/1.15    /0.01/0.00/0.02/0.01",//131CD
                                            "NA/0.15/NA/0.01        /NA/3.58/0.72      /NA/1.80/NA        /NA/NA/0.00/NA",//132N
                                            "NA/0.11/NA/NA          /NA/4.55/NA        /0.43/1.86/NA      /0.00/NA/0.00/NA",//132C
                                            "0.02/0.02/0.00/0.00    /3.52/1.70/0.60    /NA/4.19/0.85      /NA/NA/0.01/0.00",//132ND
                                            "0.00/0.03/NA/NA        /1.00/3.35/NA      /1.00/3.07/0.94    /0.00/0.00/0.00/0.00",//132CD
                                            "NA/0.07/NA/0.01        /NA/3.14/0.73      /NA/3.40/NA        /NA/NA/0.03/NA",//133N
                                            "NA/0.22/NA/NA          /NA/4.96/NA        /0.34/1.03/NA      /0.00/NA/NA/NA",//133C
                                            "0.00/0.01/0.00/0.00    /0.89/2.48/0.63    /NA/3.14/0.51      /NA/NA/0.00/0.00",//133ND
                                            "0.01/0.01/NA/NA        /1.86/3.73/NA      /0.31/1.63/0.79    /0.00/0.00/0.00/0.00",//133CD
                                            "NA/0.30/NA/0.03        /NA/5.49/0.62      /NA/1.14/NA        /NA/NA/NA/NA", //134N
                                            "0.02/0.10/0.00/0.00    /1.66/3.76/0.22    /NA /1.73/0.31     /NA/NA/0.00/0.00",//134ND
                                            "0.03/0.23/NA/NA        /1.53/4.94/NA      /0.31/0.95/0.31    /0.00/0.00/NA/0.00",//134CD
                                            "0.04/0.23/0.00/0.00    /1.54/4.85/0.21    /NA/1.03/0.32      /NA/NA/NA/0.00",//135ND
                                            },
                                    "Correction matrix for isotope distributions in percent from the Thermo data sheet (see documentation);"
                                    "Please provide 32 entries (rows), separated by comma, where each entry contains 14 values in the following format: <-C13-H2>/<-2C13>/<-N15-H2>/<-C13-N15>/<-H2>/<-C13>/<-N15>/<+N15>/<+C13>/<+H2>/<+N15+C13>/<+N15+H2>/<+2C13>/<+C13+H2> e.g. one row may look like this: 'NA/NA/NA/NA   /0.82/NA/NA/  /0.30/8.71/0.33  /0.00/0.00/0.26/0.00'. You may use whitespaces at your leisure to ease reading.");
        defaultsToParam_();
    }

void TMTThirtyTwoPlexQuantitationMethod::updateMembers_(){
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
    channels_[29].description =  param_.getValue("channel_134ND_description").toString();
    channels_[30].description =  param_.getValue("channel_134CD_description").toString();
    channels_[31].description =  param_.getValue("channel_135ND_description").toString();

    std::vector<std::string>::const_iterator t_it = std::find(TMTThirtyTwoPlexQuantitationMethod::channel_names_.begin(),
                                                            TMTThirtyTwoPlexQuantitationMethod::channel_names_.end(),
                                                            param_.getValue("reference_channel"));

    reference_channel_ = t_it - TMTThirtyTwoPlexQuantitationMethod::channel_names_.begin();
}

TMTThirtyTwoPlexQuantitationMethod::TMTThirtyTwoPlexQuantitationMethod(const TMTThirtyTwoPlexQuantitationMethod&) = default;


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


/// below is the code for gettig the iteraction vector and the same can be done for the tmt35
/*
std::vector<std::vector<int>> tmt32(32, std::vector<int>(14));// to store the interference
std::vector<std::vector<double>> after_mass32(35, std::vector<double>(14));
double o_mass32[32] = {
    126.127726, 127.124761, 127.131081, 127.134003, 128.128116, 128.134436,
    128.131038, 128.137358, 129.131471, 129.137790, 129.134393, 129.140713,
    130.134825, 130.141145, 130.137748, 130.144068, 131.138180, 131.144500,
    131.141103, 131.147423, 132.141535, 132.147855, 132.144458, 132.150778,
    133.144890, 133.151210, 133.147813, 133.154133, 134.148245, 134.151171,
    134.157491, 135.154526
};// original mass

std::vector<int[14]> tmt35(35); // same here
std::vector<double[14]> after_mass35;
double o_mass35[35] = {
    126.127726, 127.124761, 127.131081, 127.134003, 128.128116,
    128.134436, 128.131038, 128.137358, 129.131471, 129.137790,
    129.134393, 129.140713, 130.134825, 130.141145, 130.137748,
    130.144068, 131.138180, 131.144500, 131.141103, 131.147423,
    132.141535, 132.147855, 132.144458, 132.150778, 133.144890,
    133.151210, 133.147813, 133.154133, 134.148245, 134.154566,
    134.151171, 134.157491, 135.151601,135.154526 ,135.160846
};  
// <-C13-H2",<-2C13",<-N15-H2",<-C13-N15",<-H2",<-C13",<-N15",<+N15",<+C13",<+H2",<+N15+C13",<+N15+H2",<+2C13",<+C13+H2>
std::vector<std::string> channel_seq = {"-C13-H2","-2C13","-N15-H2","-C13-N15","-H2","-C13","-N15","+N15","+C13","+H2","+N15+C13","+N15+H2","+2C13","+C13+H2"};
double weights[14]; // we will calculate only once for the above string

const double c = Constants::C13C12_MASSDIFF_U;
const double n = Constants::N15N14_MASSDIFF_U;
const double h = Constants::H2H1_MASSDIFF_U;

double get_mass(std::string inp){
    double out = 0.0f;
    if(inp.find("N15") < inp.size()){
        out += n;
    }
    if(inp.find("H2") < inp.size()){
        out += h;
    }
    if(inp.find("C13") < inp.size()){
        int i = inp.find("C13");
        if(inp[i-1] == '2'){
            out += 2*c;
        }else{
            out += c;
        }
    }
    return out;
}

double diff_n[32*14];

int find_closest_index(double val, const double *o_mass32,int idx, double tolerance = 0.0015) {
    int best_idx = -1;
    double min_diff = 1e9; // large initial diff

    for (int i = 0; i < 32; ++i) {
        double diff = fabs(val - o_mass32[i]);
        if (diff < min_diff) {
            min_diff = diff;
            best_idx = i;
        }
    }

    // If the closest value is still too far away, treat it as invalid
    // std::cout << min_diff << " ";
    diff_n[idx] = min_diff;
    if (min_diff > tolerance) return -1;
    return best_idx;
}


// O(n^2 * no of channels) <-- bad
void set_interference(std::vector<std::vector<int>> &intf, std::vector<std::vector<double>>& mass){
    for(int i = 0; i < 32; i++){
        for(int j = 0; j < 14; j++){
            intf[i][j] = find_closest_index(mass[i][j], o_mass32, i*j );
        }
        // std::cout << std::endl;
    }
}

int main(){

    double min_n = 0.0f;
    for (size_t i = 0; i < 31; i++){
        min_n = std::min(min_n, std::abs(o_mass35[i]-o_mass35[i+1]));
    }
    std::cout << "mass min diff " << min_n << std::endl;
    

    int i = 0;

    for(auto val : channel_seq){
        if(val.find("-") == 0){
            // std::cout << " k\n";
            weights[i] = -1.0f * get_mass(val);
        }else{
            // std::cout << "nah\n";
            weights[i] = get_mass(val);
        }
    i++;
    }
    // for 32plex
    for(int i = 0; i < 32; i++){
        for(int j = 0; j < 14; j++){
            after_mass32[i][j] = o_mass32[i] + weights[j]; 
        }
    }

    set_interference(tmt32, after_mass32);
    for(int i = 0; i<32;i++){
        std::cout << "{";
        for(int j = 0; j<14; j++){
            std::cout << tmt32[i][j] << ", ";
        }
        std::cout << "}\n";
    }
    min_n = 1.0f;
    for(int i = 0; i< 32*14; i++){
        min_n = std::max(min_n, diff_n[i]);
    }
    std::cout << "min is " << min_n << std::endl;
}

*/