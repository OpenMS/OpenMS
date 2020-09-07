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
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//

#pragma once

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnConfig.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Info.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FeatureLCProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Feature.h>

#include <map>
#include <vector>
#include <string>

namespace OpenMS
{

  class SUPERHIRN_DLLAPI SHFeature
  {

    ////////////////////////////////////////////////
    // declaration of the private members:

private:

    /////////////////////////////////////////////
    // IDENTIFICATION PARAMETERS:
    // name of the spectra:
    std::map<double, std::vector<MS2Info> > MS2_SCANS;
    /////////////////////////////////////////////

    /////////////////////////////////////////////
    // raw MS peak parameters:
    int scan_apex;
    int scan_start;
    int scan_end;
    double total_peak_area;
    double apex_peak_intensity;
    double PEAK_SCORE;
    double SignalToNoise;
    double BackgroundNoise;

    ////////////////////////////////////////////
    // Analysis parameters:
    double alignment_error_up;
    double alignment_error_down;
    double SCORE_HOLDER;
    bool feature_match_status;
    double PI;

    ////////////////////////////////////////////
    // LC/MS run ID parameters:
    int spectrum_ID;
    int MASTER_ID;

    ///////////////////////////////////////////
    // string to store ms1 feature extra information:
    std::string featureExtraInformation;

    ///////////////////////////////////////////
    // LC elution profile:
    FeatureLCProfile * LCprofile;

    // static values:
    static double _MONO_H;
    static double _MONO_O;

    //////////////////////////////////////////
    // LC/MS matching things:
    std::map<int, SHFeature> matched_feature_list;

    // ranges of m/z and tr:
    double TR_APEX;
    double MONO_MZ_START;
    double MONO_MZ_END;
    double MONO_MZ_ORIGINAL;

    //////////////////////////////////////////
    // associated MS2 feature:
    MS2Feature * MS2TraceFeature;

    ////////////////////////////////////////////////
    // declaration of the public members:

public:

    double TR;
    double MONO_MZ;
    double TR_START;
    double TR_END;
    int charge_state;
    int feature_ID;

    // class destructor
    ~SHFeature();

    // copy constructor:
    SHFeature(const SHFeature &);

    // class constructor
    SHFeature(const SHFeature *);
    // constructor for the object SHFeature:
    SHFeature(double, double, int, int, int, int, float, float, float);
    // constructor for the object SHFeature:
    SHFeature(float, int, int);
    SHFeature(MS2Feature *);
    SHFeature();
    // copy constructor:
    SHFeature & operator=(const SHFeature &);

    // show the content of the spectra
    void show_info();
    // show MS/MS spectra info:
    void showMS2consensSpectraInfo();

    //////////////////////////////////
    // comparison operators:
    bool operator==(const SHFeature &);

    // add MS/MS info to the SHFeature:
    void add_MS2_info(MS2Info *);
    void add_MS2_info(std::map<double, std::vector<MS2Info> > *);
    bool get_MS2_info();
    bool get_MS2_info(double);

    bool check_MS2_empty();

    void removeAllMS2Information();

    int get_MS2_SCANS_SIZE();

    std::map<double, std::vector<MS2Info> > * get_MS2_SCAN_MAP();

    std::map<double, std::vector<MS2Info> >::iterator get_MS2_SCANS_START();
    std::map<double, std::vector<MS2Info> >::iterator get_MS2_SCANS_END();
    // get the best ms2 scan == closest to the apex:
    MS2Info * get_best_MS2_SCAN();
    MS2Info * get_best_MS2_SCAN(double);

    void setFeatureExtraInformation(std::string in);
    std::string getFeatureExtraInformation();

    // functions to set/access matched features:
    void add_matched_feature(SHFeature *);
    std::map<int, SHFeature> * get_match_list_REFERENCE();
    std::map<int, SHFeature> get_match_list();
    std::map<int, SHFeature>::iterator get_match_list_start();
    std::map<int, SHFeature>::iterator get_match_list_end();
    std::map<int, SHFeature>::iterator find_match_by_id(int ID);

    // get feature at a certain LC-MS by LC_MS id
    SHFeature * get_feature(int);

    // get the total peak are over all matched features:
    double get_MATCHED_peak_area();
    bool check_match_by_id(int);
    void erase_match_list();
    // get the profile over all matched features:
    std::map<int, double> get_feature_profile();

    // return number of times this feature has been seen = nb_replicates in list plus 1!
    int get_replicate_match_nb();
    int get_matching_nb();
    // return the sum of all intensities over replicates:
    double get_replicate_intensity_sum();

    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class

    // access the parent mass of feature, calculated from the SQ
    double get_MZ();
    void set_MZ(double in);
    double get_MZ_START();
    void set_MZ_START(double in);
    double get_MZ_END();
    void set_MZ_END(double in);

    double get_THEO_MZ();
    double get_THEO_MZ(double T);
    std::string get_AC();
    std::string get_AC(double T);
    bool check_AC(std::string in);
    bool check_AC(std::string in, double T);
    std::string get_SQ();
    std::string get_SQ(double T);
    std::string get_TOTAL_SQ();
    std::string get_TOTAL_SQ(double T);
    std::string get_MOD_SQ();
    std::string get_MOD_SQ(double T);
    double get_pep_prob();
    double get_pep_prob(double T);
    std::string get_MS2_TYPE_TAG();
    std::string get_MS2_TYPE_TAG(double T);
    int get_MS2_scan();
    int get_MS2_scan(double T);
    std::map<double, std::vector<MS2Info> > * get_MS2_SCAN_LIST();
    std::map<double, std::vector<MS2Info> >::iterator get_MS2_SCAN_LIST_START();
    std::map<double, std::vector<MS2Info> >::iterator get_MS2_SCAN_LIST_END();

    int get_scan_number();
    void set_scan_number(int in);
    int get_scan_start();
    void set_scan_start(int in);
    int get_scan_end();
    void set_scan_end(int in);
    int get_charge_state();
    void set_charge_state(int in);
    void set_peak_area(float in);
    double get_peak_area();
    // get peak area at a certain LC/MS:
    double get_peak_area(int);
    double get_apex_peak_intensity();
    void set_apex_peak_intensity(double in);
    void normalize_peak_area_by_factor(double factor);

    double get_alignment_error_up();
    void set_alignment_error_up(double in);
    double get_alignment_error_down();
    void set_alignment_error_down(double in);

    void set_SCORE_HOLDER(double in);
    double get_SCORE_HOLDER();

    double get_retention_time();
    void set_retention_time(double in);
    double get_retention_time_START();
    void set_retention_time_START(double in);
    double get_retention_time_END();
    void set_retention_time_END(double in);

    // original mz and Tr coordinates
    double get_raw_retention_time_apex();
    void set_raw_retention_time_apex(double in);
    double get_raw_MZ();
    void set_raw_MZ(double in);

    // feature ID:
    void set_feature_ID(int in);
    int get_feature_ID();

    void set_spectrum_ID(int in);

    int get_spectrum_ID();

    void set_MASTER_ID(int in);
    int get_MASTER_ID();

    // check how many matches
    int get_nb_common_match();

    // get/set the peak score
    double get_peak_score();
    void set_peak_score(double in);

    // get the molecular mass of the corresponding peptide!
    double get_Molecular_Mass();

    // feature PI:
    double get_FEATURE_PI();
    void set_FEATURE_PI(double in);

    // check charge states, in cases where a feature was
    // created based on a MS2 trace, charge state is unknown ( = -1 )
    // -> derive the charge state from the matched feature (if this is
    // also not -1
    void deriveChargeStates(SHFeature *);

    // LC elution profile
    void setLCelutionProfile(FeatureLCProfile * in);
    FeatureLCProfile * getLCelutionProfile();

    //////////////////////////////////////////////
    // parameters computed over matched features:
    double get_profile_retention_time();
    double get_profile_Molecular_Mass();

    /////////////////////////////////////////////
    // status if feature has been matched:
    bool get_feature_match_status();
    void set_feature_match_status(bool in);

    ///////////////////////////////////////////
    // access the MS2 feature
    void addMS2Feature(MS2Feature * in);
    void removeMS2Feature();
    MS2Feature * getMS2Feature();

    double getSignalToNoise();
    void setSignalToNoise(double in);

    double getBackgroundNoiseLevel();
    void setBackgroundNoiseLevel(double in);

    //////////////////////////////////////////////
    // get static members:
    //static double get_MZ_TOL(){return MZ_TOL;};
    static double get_MONO_H();

    // compare to masses at the PPM value and decided
    // if they fall into the m/z tolerance window
    static bool compareFeatureMassValuesAtPPMLevel(double, double);

    // get the masse error at the PPM value
    static double getFeatureMassErrorAtPPMLevel(double);

  };

} // ns

