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

#include <OpenMS/CONCEPT/Types.h>

#include <string>
#include <vector>
#include <map>

namespace OpenMS
{

  class SUPERHIRN_DLLAPI MS2Info
  {

    ////////////////////////////////////////////////
    // declaration of the private members:

private:

    int ID;

    std::string SQ;
    std::string FULL_SQ;
    std::string PREV_AA;
    std::vector<std::string> AC;
//  std::string ORIGINAL_INTERACT_FILE;
    std::string MS2_TYPE_TAG;

    // peptide prophet analysis:
    float PEP_PROB;

    // sorcerer search results:
    double XCORR;
    double DELTA_CN;

    double MONO_MZ;
    double THEO_MZ;
    double NEUTRAL_MR;

    int CHRG;
    int SCAN_START;
    int SCAN_END;

    double TR;

    std::map<int, double> MOD_LIST;

    // static values:
    static const double _MONO_H;
    static const double _MONO_O;

    ////////////////////////////////////////////////
    // declaration of the public members:

public:
    static const double mono_mass[26];
    static const char AA[20];
// static double MS2_TR_TOL;
// static bool THEO_MATCH_MODUS;
// static double MS2_MZ_PPM_TOLERANCE;

    // class destructor
    ~MS2Info();

    // class constructor
    MS2Info();
    MS2Info(int);
    MS2Info(std::string, std::string, float);
    MS2Info(std::string, std::string, int, float);
    MS2Info(std::string, std::string, float, int, int);
    // class copy constructor
    MS2Info(const MS2Info &);
    // class copy constructor
    MS2Info(const MS2Info *);

    //////////////////////////////////////////////////
    // overload operators:
    MS2Info & operator=(const MS2Info &);
    bool operator==(const MS2Info &);
    MS2Info & operator<=(const MS2Info &);
    MS2Info & operator>=(const MS2Info &);
    MS2Info & operator<(const MS2Info &);
    MS2Info & operator>(const MS2Info &);

    // add modification
    void add_modification(int, double);
    std::map<int, double>::iterator get_Modification_list_start();

    std::map<int, double>::iterator get_Modification_list_end();

    std::map<int, double>::iterator find_Modification(int pos);

    std::map<int, double> * get_Modification_list();

    bool check_MODIFICATION();

    // calculates the theoretical mass from a sequence:
    void set_THEO_MASS_from_SQ();
    double get_THEO_MZ();

    // sets modification SQ:
    void set_FULL_SQ();
    void set_SQ(std::string in);

    // show info:
    void show_info();

    // check whether proteotypic peptide:
    bool get_PROTEO_TYPE();
    // check the tryptic state:
    // 2: full tryptic
    // 1: semi tryptic
    // 0: non tryptic
    int get_TRYPTIC_STATE();

    // AC functions:
    // check if this AC or not:
    bool compare_AC(std::string);
    // search a pattern in the  AC list:
    bool search_AC_pattern(std::string);

    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class
    std::string get_SQ();

    std::string get_MOD_SQ();

    std::string get_TOTAL_SQ();

    std::string get_AC();

    std::vector<std::string> get_ALL_AC();

    std::vector<std::string>::iterator get_ALL_AC_START();

    std::vector<std::string>::iterator get_ALL_AC_END();

    bool find_AC(std::string);
    void set_AC(std::string);
    float get_PEP_PROB();

    void set_PEP_PROB(float in);

    double get_MONO_MZ();

    void set_MONO_MZ(double);

    double get_NEUTRAL_MR();

    void set_NEUTRAL_MR(double);

    int get_CHRG();

    void set_CHRG(int in);

    int get_SCAN();

    int get_SCAN_START();

    void set_SCAN_START(int in);

    int get_SCAN_END();

    void set_SCAN_END(int in);

    int get_ID();

    double get_DELTA_CN();

    void set_DELTA_CN(double in);

    double get_XCORR();

    void set_XCORR(double in);

    void set_MS2_TYPE_TAG(std::string in);

    std::string get_MS2_TYPE_TAG();

    // access the retention parameter:
    double getRetentionTime();

    void setRetentionTime(double in);

    double get_MONO_AA_MASS(int);

    std::string get_PREV_AA();

    void set_PREV_AA(std::string in);
  };

} // ns

