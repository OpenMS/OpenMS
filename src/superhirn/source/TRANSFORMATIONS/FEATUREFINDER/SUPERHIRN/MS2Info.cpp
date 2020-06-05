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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Info.h>

#include <string>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cstdio>

namespace OpenMS
{

  using namespace std;

// static values:
  const double MS2Info::_MONO_H = 1.00728;
  const double MS2Info::_MONO_O = 15.99491;

//  monoisotopic masses of all amino acids:
  const double MS2Info::mono_mass[26] = {
    71.03711, 0.0, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891,
    113.08406, 0.0, 128.09496, 113.08406, 131.04049, 114.04293, 0.0, 97.05276, 128.05858, 156.10111, 87.03203,
    101.04768, 0.0, 99.06841, 186.07931, 0.0, 163.06333, 0.0
  };
//  one letter code of all amino acids:
  const char MS2Info::AA[20] =
  {
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T',
    'V', 'W', 'Y'
  };

// double MS2Info::MS2_TR_TOL;
// bool MS2Info::THEO_MATCH_MODUS;
// double MS2Info::MS2_MZ_PPM_TOLERANCE;

////////////////////////////////////////////////
// constructor for the object MS2Info:
  MS2Info::MS2Info()
  {
    PEP_PROB = 0;
    DELTA_CN = 0;
    XCORR = 0;
    THEO_MZ = 0;
    MONO_MZ = 0;
    NEUTRAL_MR = 0;
    CHRG = 0;
    SCAN_START = 0;
    SCAN_END = 0;
    ID = -1;
    TR = -1.0;
  }

////////////////////////////////////////////////
// constructor for the object MS2Info:
  MS2Info::MS2Info(string IN_AC, string IN_SQ, float IN_PEP)
  {
    PEP_PROB = IN_PEP;
    DELTA_CN = 0;
    XCORR = 0;
    THEO_MZ = 0;
    MONO_MZ = 0;
    NEUTRAL_MR = 0;
    CHRG = 0;
    ID = -1;
    TR = -1.0;
    SQ = IN_SQ;
    set_AC(IN_AC);
    set_THEO_MASS_from_SQ();
    set_FULL_SQ();

  }

////////////////////////////////////////////////
// constructor for the object MS2Info:
  MS2Info::MS2Info(string IN_AC, string IN_SQ, int IN_CHRG, float IN_PEP)
  {
    PEP_PROB = IN_PEP;
    DELTA_CN = 0;
    XCORR = 0;
    THEO_MZ = 0;
    MONO_MZ = 0;
    NEUTRAL_MR = 0;
    ID = -1;
    TR = -1.0;
    SQ = IN_SQ;
    set_AC(IN_AC);
    CHRG = IN_CHRG;
    set_THEO_MASS_from_SQ();
    set_FULL_SQ();
  }

////////////////////////////////////////////////
// constructor for the object MS2Info:
  MS2Info::MS2Info(string IN_AC, string IN_SQ, float IN_PEP, int IN_CHRG, int IN_SCAN)
  {
    PEP_PROB = IN_PEP;
    THEO_MZ = 0;
    MONO_MZ = 0;
    DELTA_CN = 0;
    XCORR = 0;
    NEUTRAL_MR = 0;
    ID = -1;
    TR = -1.0;
    SQ = IN_SQ;
    set_AC(IN_AC);
    CHRG = IN_CHRG;
    SCAN_START = IN_SCAN;
    SCAN_END = IN_SCAN;
    set_THEO_MASS_from_SQ();
    set_FULL_SQ();

  }

////////////////////////////////////////////////
// constructor for the object MS2Info:
  MS2Info::MS2Info(int IN_ID)
  {
    ID = IN_ID;
    PEP_PROB = 0;
    THEO_MZ = 0;
    DELTA_CN = 0;
    XCORR = 0;
    MONO_MZ = 0;
    NEUTRAL_MR = 0;
    CHRG = 0;
    SCAN_START = 0;
    SCAN_END = 0;
    TR = -1.0;
  }

//////////////////////////////////////////////////
// class desctructor of MS2Info
  MS2Info::~MS2Info()
  {
    MOD_LIST.clear();
    FULL_SQ.clear();
    SQ.clear();
    AC.clear();
    TR = -1.0;
  }

//////////////////////////////////////////////////
// class copy constructor of MS2Info
  MS2Info::MS2Info(const MS2Info & tmp)
  {
    ID = tmp.ID;
    PEP_PROB = tmp.PEP_PROB;
    DELTA_CN = tmp.DELTA_CN;
    XCORR = tmp.XCORR;
    THEO_MZ = tmp.THEO_MZ;
    MONO_MZ = tmp.MONO_MZ;
    NEUTRAL_MR = tmp.NEUTRAL_MR;
    CHRG = tmp.CHRG;
    SCAN_START = tmp.SCAN_START;
    SCAN_END = tmp.SCAN_END;
    TR = tmp.TR;
    AC = tmp.AC;
    SQ = tmp.SQ;
    PREV_AA = tmp.PREV_AA;
    FULL_SQ = tmp.FULL_SQ;
    MOD_LIST = tmp.MOD_LIST;
// ORIGINAL_INTERACT_FILE = tmp.ORIGINAL_INTERACT_FILE;
    MS2_TYPE_TAG = tmp.MS2_TYPE_TAG;
  }

//////////////////////////////////////////////////
// class copy constructor of MS2Info
  MS2Info::MS2Info(const MS2Info * tmp)
  {
    ID = tmp->ID;
    PEP_PROB = tmp->PEP_PROB;
    DELTA_CN = tmp->DELTA_CN;
    XCORR = tmp->XCORR;
    THEO_MZ = tmp->THEO_MZ;
    MONO_MZ = tmp->MONO_MZ;
    NEUTRAL_MR = tmp->NEUTRAL_MR;
    CHRG = tmp->CHRG;
    SCAN_START = tmp->SCAN_START;
    SCAN_END = tmp->SCAN_END;
    TR = tmp->TR;
    AC = tmp->AC;
    SQ = tmp->SQ;
    PREV_AA = tmp->PREV_AA;
    FULL_SQ = tmp->FULL_SQ;
    MOD_LIST = tmp->MOD_LIST;
// ORIGINAL_INTERACT_FILE = tmp->ORIGINAL_INTERACT_FILE;
    MS2_TYPE_TAG = tmp->MS2_TYPE_TAG;
  }

//////////////////////////////////////////////////
// copy constructor:
  MS2Info & MS2Info::operator=(const MS2Info & tmp)
  {
    ID = tmp.ID;
    PEP_PROB = tmp.PEP_PROB;
    DELTA_CN = tmp.DELTA_CN;
    XCORR = tmp.XCORR;
    THEO_MZ = tmp.THEO_MZ;
    MONO_MZ = tmp.MONO_MZ;
    NEUTRAL_MR = tmp.NEUTRAL_MR;
    CHRG = tmp.CHRG;
    SCAN_START = tmp.SCAN_START;
    SCAN_END = tmp.SCAN_END;
    AC = tmp.AC;
    SQ = tmp.SQ;
    TR = tmp.TR;
    PREV_AA = tmp.PREV_AA;
    FULL_SQ = tmp.FULL_SQ;
    MOD_LIST = tmp.MOD_LIST;
// ORIGINAL_INTERACT_FILE = tmp.ORIGINAL_INTERACT_FILE;
    MS2_TYPE_TAG = tmp.MS2_TYPE_TAG;
    return *this;
  }

//////////////////////////////////////////////////
// equal operator:
  bool MS2Info::operator==(const MS2Info & tmp)
  {

    if (SQ == tmp.SQ)
      return true;
    else
      return false;
  }

//////////////////////////////////////////////////
// calculates the theoretical mass from a sequence:
  void MS2Info::set_THEO_MASS_from_SQ()
  {

    THEO_MZ = 0;
    int nb = 0;
    double TMP = 0;

    for (unsigned int POS = 0; POS < SQ.size(); POS++)
    {

      // check for modification:
      map<int, double>::iterator P = MOD_LIST.find(POS);
      if (P == MOD_LIST.end())
      {

        /////////////////////////////////////
        // compute THEO M/Z
        switch (SQ[POS])
        {

        case 'X':
          nb = (int) 'L' - (int) 'A';
          TMP += mono_mass[nb];
          break;

        default:
          nb = (int) SQ[POS] - (int) 'A';

          if ((nb >= 0) && (nb < 26))
          {
            TMP += mono_mass[nb];
          }
          break;
        }
        ////////////////////////////////
      }
      else
      {
        TMP += (*P).second;
      }
    }

    if (TMP > 0.0)
    {
      // add the H20 (+18) and
      // add the charge state:
      TMP += (2.0 * _MONO_H + _MONO_O);
      TMP += double(CHRG) * _MONO_H;
      TMP /= double(CHRG);
      THEO_MZ = TMP;
    }

  }

//////////////////////////////////////////////////
// get the theoretical mass from a sequence:
  double MS2Info::get_MONO_AA_MASS(int POS)
  {
    int nb = 0;
    switch (SQ[POS])
    {
    case 'X':
      nb = (int) 'L' - (int) 'A';
      break;

    default:
      nb = (int) SQ[POS] - (int) 'A';
      break;
    }

    return mono_mass[nb];
  }

//////////////////////////////////////////////////
// compute the precursor
  void MS2Info::set_MONO_MZ(double in)
  {
    MONO_MZ = in;
    NEUTRAL_MR = in;
    NEUTRAL_MR *= double(CHRG);
    NEUTRAL_MR -= double(CHRG) * _MONO_H;
  }

//////////////////////////////////////////////////
// compute the precursor
  void MS2Info::set_NEUTRAL_MR(double in)
  {
    NEUTRAL_MR = in;
    MONO_MZ = in;
    MONO_MZ += double(CHRG) * _MONO_H;
    MONO_MZ /= double(CHRG);
  }

//////////////////////////////////////////////////
// add modification
  void MS2Info::add_modification(int POS, double delta)
  {

    map<int, double>::iterator M = MOD_LIST.find(POS);
    if (M != MOD_LIST.end())
    {
      MOD_LIST.erase(M);
    }

    // add modification
    MOD_LIST.insert(pair<int, double>(POS, delta));
    // recompute the theoretical mass:
    set_THEO_MASS_from_SQ();
    // set the modified SQ:
    set_FULL_SQ();

  }

/////////////////////////////////////////////////
// show info:
  void MS2Info::show_info()
  {
    printf("\t\tMS2 ID: prec. m/z=%0.5f,theo. m/z=%0.5f,AC=%s,SQ=%s,P=%0.2f,scan=%d,tr=%0.2f,z=%d\n", MONO_MZ, THEO_MZ,
           get_AC().c_str(), get_TOTAL_SQ().c_str(), PEP_PROB, SCAN_START, TR, CHRG);
  }

/////////////////////////////////////////////////
// sets modificatied SQ:
  void MS2Info::set_FULL_SQ()
  {

    FULL_SQ.clear();

    size_t pos = 0;
    for (unsigned int i = 0; i < SQ.size(); i++)
    {
      // insert normal AA letter:
      FULL_SQ += SQ[i];
      pos++;
      // check for modifications:
      map<int, double>::iterator F = find_Modification(i);
      if (F != get_Modification_list_end())
      {
        char buffer[20];
        snprintf(buffer, 20, "[%0.4f]", (*F).second);
        FULL_SQ += buffer;
        pos += strlen(buffer);
      }
    }
  }

//////////////////////////////////////////////////
// check if AC is in here:
  bool MS2Info::find_AC(string in)
  {

    vector<string>::iterator F = find(AC.begin(), AC.end(), in);
    if (AC.end() == F)
    {
      return false;
    }
    else
    {
      return true;
    }

  }

//////////////////////////////////////////////////
// check whethere proteotypic peptide:
  bool MS2Info::get_PROTEO_TYPE()
  {

    if (AC.size() > 1)
    {
      return false;
    }
    else
    {
      return true;
    }
  }

///////////////////////////////////////////////////
// add an AC to the ms2 scan:
  void MS2Info::set_AC(string in)
  {

    vector<string>::iterator F = find(AC.begin(), AC.end(), in);
    if (F == AC.end())
    {
      AC.push_back(in);
    }
  }

///////////////////////////////////////////////////
// check if this AC or not:
  bool MS2Info::compare_AC(string in)
  {

    vector<string>::iterator F = find(AC.begin(), AC.end(), in);
    if (F != AC.end())
    {
      return true;
    }
    else
    {
      return false;
    }
  }

///////////////////////////////////////////////////
// search a pattern in the  AC list:
  bool MS2Info::search_AC_pattern(string in)
  {

    vector<string>::iterator F = AC.begin();
    while (F != AC.end())
    {
      if ((*F).find(in) != string::npos)
      {
        return true;
      }
      ++F;
    }

    return false;
  }

//////////////////////////////////////////////////
// check the tryptic state:
// 2: full tryptic
// 1: semi tryptic
// 0: non tryptic
  int MS2Info::get_TRYPTIC_STATE()
  {

    int status = 0;
    // check C-terminus:
    if ((SQ[SQ.size() - 1] == 'R') || (SQ[SQ.size() - 1] == 'K'))
    {
      status++;
    }

    // check N-terminus:
    if ((PREV_AA == "R") || (PREV_AA == "K"))
    {
      status++;
    }
    return status;
  }


  std::map<int, double>::iterator MS2Info::get_Modification_list_start()
  {
    return MOD_LIST.begin();
  }

  std::map<int, double>::iterator MS2Info::get_Modification_list_end()
  {
    return MOD_LIST.end();
  }

  std::map<int, double>::iterator MS2Info::find_Modification(int pos)
  {
    return MOD_LIST.find(pos);
  }

  std::map<int, double> * MS2Info::get_Modification_list()
  {
    return &(MOD_LIST);
  }

  bool MS2Info::check_MODIFICATION()
  {
    return !MOD_LIST.empty();
  }

  double MS2Info::get_THEO_MZ()
  {
    return THEO_MZ;
  }

  void MS2Info::set_SQ(std::string IN)
  {
    SQ = IN;
    set_THEO_MASS_from_SQ();
    set_FULL_SQ();
  }

  std::string MS2Info::get_SQ()
  {
    return SQ;
  }

  std::string MS2Info::get_MOD_SQ()
  {
    return FULL_SQ;
  }

  std::string MS2Info::get_TOTAL_SQ()
  {
    return get_PREV_AA() + "." + get_MOD_SQ();
  }

  std::string MS2Info::get_AC()
  {
    return *(AC.begin());
  }

  std::vector<std::string> MS2Info::get_ALL_AC()
  {
    return AC;
  }

  std::vector<std::string>::iterator MS2Info::get_ALL_AC_START()
  {
    return AC.begin();
  }

  std::vector<std::string>::iterator MS2Info::get_ALL_AC_END()
  {
    return AC.end();
  }

  float MS2Info::get_PEP_PROB()
  {
    return PEP_PROB;
  }

  void MS2Info::set_PEP_PROB(float IN)
  {
    PEP_PROB = IN;
  }

  double MS2Info::get_MONO_MZ()
  {
    return MONO_MZ;
  }

  double MS2Info::get_NEUTRAL_MR()
  {
    return NEUTRAL_MR;
  }

  int MS2Info::get_CHRG()
  {
    return CHRG;
  }

  void MS2Info::set_CHRG(int IN)
  {
    CHRG = IN;
  }

  int MS2Info::get_SCAN()
  {
    return SCAN_START;
  }

  int MS2Info::get_SCAN_START()
  {
    return SCAN_START;
  }

  void MS2Info::set_SCAN_START(int IN)
  {
    SCAN_START = IN;
  }

  int MS2Info::get_SCAN_END()
  {
    return SCAN_END;
  }

  void MS2Info::set_SCAN_END(int IN)
  {
    SCAN_END = IN;
  }

  int MS2Info::get_ID()
  {
    return ID;
  }

  double MS2Info::get_DELTA_CN()
  {
    return DELTA_CN;
  }

  void MS2Info::set_DELTA_CN(double IN)
  {
    DELTA_CN = IN;
  }

  double MS2Info::get_XCORR()
  {
    return XCORR;
  }

  void MS2Info::set_XCORR(double IN)
  {
    XCORR = IN;
  }

  void MS2Info::set_MS2_TYPE_TAG(std::string IN)
  {
    MS2_TYPE_TAG = IN;
  }

  std::string MS2Info::get_MS2_TYPE_TAG()
  {
    return MS2_TYPE_TAG;
  }

  double MS2Info::getRetentionTime()
  {
    return TR;
  }

  void MS2Info::setRetentionTime(double IN)
  {
    TR = IN;
  }

  std::string MS2Info::get_PREV_AA()
  {
    return PREV_AA;
  }

  void MS2Info::set_PREV_AA(std::string IN)
  {
    PREV_AA = IN;
  }

}
