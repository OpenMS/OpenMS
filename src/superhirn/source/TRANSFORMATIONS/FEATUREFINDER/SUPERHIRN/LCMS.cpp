// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMS.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdio>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Info.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FeatureLCProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SHFeature.h>

namespace OpenMS
{

  using namespace std;

  double LCMS::MINIMAL_PEP_PROPHET_THERSHOLD = -3.0;
// double LCMS::PEP_PROPHET_THERSHOLD;

////////////////////////////////////////////////
// constructor
  LCMS::LCMS()
  {
    spectrum_id = -1;
    MASTER_ID = -1;
  }

////////////////////////////////////////////////
// constructor for the object LCMS:
  LCMS::LCMS(string IN_name)
  {
    spec_name = IN_name;
    spectrum_id = -1;
    MASTER_ID = -1;
  }

//////////////////////////////////////////////////
// copy constructor:
  LCMS::LCMS(const LCMS * tmp)
  {
    spec_name = tmp->spec_name;
    spectrum_id = tmp->spectrum_id;
    raw_spec_names = tmp->raw_spec_names;
    MASTER_ID = tmp->MASTER_ID;
    ALIGNMENT_ERROR = tmp->ALIGNMENT_ERROR;
    feature_list = tmp->feature_list;
  }

//////////////////////////////////////////////////
// copy constructor:
  LCMS::LCMS(const LCMS & tmp)
  {
    spec_name = tmp.spec_name;
    spectrum_id = tmp.spectrum_id;
    raw_spec_names = tmp.raw_spec_names;
    MASTER_ID = tmp.MASTER_ID;
    ALIGNMENT_ERROR = tmp.ALIGNMENT_ERROR;
    feature_list = tmp.feature_list;
  }

//////////////////////////////////////////////////
// copy constructor:
  LCMS & LCMS::operator=(const LCMS & tmp)
  {

    spec_name = tmp.spec_name;
    spectrum_id = tmp.spectrum_id;
    raw_spec_names = tmp.raw_spec_names;
    MASTER_ID = tmp.MASTER_ID;
    ALIGNMENT_ERROR = tmp.ALIGNMENT_ERROR;
    feature_list = tmp.feature_list;

    return *this;
  }

//////////////////////////////////////////////////
// class desctructor
  LCMS::~LCMS()
  {

    //  destroy the feature list:
    if (!feature_list.empty())
      feature_list.clear();

    if (!raw_spec_names.empty())
      raw_spec_names.clear();

    ALIGNMENT_ERROR.clear();

  }

//////////////////////////////////////////////////
// order the features according to parent mass:
  void LCMS::order_by_mass()
  {
    sort(feature_list.begin(), feature_list.end(), OPERATOR_MZ());
  }

//////////////////////////////////////////////////
// show the content of the spectra:
  void LCMS::show_info()
  {

    if (get_spec_name().size() > 0)
      printf("\t\t -- LC-MS name: %s ", get_spec_name().c_str());
    else
      printf("\t\t -- LC-MS ID: %d,", get_spectrum_ID());

    if (get_nb_raw_specs() != 0)
    {
      printf("[MASTER MAP ID=%d] ", get_MASTER_ID());
    }
    else
    {
      printf("[LC-MS ID=%d] ", get_spectrum_ID());
    }

    printf(" #features: %d, #MS/MS ids: %d (no Thresholding: %d)\n", get_nb_features(), get_nb_identified_features(),
           get_nb_identified_features(MINIMAL_PEP_PROPHET_THERSHOLD));

    // show the child LC/MS runs:
    map<int, string>::iterator C = get_raw_spec_name_start();
    while (C != get_raw_spec_name_end())
    {
      printf("\t\t\t - Child LC-MS: %s [ID=%d]\n", (*C).second.c_str(), (*C).first);
      ++C;
    }
    vector<SHFeature>::iterator p = feature_list.begin();
    while (p != feature_list.end())
    {
      // if((*p).get_MS2_info()){
      //(*p).show_info();
      //}
      ++p;
    }
  }

//////////////////////////////////////////////////
// search the list of feature for the one with input ID:
  SHFeature * LCMS::find_feature_by_ID(int ID)
  {

    vector<SHFeature>::iterator p = feature_list.begin();
    while (p != feature_list.end())
    {

      if (ID == (*p).get_feature_ID())
      {
        return &(*p);
      }

      ++p;
    }

    return nullptr;
  }

//////////////////////////////////////////////////
// count the number of common peaks of a given number of LC-MS:
// ONLY count the one which appear count-times
  int LCMS::get_nb_common_peaks(int count)
  {

    // search for the others:
    int common_count = 0;

    vector<SHFeature>::iterator p = feature_list.begin();
    while (p != feature_list.end())
    {

      // get the peak at a charge state:
      SHFeature * PEAK = &(*p);

      if (PEAK != nullptr)
      {
        if (PEAK->get_nb_common_match() == count)
        {
          common_count++;
        }
      }
      // next feature
      ++p;
    }

    return common_count;
  }

//////////////////////////////////////////////////////////////////
// remove a feature from the LC/MS run:
  void LCMS::remove_feature(SHFeature * in)
  {
    vector<SHFeature>::iterator P = find(feature_list.begin(), feature_list.end(), in);
    if (P != feature_list.end())
    {
      P->show_info();
      P = feature_list.erase(P);
    }
  }

//////////////////////////////////////////////////////////////////
// remove a feature from the LC/MS run by ID:
  void LCMS::remove_feature_by_ID(SHFeature * in)
  {
    remove_feature_by_ID(in->get_feature_ID());
  }

//////////////////////////////////////////////////////////////////
// remove a feature from the LC/MS run by ID:
  void LCMS::remove_feature_by_ID(int ID)
  {
    vector<SHFeature>::iterator P = feature_list.begin();
    while (P != feature_list.end())
    {
      if (P->get_feature_ID() == ID)
      {
        feature_list.erase(P);
        break;
      }

      ++P;

    }
  }

///////////////////////////////////////////////////////////
// get alignment error at specific TR:
  void LCMS::get_alignment_error(double in, double * UP, double * DOWN)
  {

    if (!ALIGNMENT_ERROR.empty())
    {

      map<double, pair<double, double> >::iterator P = ALIGNMENT_ERROR.lower_bound(in);

      // exact match:
      if ((*P).first == in)
      {
        *UP = (*P).second.first;
        *DOWN = (*P).second.second;
      }
      // is bigger than any element:
      else if (P == ALIGNMENT_ERROR.end())
      {
        --P;
        *UP = (*P).second.first;
        *DOWN = (*P).second.second;
      }
      // is smaller than any element:
      else if (P == ALIGNMENT_ERROR.begin())
      {
        *UP = (*P).second.first;
        *DOWN = (*P).second.second;
      }
      else
      {
        // is within the range of the list but has not found an exact match
        // -> extrapolate between 2 data points:
        double up_E_u = (*P).second.first;
        double down_E_u = (*P).second.second;
        double TR_u = (*P).first;

        --P;

        double up_E_d = (*P).second.first;
        double down_E_d = (*P).second.second;
        double TR_d = (*P).first;

        double w_u = (in - TR_d) / (TR_u - TR_d);
        double w_d = (TR_u - in) / (TR_u - TR_d);

        *UP = w_u * up_E_u + w_d * up_E_d;
        *DOWN = w_u * down_E_u + w_d * down_E_d;
      }

    }
  }

///////////////////////////////////////////////////////////////
// check if this LC/MS ID is present in the raw LC/MS runs
  bool LCMS::find_LC_MS_by_ID(int ID)
  {

    map<int, string>::iterator F = raw_spec_names.find(ID);
    if (F != raw_spec_names.end())
    {
      return true;
    }
    else
    {
      return false;
    }
  }

///////////////////////////////////////////////////////////////
// compare the LC/MS runs names
  bool LCMS::check_LCMS_name(string in)
  {

    if (spec_name.find(in) != string::npos)
    {
      return true;
    }

    map<int, string>::iterator F = raw_spec_names.begin();
    while (F != raw_spec_names.end())
    {
      if ((*F).second.find(in) != string::npos)
      {
        return true;
      }
      ++F;
    }

    return false;
  }

//////////////////////////////////////////////////
// set the id of all features
  void LCMS::setFeatureLCMSID()
  {

    vector<SHFeature>::iterator p = feature_list.begin();
    while (p != feature_list.end())
    {
      (*p).set_spectrum_ID(get_spectrum_ID());
      ++p;
    }
  }

  void LCMS::tag_peaks_with_spectrum_ID()
  {
    std::vector<SHFeature>::iterator p = feature_list.begin();
    while (p != feature_list.end())
    {
      (*p).set_spectrum_ID(get_spectrum_ID());
      ++p;
    }
  }

  // get the whole feature list:
  void LCMS::clear_feature_list()
  {   return feature_list.clear(); }
  std::vector<SHFeature> LCMS::get_feature_list()
  {   return feature_list; }
  std::vector<SHFeature> * LCMS::get_feature_list_reference()
  {   return &feature_list; }
  bool LCMS::check_feature_list_empty()
  {   return feature_list.empty(); }

  // access end /start of list:
  std::vector<SHFeature>::iterator LCMS::get_feature_list_begin()
  {   return feature_list.begin(); }
  std::vector<SHFeature>::iterator LCMS::get_feature_list_end()
  {   return feature_list.end(); }

  // add a new feature to the list:
  void LCMS::add_feature(SHFeature * IN)
  {

    if (IN->get_feature_ID() == -1)
    {
      IN->set_feature_ID((int) feature_list.size());
    }
    feature_list.push_back(*IN);
  }

  void LCMS::remove_feature(int i)
  {
    if (i < int(feature_list.size()))
    {
      feature_list.erase(feature_list.begin() + i);
    }
  }

  // remove a feature by iterator and return the iterator to the next element
  std::vector<SHFeature>::iterator LCMS::remove_feature_from_list(std::vector<SHFeature>::iterator IN)
  {   return feature_list.erase(IN); }

  // get number of feature added:
  unsigned int LCMS::get_nb_features()
  {   return (unsigned int) feature_list.size(); }

  std::string LCMS::get_spec_name()
  {   return spec_name; }
  void LCMS::set_spec_name(std::string IN)
  {   spec_name = IN; }

  // set / get spectrum id:
  int LCMS::get_spectrum_ID()
  {   return spectrum_id; }
  void LCMS::set_spectrum_ID(int IN)
  {   spectrum_id = IN; }

  // access the profile data names:
  void LCMS::remove_raw_spec_name(int ID)
  {   raw_spec_names.erase(ID); }
  void LCMS::add_raw_spec_name(int ID, std::string name)
  {   raw_spec_names.insert(make_pair(ID, name)); }
  bool LCMS::check_raw_spec_name_empty()
  {   return raw_spec_names.empty(); }
  std::map<int, std::string>::iterator LCMS::get_raw_spec_name_start()
  {   return raw_spec_names.begin(); }
  std::map<int, std::string>::iterator LCMS::get_raw_spec_name_end()
  {   return raw_spec_names.end(); }
  std::map<int, std::string> LCMS::get_raw_spec_name_map()
  {   return raw_spec_names; }
  int LCMS::get_nb_raw_specs()
  {   return (int) raw_spec_names.size(); }
  std::string LCMS::get_raw_spec_name(int ID)
  {
    std::map<int, std::string>::iterator p = raw_spec_names.find(ID);
    if (p == raw_spec_names.end())
    {
      return "";
    }
    return (*p).second;
  }

  // add the raw spectrum map:
  void LCMS::add_raw_spec_name_map(std::map<int, std::string> IN)
  {
    std::map<int, std::string>::iterator p = IN.begin();
    while (p != IN.end())
    {
      int ID = (*p).first;
      std::map<int, std::string>::iterator F = raw_spec_names.find(ID);
      if (F != raw_spec_names.end())
      {
        ID += (int) raw_spec_names.size();
      }
      raw_spec_names.insert(make_pair(ID, (*p).second));
      ++p;
    }
  }

  // counts the number of ms features, which contain MS2 info:
  int LCMS::get_nb_identified_features()
  {
    int count = 0;
    std::vector<SHFeature>::iterator P = get_feature_list_begin();
    while (P != get_feature_list_end())
    {
      if ((*P).get_MS2_info())
        count++;
      ++P;
    }
    return count;
  }

  // counts the number of ms features, which contain MS2 info (no thresholding)
  int LCMS::get_nb_identified_features(double PepProb_T)
  {
    int count = 0;
    std::vector<SHFeature>::iterator P = get_feature_list_begin();
    while (P != get_feature_list_end())
    {
      if ((*P).get_MS2_info(PepProb_T))
        ++count;
      ++P;
    }
    return count;
  }

  //////////////////////////////////
  // access the alignment error:
  // save an error:
  void LCMS::add_alignment_error(double TR, double ERROR_UP, double ERROR_DOWN)
  {
    std::pair<double, double> tmp(ERROR_UP, ERROR_DOWN);
    ALIGNMENT_ERROR.insert(std::pair<double, std::pair<double, double> >(TR, tmp));
  }

  // access MASTER run ID:
  void LCMS::set_MASTER_ID(int IN)
  {   MASTER_ID = IN; }
  int LCMS::get_MASTER_ID()
  {   return MASTER_ID; }

}
