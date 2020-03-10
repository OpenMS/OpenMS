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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SHFeature.h>

#include <string>
#include <vector>
#include <map>

namespace OpenMS
{

  class SUPERHIRN_DLLAPI LCMS
  {

    ////////////////////////////////////////////////
    // declaration of the private members:

private:

    // name of the spectra:
    std::string spec_name;

    // vector of object feature:
    std::vector<SHFeature> feature_list;

    // a unique spectrum id to identify a spectrum:
    int spectrum_id;

    // MASTER RUN ID:
    int MASTER_ID;

    // the LC-MS raw data names and their IDs
    std::map<int, std::string> raw_spec_names;

    // alignment error:
    std::map<double, std::pair<double, double> > ALIGNMENT_ERROR;

    ////////////////////////////////////////////////
    // declaration of the public members:

public:

    static double MINIMAL_PEP_PROPHET_THERSHOLD;
//  static double PEP_PROPHET_THERSHOLD;

    // class destructor
    ~LCMS();

    // class constructor
    LCMS(std::string);
    LCMS();
    // copy constructor
    LCMS(const LCMS *);

    // copy constructor
    LCMS(const LCMS &);

    // show the content of the spectra
    void show_info();

    // copy constructor:
    LCMS & operator=(const LCMS &);

    // sort the features according their parent mass:
    void order_by_mass();

    // function to compare the feature mass:
    float compare_feature_mass(const void *, const void *);

    // this structure provides the function to compare
    // in the sorting algorithm:
    struct OPERATOR_MZ
    {
      // provide the compare function for sort:
      bool operator()(const SHFeature A, const SHFeature B) const
      {
        // check if they have same mass
        if (A.MONO_MZ == B.MONO_MZ)
        {
          return A.TR < B.TR;
        }
        else
        {
          return A.MONO_MZ < B.MONO_MZ;
        }
      }

    };

    // this structure provides the function to compare
    // in the sorting algorithm:
    struct OPERATOR_FeatureCompare
    {
      // provide the compare function for sort:
      bool operator()(const SHFeature A, const SHFeature B) const
      {
        // check if they have same mass
        if (A.feature_ID == B.feature_ID)
        {
          return true;
        }
        else
        {
          return false;
        }
      }

    };

    // tag the feature with the spectrum id:
    void tag_peaks_with_spectrum_ID();

    // count the number of common peaks of a given number of LC-MS:
    int get_nb_common_peaks(int);

    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class

    // get the whole feature list:
    void clear_feature_list();
    std::vector<SHFeature> get_feature_list();
    std::vector<SHFeature> * get_feature_list_reference();
    bool check_feature_list_empty();

    // access end /start of list:
    std::vector<SHFeature>::iterator get_feature_list_begin();
    std::vector<SHFeature>::iterator get_feature_list_end();

    // add a new feature to the list:
    void add_feature(SHFeature * in);

    // remove a feature from the LC/MS run by ID:
    void remove_feature_by_ID(SHFeature *);
    void remove_feature_by_ID(int);
    // remove a feature from the LC/MS run:
    void remove_feature(SHFeature *);
    void remove_feature(int i);

    // remove a feature by iterator and return the iterator to the next element
    std::vector<SHFeature>::iterator remove_feature_from_list(std::vector<SHFeature>::iterator in);

    // get number of feature added:
    unsigned int get_nb_features();

    std::string get_spec_name();
    void set_spec_name(std::string in);

    // set / get spectrum id:
    int get_spectrum_ID();
    void set_spectrum_ID(int in);

    // set the id of all features
    void setFeatureLCMSID();

    // search the list of feature for the one with input ID:
    SHFeature * find_feature_by_ID(int);

    // access the raw data names:
    void remove_raw_spec_name(int ID);
    void add_raw_spec_name(int ID, std::string name);
    bool check_raw_spec_name_empty();
    std::map<int, std::string>::iterator get_raw_spec_name_start();
    std::map<int, std::string>::iterator get_raw_spec_name_end();
    std::map<int, std::string> get_raw_spec_name_map();
    int get_nb_raw_specs();

    std::string get_raw_spec_name(int ID);

    // compare the LC/MS runs names
    bool check_LCMS_name(std::string);

    // check if this LC/MS ID is present in the raw LC/MS runs
    bool find_LC_MS_by_ID(int);

    // add the raw spectrum map:
    void add_raw_spec_name_map(std::map<int, std::string> in);

    // counts the number of ms features, which contain MS2 info:
    int get_nb_identified_features();

    // counts the number of ms features, which contain MS2 info (no thresholding)
    int get_nb_identified_features(double PepProb_T);

    //////////////////////////////////
    // access the alignment error:
    // save an error:
    void add_alignment_error(double TR, double ERROR_UP, double ERROR_DOWN);

    // get alignment error at specific TR:
    void get_alignment_error(double, double *, double *);

    // access MASTER run ID:
    void set_MASTER_ID(int in);
    int get_MASTER_ID();
  };

}

