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
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMSCData.h>

#include <map>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCElutionPeak.h>

using namespace std;

namespace OpenMS
{

// minimal MS peak intensity allowed:
//float LCMSCData::intensity_min_threshold;

// Constructor & destructor
///////////////////////////////////////////////////////////////////////////////

  LCMSCData::MZ_LIST_ITERATOR LCMSCData::get_DATA_start(){ return DATA.begin(); }
  LCMSCData::MZ_LIST_ITERATOR LCMSCData::get_DATA_end(){ return DATA.end(); }

  LCMSCData::LCMSCData()
  {
  }

///////////////////////////////////////////////////////////////////////////////

  LCMSCData::~LCMSCData()
  {
    DATA.clear();
  }

////////////////////////////////////////////////
// constructor for the object :
  LCMSCData::LCMSCData(const LCMSCData * tmp)
  {
    DATA = tmp->DATA;
  }

////////////////////////////////////////////////
// constructor for the object :
  LCMSCData::LCMSCData(const LCMSCData & tmp)
  {
    DATA = tmp.DATA;
  }

////////////////////////////////////////////////
// constructor for the object :
  LCMSCData & LCMSCData::operator=(const LCMSCData & tmp)
  {
    DATA = tmp.DATA;
    return *this;
  }

// Public methods
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// find data of a specific m/z:
  LCMSCData::MZ_LIST_ITERATOR LCMSCData::get_MZ_by_iterator(double MZ)
  {
    MZ_LIST_ITERATOR P = DATA.find(MZ);
    return P;
  }

///////////////////////////////////////////////////////////////////////////////
// add data into the structure:
  void LCMSCData::add_LC_elution_peak(double MZ, LCElutionPeak * in)
  {
    // get the scan apex:
    int APEX = in->get_scan_apex();

    // check if this mass is already stored:
    LCMSCData::MZ_LIST_ITERATOR P = get_MZ_by_iterator(MZ);

    if (P == get_DATA_end())
    {
      // ok, not yet inserted this mass:
      elution_peak_list tmp;
      tmp.insert(pair<int, LCElutionPeak>(APEX, *in));

      // insert into m/z list:
      DATA.insert(pair<double, elution_peak_list>(MZ, tmp));
    }
    else
    {
      // ok, already existing:
      (*P).second.insert(pair<int, LCElutionPeak>(APEX, *in));
    }
  }

///////////////////////////////////////////////////////////////////////////////
// get a list of m/z observed in a scan +/- scan_tolerance:
// return the area of the LC elution peaks
  vector<LCElutionPeak> LCMSCData::get_MZ_list(int SCAN)
  {

    int start_scan = SCAN;
    int end_scan = SCAN;

    vector<LCElutionPeak> out;

    // go through the structure and find all m/z at this scan:
    LCMSCData::MZ_LIST_ITERATOR P = get_DATA_start();

    while (P != get_DATA_end())
    {

      double this_INT = 0;
      LCElutionPeak * TMP = nullptr;

      // search around some scan region:
      for (int this_scan = start_scan; this_scan < end_scan; this_scan++)
      {

        // find elution peaks by the apex:
        elution_peak_list_ITERATOR Q = (*P).second.find(this_scan);

        if (Q != (*P).second.end())
        {
          double tmp = (*Q).second.get_total_peak_area();
          if (this_INT < tmp)
          {
            this_INT = tmp;
            TMP = &((*Q).second);
          }
        }
      }

      if ((this_INT > 0) && (this_INT >= SuperHirnParameters::instance()->getIntensityThreshold()) && (TMP != nullptr))
      {
        out.push_back(*TMP);
      }

      // next m/z:
      P++;
    }

    return out;
  }

///////////////////////////////////////////////////////////////////////////////
// get a list of m/z observed in a scan +/- scan_tolerance:
// return the area of the LC elution peaks
  vector<LCElutionPeak> LCMSCData::get_MZ_list(int SCAN, int TOL)
  {

    int start_scan = SCAN - TOL;
    int end_scan = SCAN + TOL;

    LCElutionPeak * TMP = nullptr;
    vector<LCElutionPeak> out;

    // go through the structure and find all m/z at this scan:
    MZ_LIST_ITERATOR P = get_DATA_start();

    while (P != get_DATA_end())
    {

      double this_INT = 0;

      // search around some scan region:
      for (int this_scan = start_scan; this_scan < end_scan; this_scan++)
      {

        // find elution peaks by the apex:
        elution_peak_list_ITERATOR Q = (*P).second.find(this_scan);

        if (Q != (*P).second.end())
        {
          double tmp = (*Q).second.get_total_peak_area();
          if (this_INT < tmp)
          {
            this_INT = tmp;
            TMP = &((*Q).second);
          }
        }
      }

      if ((this_INT > 0) && (this_INT >= SuperHirnParameters::instance()->getIntensityThreshold()) && (TMP != nullptr))
      {
        out.push_back(*TMP);
      }

      // next m/z:
      P++;
    }

    return out;
  }

///////////////////////////////////////////////////////////////////////////////
// get all extracted LC peaks:
  vector<LCElutionPeak *> LCMSCData::get_ALL_peak()
  {

    LCElutionPeak * TMP = nullptr;
    vector<LCElutionPeak *> out;

    // go through the structure and find all m/z at this scan:
    MZ_LIST_ITERATOR P = get_DATA_start();

    while (P != get_DATA_end())
    {

      // find elution peaks by the apex:
      elution_peak_list_ITERATOR Q = (*P).second.begin();

      while (Q != (*P).second.end())
      {
        TMP = &((*Q).second);
        // cout<<TMP->get_apex_MZ()<<endl;
        out.push_back(TMP);
        Q++;
      }

      // next m/z:
      P++;
    }

    return out;
  }

/////////////////////////////////////////////////////////////////////////////
// get a vector with all LC peaks ordered by their score:
  vector<LCElutionPeak *> LCMSCData::get_ALL_peak_ordered()
  {
    vector<LCElutionPeak *> tmp_DATA = get_ALL_peak();
    return tmp_DATA;
  }

}
