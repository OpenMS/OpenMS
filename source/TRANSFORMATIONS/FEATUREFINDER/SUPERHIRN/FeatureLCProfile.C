// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Florian Zeller $
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


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FeatureLCProfile.h>

namespace OpenMS
{

  using namespace std;

////////////////////////////////////////////////
// constructor for the object FeatureLCProfile:
  FeatureLCProfile::FeatureLCProfile()
  {
  }

////////////////////////////////////////////////
// constructor for the object FeatureLCProfile:
  FeatureLCProfile::FeatureLCProfile(double apex_MZ, double apex_TR, double apex_Intensity, int apex_scan, int charge_state, double peak_area)
  {

    // set the apex:
    apexMS1Signal.mass = apex_MZ;
    apexMS1Signal.TR = apex_TR;
    apexMS1Signal.intensity = apex_Intensity;
    apexMS1Signal.scan = apex_scan;
    apexMS1Signal.charge = charge_state;

    // set the computed peak area:
    LCelutionArea = peak_area;

  }

////////////////////////////////////////////////
// constructor for the object FeatureLCProfile:
  FeatureLCProfile::FeatureLCProfile(double apex_MZ, double apex_TR, int charge_state, double peak_area)
  {

    // set the apex:
    apexMS1Signal.mass = apex_MZ;
    apexMS1Signal.TR = apex_TR;
    apexMS1Signal.intensity = -1;
    apexMS1Signal.scan = -1;
    apexMS1Signal.charge = charge_state;

    // set the computed peak area:
    LCelutionArea = peak_area;
  }

//////////////////////////////////////////////////
// class desctructor of FeatureLCProfile
  FeatureLCProfile::~FeatureLCProfile()
  {
    LCelutionSignals.clear();
    if (!outsideLCelutionSignals.empty())
    {
      outsideLCelutionSignals.clear();
    }
  }

//////////////////////////////////////////////////
// class copy constructor of FeatureLCProfile
  FeatureLCProfile::FeatureLCProfile(const FeatureLCProfile & tmp)
  {
    LCelutionSignals = tmp.LCelutionSignals;
    outsideLCelutionSignals = tmp.outsideLCelutionSignals;
    apexMS1Signal = tmp.apexMS1Signal;
    LCelutionArea = tmp.LCelutionArea;
  }

//////////////////////////////////////////////////
// class copy constructor of FeatureLCProfile
  FeatureLCProfile::FeatureLCProfile(const FeatureLCProfile * tmp)
  {
    LCelutionSignals = tmp->LCelutionSignals;
    outsideLCelutionSignals = tmp->outsideLCelutionSignals;
    apexMS1Signal = tmp->apexMS1Signal;
    LCelutionArea = tmp->LCelutionArea;
  }

//////////////////////////////////////////////////
// copy constructor:
  FeatureLCProfile & FeatureLCProfile::operator=(const FeatureLCProfile & tmp)
  {
    LCelutionSignals = tmp.LCelutionSignals;
    outsideLCelutionSignals = tmp.outsideLCelutionSignals;
    apexMS1Signal = tmp.apexMS1Signal;
    LCelutionArea = tmp.LCelutionArea;
    return *this;
  }

/////////////////////////////////////////////////
// add / get signals:
  void FeatureLCProfile::addMS1elutionSignal(double mass, double intensity, int scan, int charge, double TR)
  {
    MS1Signal tmp;
    tmp.mass = mass;
    tmp.intensity = intensity;
    tmp.scan = scan;
    tmp.charge = charge;
    tmp.TR = TR;
    LCelutionSignals.insert(std::make_pair(scan, tmp));

  }

/////////////////////////////////////////////////
// add / get signals:
  void FeatureLCProfile::addMS1elutionSignal(MS1Signal * in)
  {
    LCelutionSignals.insert(std::make_pair(in->scan, *in));
  }

/////////////////////////////////////////////////
// add / get signals:
  void FeatureLCProfile::addOutsideMS1elutionSignal(double mass, double intensity, int scan, int charge, double TR)
  {
    MS1Signal tmp;
    tmp.mass = mass;
    tmp.intensity = intensity;
    tmp.scan = scan;
    tmp.charge = charge;
    tmp.TR = TR;
    outsideLCelutionSignals.insert(std::make_pair(scan, tmp));

  }

/////////////////////////////////////////////////
// change all elution times by a factor:
  void FeatureLCProfile::changeElutionTimesByFactor(double factor)
  {
    apexMS1Signal.TR += factor;
    map<int, MS1Signal>::iterator P = getLCelutionSignalsStart();
    while (P != getLCelutionSignalsEnd())
    {

      P->second.TR += factor;
      P++;
    }

  }

}
