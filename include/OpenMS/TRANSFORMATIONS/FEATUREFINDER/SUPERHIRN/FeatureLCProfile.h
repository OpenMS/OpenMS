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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_FEATURELCPROFILE_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_FEATURELCPROFILE_H

#include <OpenMS/CONCEPT/Types.h>

#include <map>

namespace OpenMS
{

// this structure stores the monoisotopic
// signals of LC elution peak:
  struct MS1Signal
  {
    double mass;
    double TR;
    double intensity;
    int scan;
    int charge;
  };

  class OPENMS_DLLAPI FeatureLCProfile
  {


    ////////////////////////////////////////////////
    // declaration of the private members:
private:

    std::map<int, MS1Signal> LCelutionSignals;
    std::map<int, MS1Signal> outsideLCelutionSignals;


    // elution area
    double LCelutionArea;

    // apex signale:
    MS1Signal apexMS1Signal;

    ////////////////////////////////////////////////
    // declaration of the public members:

public:

    // class destructor
    ~FeatureLCProfile();

    // class constructor
    FeatureLCProfile();
    FeatureLCProfile(double, double, int, double);
    FeatureLCProfile(double, double, double, int, int, double);

    // class copy constructor
    FeatureLCProfile(const FeatureLCProfile &);
    // class copy constructor
    FeatureLCProfile(const FeatureLCProfile *);


    //////////////////////////////////////////////////
    // overload operators:
    FeatureLCProfile & operator=(const FeatureLCProfile &);
    bool operator==(const FeatureLCProfile &);
    FeatureLCProfile & operator<=(const FeatureLCProfile &);
    FeatureLCProfile & operator>=(const FeatureLCProfile &);
    FeatureLCProfile & operator<(const FeatureLCProfile &);
    FeatureLCProfile & operator>(const FeatureLCProfile &);

    // change all elution times by a factor:
    void changeElutionTimesByFactor(double);


    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class

    // add / get signals:
    void addMS1elutionSignal(double, double, int, int, double);
    void addOutsideMS1elutionSignal(double, double, int, int, double);
    void addMS1elutionSignal(MS1Signal * in);

    std::map<int, MS1Signal> * getLCelutionSignalMap(){ return &LCelutionSignals; }
    std::map<int, MS1Signal>::iterator getLCelutionSignalsStart(){ return LCelutionSignals.begin(); }
    std::map<int, MS1Signal>::reverse_iterator getLastLCelutionSignal(){ return LCelutionSignals.rbegin(); }
    std::map<int, MS1Signal>::iterator getLCelutionSignalsEnd(){ return LCelutionSignals.end(); }
    int getNbLCelutionSignals(){ return (int) LCelutionSignals.size(); }


  };

} // ns

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_FEATURELCPROFILE_H
