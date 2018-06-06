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
/*
 *  CentroidPeak.cpp
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#include <list>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <map>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>

namespace OpenMS
{

  using namespace std;

// int CentroidPeak::sfCentroidWindowWidth = 5; // Centroid window width

// public methods

// Default constructor
  CentroidPeak::CentroidPeak()
  {
    mass_ = 0.0;
    intensity_ = 0.0;
    isotopIdx_ = 0;
    fittedIntensity_ = 0.0;
    orgIntensity_ = 0.0;
  }

// Constructs a centroid peak object with mass / intensity values
  CentroidPeak::CentroidPeak(double pMass, double pIntensity)
  {
    mass_ = pMass;
    intensity_ = pIntensity;
    isotopIdx_ = 0;
    fittedIntensity_ = 0.0;
    orgIntensity_ = pIntensity;
  }

// Constructs a centroid peak object with mass / intensity values
  CentroidPeak::CentroidPeak(double pMass, double pIntensity, double iRT)
  {
    mass_ = pMass;
    intensity_ = pIntensity;
    isotopIdx_ = 0;
    fittedIntensity_ = 0.0;
    orgIntensity_ = pIntensity;
    tr_ = iRT;
  }

// Copy constructor
  CentroidPeak::CentroidPeak(const CentroidPeak & pCentroidPeak)  // Object to copy
  {
    mass_ = pCentroidPeak.mass_;
    intensity_ = pCentroidPeak.intensity_;
    isotopIdx_ = pCentroidPeak.isotopIdx_;
    fittedIntensity_ = pCentroidPeak.fittedIntensity_;
    signalToNoise_ = pCentroidPeak.signalToNoise_;
    orgIntensity_ = pCentroidPeak.orgIntensity_;
    extraPeakInfo_ = pCentroidPeak.extraPeakInfo_;
    tr_ = pCentroidPeak.tr_;
  }

// Copy constructor
  CentroidPeak::CentroidPeak(const CentroidPeak * pCentroidPeak)  // Object to copy
  {
    mass_ = pCentroidPeak->mass_;
    intensity_ = pCentroidPeak->intensity_;
    isotopIdx_ = pCentroidPeak->isotopIdx_;
    fittedIntensity_ = pCentroidPeak->fittedIntensity_;
    signalToNoise_ = pCentroidPeak->signalToNoise_;
    orgIntensity_ = pCentroidPeak->orgIntensity_;
    extraPeakInfo_ = pCentroidPeak->extraPeakInfo_;
    tr_ = pCentroidPeak->tr_;
  }

// Destructor
  CentroidPeak::~CentroidPeak()
  {
  }

// Operators

// Copies values by assignemnt = operator
  CentroidPeak & CentroidPeak::operator=(const CentroidPeak & pCentroidPeak) // Object to be assigned
  {
    mass_ = pCentroidPeak.mass_;
    intensity_ = pCentroidPeak.intensity_;
    isotopIdx_ = pCentroidPeak.isotopIdx_;
    fittedIntensity_ = pCentroidPeak.fittedIntensity_;
    signalToNoise_ = pCentroidPeak.signalToNoise_;
    orgIntensity_ = pCentroidPeak.orgIntensity_;
    extraPeakInfo_ = pCentroidPeak.extraPeakInfo_;
    tr_ = pCentroidPeak.tr_;
    return *this;
  }

// Allows sorting objects in order of ascending mass
  bool CentroidPeak::operator<(const CentroidPeak & pCentroidPeak)  // Object to be assigned
  {
    return mass_ < pCentroidPeak.mass_;
  }

// subtract intensity
  void CentroidPeak::subtractIntensity(double pIntensity)   // intensity to be subtracted
  {
    if (intensity_ < 0.0)
      return;        // do nothing for small intensities

    if (abs(intensity_ - pIntensity) / intensity_ > SuperHirnParameters::instance()->getIntensityCV())
    {
      intensity_ -= pIntensity;       // subtract if difference is larger than stat variation (CV)
    }
    else
    {
      intensity_ = 0.0;       // if difference not stat. significant, set to zero
    }
  }

// Writes data to out stream using the << operator
  std::ostream & operator<<(std::ostream & pOut, // output stream
                            CentroidPeak & pCentroidPeak) //
  {
    pOut << std::fixed << std::setprecision(4) << pCentroidPeak.getMass() << " " << std::fixed << std::setprecision(2)
    << pCentroidPeak.getIntensity();
    return pOut;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DeconvPeak class methods

// public methods

// Default constructor
  DeconvPeak::DeconvPeak()
  {
    mass_ = 0.0;
    intensity_ = 0.0;
    isotopIdx_ = 0;
    charge_ = 0;
    nrIsotopes_ = 0;
    score_ = 0.0;
    c13MassError_ = 0.0;
  }

// Constructs a centroid peak object with mass / intensity values
  DeconvPeak::DeconvPeak(double pMass, double pIntensity, int pCharge, int pNrIsotopes, double pScore,
                         double pC13MassError)
  {
    mass_ = pMass;
    intensity_ = pIntensity;
    isotopIdx_ = 0;
    charge_ = pCharge;
    nrIsotopes_ = pNrIsotopes;
    score_ = pScore;
    c13MassError_ = pC13MassError;
  }

// Copy constructor
  DeconvPeak::DeconvPeak(const DeconvPeak & pDeconvPeak) :
    CentroidPeak(0, 0)         // Object to copy
  {
    mass_ = pDeconvPeak.mass_;
    intensity_ = pDeconvPeak.intensity_;
    isotopIdx_ = pDeconvPeak.isotopIdx_;
    charge_ = pDeconvPeak.charge_;
    nrIsotopes_ = pDeconvPeak.nrIsotopes_;
    score_ = pDeconvPeak.score_;
    c13MassError_ = pDeconvPeak.c13MassError_;
    isotopicPeaks_ = pDeconvPeak.isotopicPeaks_;
    extraPeakInfo_ = pDeconvPeak.extraPeakInfo_;

  }

// Copy constructor
  DeconvPeak::DeconvPeak(const DeconvPeak * pDeconvPeak)  // Object to copy
  {
    mass_ = pDeconvPeak->mass_;
    intensity_ = pDeconvPeak->intensity_;
    isotopIdx_ = pDeconvPeak->isotopIdx_;
    charge_ = pDeconvPeak->charge_;
    nrIsotopes_ = pDeconvPeak->nrIsotopes_;
    score_ = pDeconvPeak->score_;
    c13MassError_ = pDeconvPeak->c13MassError_;
    isotopicPeaks_ = pDeconvPeak->isotopicPeaks_;
    extraPeakInfo_ = pDeconvPeak->extraPeakInfo_;
  }

// Destructor
  DeconvPeak::~DeconvPeak()
  {
  }

// Operators

// Copies values by assignemnt = operator
  DeconvPeak & DeconvPeak::operator=(const DeconvPeak & pDeconvPeak) // Object to be assigned
  {
    mass_ = pDeconvPeak.mass_;
    intensity_ = pDeconvPeak.intensity_;
    isotopIdx_ = pDeconvPeak.isotopIdx_;
    charge_ = pDeconvPeak.charge_;
    nrIsotopes_ = pDeconvPeak.nrIsotopes_;
    score_ = pDeconvPeak.score_;
    c13MassError_ = pDeconvPeak.c13MassError_;
    isotopicPeaks_ = pDeconvPeak.isotopicPeaks_;
    extraPeakInfo_ = pDeconvPeak.extraPeakInfo_;
    return *this;
  }

// Writes data to out stream using the << operator
  ostream & operator<<(ostream & pOut, // output stream
                       DeconvPeak & pDeconvPeak) //
  {
    pOut << (CentroidPeak &)pDeconvPeak;
    pOut << " " << pDeconvPeak.getCharge() << " " << fixed << setprecision(5) << pDeconvPeak.getC13MassError();
    pOut << " " << fixed << setprecision(2) << pDeconvPeak.getScore();
    return pOut;
  }

// shows the info of the peak:
  void DeconvPeak::show_info()
  {
    printf("\tDeconvoluted Peak: mz=%.4f,I=%.4f\n", mass_, intensity_);

    if (!extraPeakInfo_.empty())
    {
      //cout<<"\t"<<extraPeakInfo<<endl;
    }

    if (!isotopicPeaks_.empty())
    {
      printf("\t");
      vector<CentroidPeak>::iterator I = isotopicPeaks_.begin();
      while (I != isotopicPeaks_.end())
      {
        printf("%0.4f(%0.0f[%0.0f]) ", (*I).getMass(), (*I).getFittedIntensity(), (*I).getOrgIntensity());
        ++I;
      }
      printf("\n");
    }
  }

// shows the info of the peak:
  void CentroidPeak::show_info()
  {
    printf("\tCentroidPeak: m/z=%.3f,I=%.4f\n", mass_, intensity_);

    if (!extraPeakInfo_.empty())
    {
      //cout<<"\t"<<extraPeakInfo<<endl;
    }
  }

}
