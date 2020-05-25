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
/*
 *  CentroidData.cpp
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>

namespace OpenMS
{

  using namespace std;

  /*
   double CentroidData::sfMassTolPpm = 10.0; // mass tolerance in ppm between isotopes
   double CentroidData::sfMassTolDa = 0.01; // mass tolerance in Da between isotopes - total mass to = mass*fMassTolPpm/1000000 + fMassTolDa
   double CentroidData::sfMinIntensity = 0.0; // peak below this values are not considered as monoisotopic peaks
   double CentroidData::sfIntensityFloor = 1.0; // intensities below this value are considered as 0
   */

  /*
   bool CentroidData::MonoIsoDebugging = false;
   double CentroidData::DebugMonoIsoMassMin = 318.8;
   double CentroidData::DebugMonoIsoMassMax = 319;
   */

// Constructor for raw data, which is then centroided. Converts it to CentroidPeak objects
  CentroidData::CentroidData(int pWindowWidth,   // width of window for quadratic fit
                             boost::shared_ptr<RawData> pRawData, // Profile data object
                             bool centMode // if data are in centroid modus
                             )
  {
    fWindowWidth = pWindowWidth;
    fNoise = 0.0;
    centroidDataModus_ = centMode;
    set(pRawData);
  }

// Constructor for raw data, which is then centroided. Converts it to CentroidPeak objects
  CentroidData::CentroidData(int pWindowWidth,   // width of window for quadratic fit
                             boost::shared_ptr<RawData> pRawData, // Profile data object
                             double iRT, bool centMode // if data are in centroid modus
                             )
  {
    fScanRetentionTime = iRT;
    fWindowWidth = pWindowWidth;
    fNoise = 0.0;
    centroidDataModus_ = centMode;
    set(pRawData);
  }

  CentroidData::~CentroidData()
  {
    fCentroidPeaks.clear();
  }

// Writes data to out stream using the << operator
  ostream & operator<<(ostream & pOut, // output stream
                       CentroidData & pCentroidData) //
  {
    list<CentroidPeak> p;
    list<CentroidPeak>::iterator pi;

    pCentroidData.get(p);
    for (pi = p.begin(); pi != p.end(); ++pi)
    {
      pOut << *pi << endl;
    }

    return pOut;
  }

// Public methods

// get centroid peak list
  void CentroidData::get(list<CentroidPeak> & pCentroidPeaks)
  {
    pCentroidPeaks = fCentroidPeaks;
  }

// Sets list of already centroided mass and intensity values. Converts them to CentroidPeak objects
// These values must be sorted with respect to increasing mass
  void CentroidData::set(vector<double> & pCentroidMasses,  // Centroided masses
                         vector<double> & pCentroidIntens // Centroided intensities
                         )
  {
    vector<double>::iterator mi, hi;

    fCentroidPeaks.clear();

    for (mi = pCentroidMasses.begin(), hi = pCentroidIntens.begin(); mi != pCentroidMasses.end(); ++mi, ++hi)
    {
      CentroidPeak peak(*mi, *hi);
      fCentroidPeaks.push_back(peak);
    }

    resetPeakGroupIter();
  }

// Sets raw profile data object. Converts it to CentroidPeak objects
  void CentroidData::set(boost::shared_ptr<RawData> pRawData   // Profile data
                         )
  {
    calcCentroids(pRawData);
    resetPeakGroupIter();
  }

// Calculates the precentile of all intensities and sets fNoise to it
  void CentroidData::setNoise(double pPrctile   // percentile, which noiselevel is set to
                              )
  {
    vector<double> intens;
    list<CentroidPeak>::iterator pi;

    for (pi = fCentroidPeaks.begin(); pi != fCentroidPeaks.end(); ++pi)
    {
      intens.push_back(pi->getIntensity());
    }

    sort(intens.begin(), intens.end());     // ascending order

    int len = (int) intens.size();

    /////////////////////////////////////////
    // modification Lukas:
    // only do if the vector length is longer than 0!
    if (len > 0)
    {
      double prt = pPrctile * len / 100.0;
      int i1 = (int) prt;
      int i2 = i1 + 1;
      if (i2 == len)
        i2--;

      // interpolate linearly between values if prt is not integer
      fNoise = (double) ((prt - i1) * intens[i1] + (1 - prt + i1) * intens[i2]);
    }
  }

  
  void CentroidData::setWidth(int pWidth)
  {
    fWindowWidth = pWidth;
  }

  int CentroidData::getWidth()
  {
    return fWindowWidth;
  }

  double CentroidData::getNoise()
  {
    return fNoise;
  }

// Removes and deletes all CentroidePeak objects with a intensity below fNoise
// setNoise has to be called first
  void CentroidData::removeNoise()
  {
    list<CentroidPeak>::iterator pi;

    for (pi = fCentroidPeaks.begin(); pi != fCentroidPeaks.end(); ++pi)
    {
      if (fNoise > pi->getIntensity())
      {
        pi = fCentroidPeaks.erase(pi);
        --pi;
      }
    }
  }

// A peak group is a set of peaks (ordered by mass) with maximal spacing of 1 + eps Da. If the gap is higher than
//  1 + eps Da, then a new peak group is formed. Deisotoping only within one peak group. Start and end iterators
// are set for the peakgroup
  bool CentroidData::getNextPeakGroup(list<CentroidPeak>::iterator & pStart,  // start position of peak group
                                      list<CentroidPeak>::iterator & pEnd) // end position of peak group
  {
    list<CentroidPeak>::iterator pi, prev;

    pi = fPeakGroupStart;
    prev = fPeakGroupStart;
    if (pi != fCentroidPeaks.end())     // FLO: Windows Fix (crash: "could not be incremented")
    {
      ++pi;
    }
    for (; pi != fCentroidPeaks.end(); ++pi, ++prev)
    {
      double eps;
      eps = SuperHirnParameters::instance()->getMassTolPpm() * pi->getMass() / 1.0e6
            + SuperHirnParameters::instance()->getMassTolDa();
      if (abs(pi->getMass() - prev->getMass()) > 1.0 + eps)
      {
        break;
      }
    }

    pStart = fPeakGroupStart;
    pEnd = pi;
    fPeakGroupStart = pi;     // store for next call

    return fPeakGroupStart != fCentroidPeaks.end();
  }

// Resets the peak group iterator to start
  void CentroidData::resetPeakGroupIter()
  {
    fPeakGroupStart = fCentroidPeaks.begin();
  }

// Private methods

// Calculates centroides of peaks from raw data
  void CentroidData::calcCentroids(boost::shared_ptr<RawData> pRawData)   // Profile data object
  { // Calculates centroide data from profile data
    // int i, hw, j;
    // double cm, toti, min_dh;
    vector<double> masses, intens;

    pRawData->get(masses, intens);
    fCentroidPeaks.clear();

    ////////////////////////////////////////////

    if (centroidDataModus_) // if data alread centroided in mzXML file
    {
      for (int i = 0; i < (int) masses.size(); i++)
      {

        double inte = intens[i];
        double mz = masses[i];

        if (inte >= SuperHirnParameters::instance()->getLowIntensityMSSignalThreshold())
        {
          CentroidPeak peak(mz, inte, fScanRetentionTime);
          fCentroidPeaks.push_back(peak);
        }

      }
    }
    else
    {

      ////////////////////////////////////////////
      // centroid raw data
      double min_dh = SuperHirnParameters::instance()->getLowIntensityMSSignalThreshold();       // min height
      int hw = fWindowWidth / 2;

      for (int i = 2; i < (int) masses.size() - 2; i++)
      {

        // Peak must be concave in the interval [i-2 .. i+2]
        // BEFORE MARKUSFIX:
        // if (intens[i]>min_dh && intens[i]>intens[i-1]+min_dh && intens[i]>=intens[i+1] && intens[i-1]>intens[i-2]+min_dh && intens[i+1]>=intens[i+2]) {
        // AFTER MARKUSFIX:
        if (intens[i] > min_dh && intens[i] > intens[i - 1] && intens[i] >= intens[i + 1]
           && intens[i - 1] > intens[i - 2] && intens[i + 1] >= intens[i + 2])
        {

          // centroid mass:
          double cm = 0.0;
          // total intensity:
          double toti = 0.0;

          // double Tinte = intens[i];
          double Tmz = masses[i];

          for (Int j = -hw; j <= hw; j++)
          {
            double inte = intens[i - j];
            double mz = masses[i - j];

            // BEFORE MARKUSFIX:
            // cm += inte*mz;
            // toti += (double) intens[i-j];
            // AFTER MARKUSFIX
            if (abs(Tmz - mz) < 0.03)
            {
              cm += inte * mz;
              toti += inte;
            }
          }
          cm = cm / toti;           // Centre of gravity = centroid

          // take the intensity at the apex of the profile peak:
          CentroidPeak peak(cm, intens[i], fScanRetentionTime);

          // Lukas: take summed intensity over all peaks:
          // CentroidPeak peak( cm, toti);

          fCentroidPeaks.push_back(peak);
        }
      }
    }
  }

}
