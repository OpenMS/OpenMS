// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
/*
 *  IsotopicDist.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_ISOTOPICDIST_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_ISOTOPICDIST_H

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>

namespace OpenMS
{

  class OPENMS_DLLAPI IsotopicDist
  {
public:

    //  static doubles fDetectableIsoFact;
    //  static doubles fIntensityCV;

    static void init();
    static bool getMatchingPeaks(std::list<CentroidPeak>::iterator, std::list<CentroidPeak>::iterator, int, double &,
                                 double, std::list<std::list<CentroidPeak>::iterator> &);
    static void subtractMatchingPeaks(std::list<std::list<CentroidPeak>::iterator> &, int, double, DeconvPeak &);
    // static void getDistribution(double,double*&,double*&);
    // static void getMassBounds(double,int,int,int,double,double&,double&);

    // static bool getDebug() {return sfDebug;}
    // static std::ostream* getDebugStream() {return sfStream;}

    // static void setDebug(bool pDebug) {sfDebug = pDebug;}
    // static void setDebugStream(std::ostream* pStream) {sfStream = pStream;}

private:

    static int getIndex(double, int);

    static double sfIsoDist10[96][20];
    static double sfIsoDist50[96][20];
    static double sfIsoDist90[96][20];
    static double sfIsoMass10[96][20];
    static double sfIsoMass50[96][20];
    static double sfIsoMass90[96][20];
    static int sfNrIsotopes[96];
    static int sfMaxMassIndex;
    static int sfMaxIsotopeIndex;
    static double sfMinMass;
    static double sfMaxMass;
    static double sfMassStep;

    //static bool sfDebug;
    //static std::ostream* sfStream;
  };

// Returns lower and upper mass bounds for all isotopic peaks
//NEVER USED
  /*
   inline void IsotopicDist::getMassBounds(
   double pMass, // m/z of monoisotopic peak
   int pCharge, // charge
   int pIdx, // index in isodist table
   int pIsotope, // 0 = mono isotopic peaks, 1 = C13 peak, ....
   double pTol, // Mass error (without calibration) of centroids
   double& pLower, // Iso dist masses as C array, monoisotopic peak is set to 0 mass
   double& pUpper) // Iso dist intensities as C array
   {
   pLower = pMass + sfIsoDist10[pIdx][pIsotope]/pCharge - pTol;
   pUpper = pMass + sfIsoMass90[pIdx][pIsotope]/pCharge + pTol;
   }
   */

// Returns index in isotopic tables
  inline int IsotopicDist::getIndex(double pMass,   // m/z of monoisotopic peak
                                    int pCharge) // charge
  {
    double diff;
    int idx;

    diff = (pMass * pCharge - sfMinMass) / sfMassStep;
    if (diff < 0)
      idx = 0;
    else if (diff < sfMaxMassIndex)
      idx = (int) ((pMass * pCharge - sfMinMass) / sfMassStep);
    else
      idx = sfMaxMassIndex;

    return idx;
  }

} // ns

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SUPERHIRN_ISOTOPICDIST_H
