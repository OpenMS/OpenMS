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
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//

#include <list>
#include <map>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundControl.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>

namespace OpenMS
{

  using namespace std;

  BackgroundControl::BackgroundControl()
  {
    init();
  }

  BackgroundControl::~BackgroundControl()
  {
    intensityBinMap.clear();
  }

  void BackgroundControl::init()
  {

    // create a vector of intensity bin objects

    // first in the tr dimension:
    double trStart = SuperHirnParameters::instance()->getMinTR();
    while (trStart <= SuperHirnParameters::instance()->getMaxTR())
    {

      // inner loop is the mzBins:
      map<double, BackgroundIntensityBin> mzArray;
      double mzStart = SuperHirnParameters::instance()->getMinFeatureMZ();
      while (mzStart <= SuperHirnParameters::instance()->getMaxFeatureMZ())
      {

        BackgroundIntensityBin * bin = new BackgroundIntensityBin(mzStart, trStart);
        mzArray.insert(make_pair(mzStart, *bin));
        delete bin;

        mzStart += SuperHirnParameters::instance()->getBackgroundIntensityBinsMZ();
      }

      intensityBinMap.insert(make_pair(trStart, mzArray));
      trStart += SuperHirnParameters::instance()->getBackgroundIntensityBinsTR();
    }

  }

  void BackgroundControl::addPeakMSScan(double TR, list<CentroidPeak> * peakList)
  {

    map<double, map<double, BackgroundIntensityBin> >::iterator F = findTrKey(TR);
    if (F != intensityBinMap.end())
    {

      // find the mz bins:
      map<double, BackgroundIntensityBin> * mzMap = &(F->second);

      list<CentroidPeak>::iterator mpi;
      for (mpi = peakList->begin(); mpi != peakList->end(); ++mpi)
      {

        map<double, BackgroundIntensityBin>::iterator F_mz = findMzKey(mpi->getMass(), mzMap);
        if (F_mz != mzMap->end())
        {
          F_mz->second.addIntensity(mpi->getIntensity());
        }
      }
    }
  }

  map<double, BackgroundIntensityBin>::iterator BackgroundControl::findMzKey(double mz,
                                                                             map<double, BackgroundIntensityBin> * mzMap)
  {

    double constraint = SuperHirnParameters::instance()->getBackgroundIntensityBinsMZ() / 2.0;
    map<double, map<double, BackgroundIntensityBin>::iterator> outMap;

    map<double, BackgroundIntensityBin>::iterator F = mzMap->lower_bound(mz);
    if (F != mzMap->end())
    {
      double delta = fabs(F->first - mz);
      if (delta <= constraint)
      {
        outMap.insert(make_pair(delta, F));
      }
    }

    if (F != mzMap->begin())
    {
      --F;
      double delta = fabs(mz - F->first);
      if (delta <= constraint)
      {
        outMap.insert(make_pair(delta, F));
      }
    }

    if (!outMap.empty())
    {
      return outMap.begin()->second;
    }

    return mzMap->end();
  }

  map<double, map<double, BackgroundIntensityBin> >::iterator BackgroundControl::findTrKey(double Tr)
  {

    double constraint = SuperHirnParameters::instance()->getBackgroundIntensityBinsTR() * 2;

    map<double, map<double, map<double, BackgroundIntensityBin> >::iterator> outMap;
    map<double, map<double, BackgroundIntensityBin> >::iterator F = intensityBinMap.lower_bound(Tr);
    if (F != intensityBinMap.end())
    {
      double delta = fabs(Tr - F->first);
      if (delta <= constraint)
      {
        outMap.insert(make_pair(delta, F));
      }

    }

    if (F != intensityBinMap.begin())
    {
      --F;

      double delta = fabs(Tr - F->first);
      if (delta <= constraint)
      {
        outMap.insert(make_pair(delta, F));
      }
    }

    if (!outMap.empty())
    {
      return outMap.begin()->second;
    }

    return intensityBinMap.end();
  }

  double BackgroundControl::getBackgroundLevel(double mz, double tr)
  {
    // find the corresponding retention time bin:
    map<double, map<double, BackgroundIntensityBin> >::iterator F = findTrKey(tr);
    if (F != intensityBinMap.end())
    {
      map<double, BackgroundIntensityBin>::iterator F2 = findMzKey(mz, &(F->second));
      if (F2 != F->second.end())
      {
        return F2->second.getMean();
      }
    }
    return -1.0;
  }

  void BackgroundControl::processIntensityMaps()
  {

    map<double, map<double, BackgroundIntensityBin> >::iterator P1 = intensityBinMap.begin();
    while (P1 != intensityBinMap.end())
    {

      map<double, BackgroundIntensityBin>::iterator P2 = P1->second.begin();
      while (P2 != P1->second.end())
      {

        P2->second.processIntensities();
        ++P2;
      }

      ++P1;
    }
  }

} // ns
