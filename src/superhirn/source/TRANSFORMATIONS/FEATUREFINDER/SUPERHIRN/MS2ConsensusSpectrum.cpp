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
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2ConsensusSpectrum.h>

#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnUtil.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Fragment.h>

using namespace std;

namespace OpenMS
{

// mass to charge tolerance for MS2 trace level:
  double MS2ConsensusSpectrum::MS2_MZ_TOLERANCE;

////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
  MS2ConsensusSpectrum::MS2ConsensusSpectrum()
  {
  }

////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
  MS2ConsensusSpectrum::MS2ConsensusSpectrum(MS2Fragment * in)
  {
    addMS2Fragment(in);
  }

////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
  MS2ConsensusSpectrum::MS2ConsensusSpectrum(double iPrecursorMZ, double iTR, int iChrg, int iApexScan)
  {

    precursorMZ = iPrecursorMZ;
    TR = iTR;
    startTR = TR;
    endTR = TR;
    z = iChrg;
    apexScan = iApexScan;

  }

//////////////////////////////////////////////////
// class desctructor of MS2ConsensusSpectrum
  MS2ConsensusSpectrum::~MS2ConsensusSpectrum()
  {
    MS2FragmentPeaks.clear();
  }

//////////////////////////////////////////////////
// class copy constructor of MS2ConsensusSpectrum
  MS2ConsensusSpectrum::MS2ConsensusSpectrum(const MS2ConsensusSpectrum & tmp)
  {
    TR = tmp.TR;
    startTR = tmp.startTR;
    endTR = tmp.endTR;
    z = tmp.z;
    apexScan = tmp.apexScan;
    startScan = tmp.startScan;
    endScan = tmp.endScan;
    precursorMZ = tmp.precursorMZ;
    MS2FragmentPeaks.clear();
    MS2FragmentPeaks = tmp.MS2FragmentPeaks;
  }

//////////////////////////////////////////////////
// class copy constructor of MS2ConsensusSpectrum
  MS2ConsensusSpectrum::MS2ConsensusSpectrum(const MS2ConsensusSpectrum * tmp)
  {
    TR = tmp->TR;
    startTR = tmp->startTR;
    endTR = tmp->endTR;
    z = tmp->z;
    apexScan = tmp->apexScan;
    startScan = tmp->startScan;
    endScan = tmp->endScan;
    precursorMZ = tmp->precursorMZ;
    MS2FragmentPeaks.clear();
    MS2FragmentPeaks = tmp->MS2FragmentPeaks;
  }

//////////////////////////////////////////////////
// copy constructor:
  MS2ConsensusSpectrum & MS2ConsensusSpectrum::operator=(const MS2ConsensusSpectrum & tmp)
  {
    TR = tmp.TR;
    startTR = tmp.startTR;
    endTR = tmp.endTR;
    z = tmp.z;
    apexScan = tmp.apexScan;
    startScan = tmp.startScan;
    endScan = tmp.endScan;
    precursorMZ = tmp.precursorMZ;
    MS2FragmentPeaks.clear();
    MS2FragmentPeaks = tmp.MS2FragmentPeaks;
    return *this;
  }

  /****** PK Never used

//////////////////////////////////////////////////
// remove outlier fragments based on their:
// MS2Fragment::OutlierAttribute = ...
// 1: retention time
// 2: precursor mass
// etc.
  void MS2ConsensusSpectrum::removeOutlierFragments()
  {

      // store in this vector the fragments:
      vector<pair<double, void*> > ValueVector;

      // store fragments by the desired attribute as outlier detection value:
      multimap<double, MS2Fragment>::iterator P = MS2FragmentPeaks.begin();
      while (P != MS2FragmentPeaks.end())
      {

          // get the attribute:
          double value = (*P).second.getOutlierDetectionAttribute();
          // store in the vector
          ValueVector.push_back(pair<double, void*>(value, &(*P).second));

          P++;
      }

      // start the iterative outlier removal by the set attribut of a ms2 fragment:
      SuperHirnUtil myMath;
      myMath.ITERATIVE_OUTLIER_DETECTION_BY_DIXON(&ValueVector);

      // convert back after the oulier detected vector:
      multimap<double, MS2Fragment> newFragments;
      vector<pair<double, void*> >::iterator I = ValueVector.begin();
      while (I != ValueVector.end())
      {
          pair<double, void*> p = (*I);
          MS2Fragment* frag = (MS2Fragment*) p.second;
          newFragments.insert(make_pair(frag->getFragmentMz(), *frag));
          I++;
      }

      // copy back:
      MS2FragmentPeaks.clear();
      MS2FragmentPeaks = newFragments;
      newFragments.clear();

  }

//////////////////////////////////////////////////
// process the stored fragments:
  void MS2ConsensusSpectrum::processConsenusSpectraFragments()
  {

      if (MS2FragmentPeaks.size() > 1)
      {
          //////////////////////////
          // remove outlier fragments based on their:
          // MS2Fragment::OutlierAttribute = ...
          // 1: retention time
          // 2: precursor mass
          // etc.
          MS2Fragment::OutlierAttribute = 1;
          removeOutlierFragments();
          computeMS2SpectrumParameters();
      }
  }
*** PK Never Used */

//////////////////////////////////////////////////
// compute MS2 parameters
  void MS2ConsensusSpectrum::computeMS2SpectrumParameters()
  {

    if (MS2FragmentPeaks.size() > 1)
    {

      double totArea = 0;
      TR = 0;
      startTR = 0;
      endTR = 0;
      precursorMZ = 0;

      double iz = 0;
      double iapexScan = 0;
      double istartScan = 0;
      double iendScan = 0;

      multimap<double, MS2Fragment>::iterator I = MS2FragmentPeaks.begin();
      while (I != MS2FragmentPeaks.end())
      {

        double thisArea = (*I).second.getFragmentPeakArea();
        totArea += thisArea;
        TR += thisArea * (*I).second.getTR();
        startTR += thisArea * (*I).second.getStartTR();
        endTR += thisArea * (*I).second.getEndTR();
        precursorMZ += thisArea * (*I).second.getPrecursorMZ();
        istartScan += thisArea * (*I).second.getStartScan();
        iendScan += thisArea * (*I).second.getEndScan();
        iapexScan += thisArea * (*I).second.getApexScan();
        iz += thisArea * (*I).second.getCHRG();

        ++I;
      }

      TR /= totArea;
      startTR /= totArea;
      endTR /= totArea;
      precursorMZ /= totArea;

      istartScan /= totArea;
      startScan = (int) istartScan;

      iendScan /= totArea;
      endScan = (int) iendScan;

      iz /= totArea;
      z = (int) iz;

      iapexScan /= totArea;
      apexScan = (int) iapexScan;

      // rearange

    }
    else
    {

      MS2Fragment * in = &(MS2FragmentPeaks.begin()->second);

      // start / end scane:
      startScan = in->getStartScan();
      endScan = in->getEndScan();

      // start /end TR borders:
      startTR = in->getStartTR();
      endTR = in->getEndTR();

      // precursor and apex:
      precursorMZ = in->getPrecursorMZ();
      TR = in->getTR();
      z = in->getCHRG();
      apexScan = in->getApexScan();

    }
  }

//////////////////////////////////////////////////
// add a MS2 fragment:
  void MS2ConsensusSpectrum::addMS2Fragment(MS2Fragment * in)
  {

    // make a mz map:
    MS2FragmentPeaks.insert(make_pair(in->getFragmentMz(), *in));

    // compute online the average retention time
    // and the precursor mass:
    computeMS2SpectrumParameters();

  }

//////////////////////////////////////////////////////
// show MS2 spectrum info:
  void MS2ConsensusSpectrum::show_info()
  {

    printf("\tMS2 consenus spectrum: m/z=%0.3f,Tr=%0.2f,scan=%d\n", precursorMZ, TR, apexScan);
  }

////////////////////////////////////////////////////////
// find a corresponding MS2 fragment
  MS2Fragment * MS2ConsensusSpectrum::findMS2Fragment(double mass)
  {

    ///////////////////////
    // collect a list of iterators with potential candidates:
    map<double, multimap<double, MS2Fragment>::iterator> candidates;

    // scan lower mass tolerance:
    multimap<double, MS2Fragment>::iterator F = MS2FragmentPeaks.lower_bound(mass);
    multimap<double, MS2Fragment>::iterator I = F;
    if (I != MS2FragmentPeaks.begin())
    {
      --I;
    }

    while (SuperHirnUtil::compareMassValuesAtPPMLevel(I->second.getFragmentMz(), mass,
                                                      MS2ConsensusSpectrum::MS2_MZ_TOLERANCE))
    {

      candidates.insert(make_pair(fabs(I->second.getFragmentMz() - mass), I));
      if (I == MS2FragmentPeaks.begin())
      {
        break;
      }

      // next:
      --I;
    }

    // scan upper mass tolerance:
    I = F;
    if ((I != MS2FragmentPeaks.end()) && (I != MS2FragmentPeaks.begin()))
    {

      while (SuperHirnUtil::compareMassValuesAtPPMLevel(I->second.getFragmentMz(), mass,
                                                        MS2ConsensusSpectrum::MS2_MZ_TOLERANCE))
      {

        candidates.insert(make_pair(fabs(I->second.getFragmentMz() - mass), I));
        ++I;
        if (I == MS2FragmentPeaks.end())
        {
          break;
        }
      }
    }

    /////////////////////////////////////////////////////
    // find now the one with the best match:
    // i.e. take the one with the smallest Mz difference:
    if (!candidates.empty())
    {
      return &((candidates.begin())->second->second);
    }

    return nullptr;

  }

  /*
//////////////////////////////////////////////////
// remove H2O loss region of the MS2 spectra
  void MS2ConsensusSpectrum::removeWaterLossRegion()
  {

      // define water loss region:
      // 3 times H2O / by charge state:
      // max = precursor mass:
      double minLossMZRegion = precursorMZ - 30;
      double maxLossMZRegion = precursorMZ;

      multimap<double, MS2Fragment>::iterator I = getMS2FragmentPeakStart();
      while (I != getMS2FragmentPeakEnd())
      {

          if ((I->second.getFragmentMz() >= minLossMZRegion) && (I->second.getFragmentMz() < maxLossMZRegion))
          {
              getMS2FragmentMap()->erase(I++);
          }
          else
          {
              I++;
          }

      }

  }
*/

//////////////////////////////////////////////////////
// copmute the similarity of the elution shape of the
// MS2 fragment to this MS2 consensus spectrum
  double MS2ConsensusSpectrum::getLCElutionPeakSimilarity(MS2Fragment * frag)
  {

    //double startTR = frag->getStartTR();
    //if (startTR > getStartTR())
    //{
    //  startTR = getStartTR();
    //}

    //double totLCSpec = getEndTR() - startTR;
    //double startLCSpec = getTR() - startTR;
    //double corSpec = startLCSpec / totLCSpec;

    //double totLCMS2 = frag->getEndTR() - startTR;
    //double startLCMS2 = frag->getTR() - startTR;
    //double corMS2 = startLCMS2 / totLCMS2;

    ///////////
    double av = fabs(getEndTR() - frag->getEndTR());
    av += fabs(getTR() - frag->getTR());
    av += fabs(getStartTR() - frag->getStartTR());
    return av;

    // return corMS2 / corSpec;
  }


  // precursor mass:
  double MS2ConsensusSpectrum::getPrecursorMZ(){return precursorMZ; }

  // TR:
  double MS2ConsensusSpectrum::getTR(){return TR; }
  // start TR
  double MS2ConsensusSpectrum::getStartTR(){return startTR; }
  // end TR
  double MS2ConsensusSpectrum::getEndTR(){return endTR; }


  // set / get  the charge state of the precursor MZ:
  void MS2ConsensusSpectrum::setPrecursorChrg(int IN){ z = IN; }
  int MS2ConsensusSpectrum::getPrecursorChrg(){ return z; }
  // apex scan:
  int MS2ConsensusSpectrum::getApexScan(){return apexScan; }
  // start scan
  int MS2ConsensusSpectrum::getStartScan(){return startScan; }
  // end scan
  int MS2ConsensusSpectrum::getEndScan(){return endScan; }
  // get the number of consensus fragments:
  int MS2ConsensusSpectrum::getNbMS2Fragments(){return (int) MS2FragmentPeaks.size(); }

  // get the MS2 fragments list iterator:
  std::multimap<double, MS2Fragment>::iterator MS2ConsensusSpectrum::getMS2FragmentPeakStart(){return MS2FragmentPeaks.begin(); }
  std::multimap<double, MS2Fragment>::iterator MS2ConsensusSpectrum::getMS2FragmentPeakEnd(){return MS2FragmentPeaks.end(); }
  std::multimap<double, MS2Fragment> * MS2ConsensusSpectrum::getMS2FragmentMap(){return &MS2FragmentPeaks; }
}
