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

#include <vector>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Fragment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2ConsensusSpectrum.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ClusteredMS2ConsensusSpectrum.h>

using namespace std;

namespace OpenMS
{

////////////////////////////////////////////////
// constructor for the object ClusteredMS2ConsensusSpectrum:
  ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum(MS2ConsensusSpectrum * in) :
    MS2ConsensusSpectrum(in)
  {

    precursorMZ = in->getPrecursorMZ();
    TR = in->getTR();
    z = in->getPrecursorChrg();
    apexScan = in->getApexScan();

    this->addMS2ConsensusSpectrum(in);
  }

////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
  ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum(double iPrecursorMZ, double iTR, int iChrg,
                                                               int iApexScan) :
    MS2ConsensusSpectrum(iPrecursorMZ, iTR, iChrg, iApexScan)
  {
  }

////////////////////////////////////////////////
// constructor for the object ClusteredMS2ConsensusSpectrum:
  ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum(MS2Fragment * in) :
    MS2ConsensusSpectrum(in)
  {
    // store the scan numbers:
    MS2Scans.push_back(in->getApexScan());
  }

//////////////////////////////////////////////////
// class desctructor of ClusteredMS2ConsensusSpectrum
  ClusteredMS2ConsensusSpectrum::~ClusteredMS2ConsensusSpectrum()
  {
    // ClusteredSpectra.clear();
  }

//////////////////////////////////////////////////
// class copy constructor of ClusteredMS2ConsensusSpectrum
  ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum(const ClusteredMS2ConsensusSpectrum & tmp) :
    MS2ConsensusSpectrum(tmp)
  {
    MS2Scans = tmp.MS2Scans;
  }

//////////////////////////////////////////////////
// class copy constructor of ClusteredMS2ConsensusSpectrum
  ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum(const ClusteredMS2ConsensusSpectrum * tmp) :
    MS2ConsensusSpectrum(tmp)
  {
    MS2Scans = tmp->MS2Scans;
  }

//////////////////////////////////////////////////
// extracts fragments from a MS/MS spectra and inserts
// them into the Clustered MS/MS spectrum:
  void ClusteredMS2ConsensusSpectrum::extractFragmentsFromSpectra(MS2ConsensusSpectrum * in)
  {

    // go through the MS2 fragments and find the common one:
    multimap<double, MS2Fragment>::iterator P = in->getMS2FragmentPeakStart();
    while (P != in->getMS2FragmentPeakEnd())
    {

      // fragment search mass:
      MS2Fragment * frag = &(P->second);
      double searchMZ = frag->getFragmentMz();

      MS2Fragment * matchedFrag = this->findMS2Fragment(searchMZ);

      // in case where the fragment has been detected in a previous MS/MS
      // merge thes in one MS2 fragment
      if (matchedFrag != nullptr)
      {
        this->mergeMS2Fragments(matchedFrag, frag);
      }
      else
      {
        this->addMS2Fragment(frag);
      }

      ++P;
    }
  }

//////////////////////////////////////////////////
// add a MS2 fragment:
  void ClusteredMS2ConsensusSpectrum::addMS2ConsensusSpectrum(MS2ConsensusSpectrum * in)
  {

    // extract directly the MS/MS clustered spectrum
    extractFragmentsFromSpectra(in);

    // store the scan numbers:
    MS2Scans.push_back(in->getApexScan());

  }

//////////////////////////////////////////////////////
// merge a MS2 fragment into the target MS2 fragment:
  void ClusteredMS2ConsensusSpectrum::mergeMS2Fragments(MS2Fragment * target, MS2Fragment * toMerge)
  {

    // sum up intensities:
    target->setFragmentPeakArea(target->getFragmentPeakArea() + toMerge->getFragmentPeakArea());

    // average m/z:
    target->setFragmentMz((target->getFragmentMz() + toMerge->getFragmentMz()) / 2.0);

    // average m/z:
    target->setPrecursorMZ((target->getPrecursorMZ() + toMerge->getPrecursorMZ()) / 2.0);

  }

  int ClusteredMS2ConsensusSpectrum::getNumberOfSpectraScan()
  {
   return (int) MS2Scans.size();
  }

  std::vector<int>::iterator ClusteredMS2ConsensusSpectrum::getSpectraScanNumberStart()
  {
   return MS2Scans.begin();
  }

  std::vector<int>::iterator ClusteredMS2ConsensusSpectrum::getSpectraScanNumberEnd()
  {
    return MS2Scans.end();
  }

}
