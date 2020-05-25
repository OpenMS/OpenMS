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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Fragment.h>

#include <cstdio>

namespace OpenMS
{

// outlier selection attribute:
// 1 -> retention time
// 2 -> preciursor MZ
//
  int MS2Fragment::OutlierAttribute = 1;

  MS2Fragment::MS2Fragment()
  {}

////////////////////////////////////////////////
// constructor for the object MS2Fragment:
  MS2Fragment::MS2Fragment(double iPrecursorMZ, int iPrecursorCHRG, double iTR, int iScan, int iZ, double iFragmentMZ,
                           double iIntensityArea, int iScanStart, int iScanEnd, double iTrStart, double iTrEnd)
  {

    precursorMZ = iPrecursorMZ;
    precursorCHRG = iPrecursorCHRG;
    TR = iTR;
    scan = iScan;
    z = iZ;
    fragmentMZ = iFragmentMZ;
    intensityArea = iIntensityArea;
    scanStart = iScanStart;
    scanEnd = iScanEnd;
    trStart = iTrStart;
    trEnd = iTrEnd;

  }

////////////////////////////////////////////////
// constructor for the object MS2Fragment:
  MS2Fragment::MS2Fragment(double iPrecursorMZ, int iPrecursorCHRG, double iTR, int iScan, int iZ, double iFragmentMZ,
                           double iIntensityArea)
  {

    precursorMZ = iPrecursorMZ;
    precursorCHRG = iPrecursorCHRG;
    TR = iTR;
    scan = iScan;
    z = iZ;
    fragmentMZ = iFragmentMZ;
    intensityArea = iIntensityArea;
    scanStart = -1;
    scanEnd = -1;
    trStart = -1;
    trEnd = -1;

  }

//////////////////////////////////////////////////
// class desctructor of MS2Fragment
  MS2Fragment::~MS2Fragment()
  {
  }

//////////////////////////////////////////////////
// class copy constructor of MS2Fragment
  MS2Fragment::MS2Fragment(const MS2Fragment & tmp)
  {
    precursorMZ = tmp.precursorMZ;
    precursorCHRG = tmp.precursorCHRG;
    TR = tmp.TR;
    scan = tmp.scan;
    z = tmp.z;
    fragmentMZ = tmp.fragmentMZ;
    intensityArea = tmp.intensityArea;
    scanStart = tmp.scanStart;
    scanEnd = tmp.scanEnd;
    trStart = tmp.trStart;
    trEnd = tmp.trEnd;
  }

//////////////////////////////////////////////////
// class copy constructor of MS2Fragment
  MS2Fragment::MS2Fragment(const MS2Fragment * tmp)
  {
    precursorMZ = tmp->precursorMZ;
    precursorCHRG = tmp->precursorCHRG;
    TR = tmp->TR;
    scan = tmp->scan;
    z = tmp->z;
    fragmentMZ = tmp->fragmentMZ;
    intensityArea = tmp->intensityArea;
    scanStart = tmp->scanStart;
    scanEnd = tmp->scanEnd;
    trStart = tmp->trStart;
    trEnd = tmp->trEnd;

  }

//////////////////////////////////////////////////
// copy constructor:
  MS2Fragment & MS2Fragment::operator=(const MS2Fragment & tmp)
  {
    precursorMZ = tmp.precursorMZ;
    precursorCHRG = tmp.precursorCHRG;
    TR = tmp.TR;
    scan = tmp.scan;
    z = tmp.z;
    fragmentMZ = tmp.fragmentMZ;
    intensityArea = tmp.intensityArea;
    scanStart = tmp.scanStart;
    scanEnd = tmp.scanEnd;
    trStart = tmp.trStart;
    trEnd = tmp.trEnd;
    return *this;
  }

/////////////////////////////////////////////
// show info of the MS2 fragment
  void MS2Fragment::show_info()
  {

    // print AMRT tag info:
    printf("\tm/z=%0.2f|precursor=%0.4f|TR=%0.2f:", getFragmentMz(), getPrecursorMZ(), getTR());

    // scan/tr range:
    printf("[%d-%d],[%0.2f-%0.2f],", scanStart, scanEnd, trStart, trEnd);

    // area/intensity:
    printf("A=%0.1f", getFragmentPeakArea());

    printf("\n");

  }

////////////////////////////////////////////
// get the attribute of the fragment
// according to which outliers are removed
  double MS2Fragment::getOutlierDetectionAttribute()
  {

    switch (MS2Fragment::OutlierAttribute)
    {
    // retention time:
    case 1:
      return getTR();

    // precursor mass
    case 2:
      return getPrecursorMZ();

    default:
      break;
    }

    // otherwise use retention time:
    return getTR();
  }
  

  double MS2Fragment::getPrecursorMZ() { return precursorMZ; }
  void MS2Fragment::setPrecursorMZ(double iMZ) { precursorMZ = iMZ; }
  int MS2Fragment::getPrecursorCHRG() { return precursorCHRG; }
  double MS2Fragment::getTR() {return TR; }
  double MS2Fragment::getStartTR() {return trStart; }
  double MS2Fragment::getEndTR() {return trEnd; }
  double MS2Fragment::getFragmentMz() {return fragmentMZ; }
  void MS2Fragment::setFragmentMz(double iMz) {fragmentMZ = iMz; }
  int MS2Fragment::getCHRG() {return z; }
  int MS2Fragment::getApexScan() {return scan; }
  int MS2Fragment::getStartScan() {return scanStart; }
  int MS2Fragment::getEndScan() {return scanEnd; }
  double MS2Fragment::getFragmentPeakArea() {return intensityArea; }
  void MS2Fragment::setFragmentPeakArea(double iIntens) {intensityArea = iIntens; }
}
