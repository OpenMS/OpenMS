// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectrumDeconvolution.h>

namespace OpenMS
{
  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm(FLASHDeconvHelperStructs::PrecalculatedAveragine &a) :
      DefaultParamHandler("FLASHDeconvAlgorithm"), avg(a)
  {
    prevMassBinVector = std::vector<std::vector<Size>>();
    prevMinBinLogMassVector = std::vector<double>();
  }

  FLASHDeconvAlgorithm &FLASHDeconvAlgorithm::operator=(const FLASHDeconvAlgorithm &fd)
  {
    if (this == &fd)
    {
      return *this;
    }
    //...
    return *this;
  }

  //Calcualte the nominla mass from double mass. Mutiply 0.999497 reduces the rounding error.
  int FLASHDeconvAlgorithm::getNominalMass(double m)
  {
    return (int) (m * 0.999497 + .5);
  }


  //This function is the main function for the deconvolution. Takes empty DeconvolutedSpectrum and fill it up with peakGroups.
  // DeconvolutedSpectrum contains the recursor peak group for MSn.
  //A peakGroup is the collection of peaks from a single mass (monoisotopic mass). Thus it contains peaks from different charges and iostope indices.
  void FLASHDeconvAlgorithm::getPeakGroups(DeconvolutedSpectrum &dspec,
                                           int scanNumber,
                                           int &specIndex,
                                           int &massIndex,
                                           int numOverlappedScans)
  {
    auto *spec = &(dspec.getOriginalSpectrum());
    int msLevel = spec->getMSLevel();

    //For MS2,3,.. max mass and charge ranges should be determined by precursors
    auto currentChargeRange = dspec.getCurrentMaxCharge(maxCharge) - minCharge;
    auto currentMaxMass = dspec.getCurrentMaxMass(maxMass);

    //Prepare spectrum deconvolution
    auto sd = SpectrumDeconvolution(*spec, minCharge, currentMaxMass, minMass, currentMaxMass);

    //Perform deconvolution and fill in deconvolutedSpectrum
    dspec.setPeakGroups(sd.getPeakGroupsFromSpectrum(prevMassBinVector,
                                                     prevMinBinLogMassVector,
                                                     currentChargeRange,
                                                     currentMaxMass,
                                                     numOverlappedScans,
                                                     avg, msLevel));

    //TODO positiveMode should be passed.. to everywhere.

    if (dspec.empty())
    {
      return;
    }

    //Update peakGroup information after deconvolution
    for (auto &pg : dspec)
    {
      sort(pg.peaks.begin(), pg.peaks.end());
      pg.spec = spec;
      pg.specIndex = specIndex;
      pg.scanNumber = scanNumber;
      pg.massIndex = massIndex++;
    }

    specIndex++;
  }

  void FLASHDeconvAlgorithm::updateMembers_()
  {
    minCharge = param_.getValue("min_charge");
    maxCharge = param_.getValue("max_charge");

    if (minCharge > maxCharge)
    {
      int tmp = minCharge;
      minCharge = maxCharge;
      maxCharge = tmp;
    }

    maxMass = param_.getValue("max_mass");
    minMass = param_.getValue("min_mass");
  }
}


