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

  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm() :
      DefaultParamHandler("FLASHDeconvAlgorithm")
  {
    prevMassBinVector = std::vector<std::vector<Size>>();
    prevMinBinLogMassVector = std::vector<double>();
    defaults_.setValue("min_charge", 1, "minimum charge state (can be negative for negative mode)");
    defaults_.setValue("max_charge", 100, "maximum charge state (can be negative for negative mode)");
    defaults_.setValue("min_mass", 50.0, "minimum mass (Da)");
    defaults_.setValue("max_mass", 100000.0, "maximum mass (Da)");
    defaults_.setValue("min_peaks",
                       IntList{3, 1},
                       "minimum number of supporting peaks for MS1, 2, ...  (e.g., -min_peaks 3 2 to specify 3 and 2 for MS1 and MS2, respectively)");
    defaults_.setValue("tol",
                       DoubleList{10.0, 5.0},
                       "ppm tolerance for MS1, 2, ... (e.g., -tol 10.0 5.0 to specify 10.0 and 5.0 ppm for MS1 and MS2, respectively)");

    defaults_.setValue("min_isotope_cosine",
                       DoubleList{.75, .85},
                       "cosine threshold between avg. and observed isotope pattern for MS1, 2, ... (e.g., -min_isotope_cosine 0.8 0.6 to specify 0.8 and 0.6 for MS1 and MS2, respectively)");
    defaults_.setValue("min_charge_cosine",
                       .5,
                       "cosine threshold between per-charge-intensity and fitted gaussian distribution (applies only to MS1)");
    defaults_.setValue("max_mass_count",
                       IntList{-1, -1},
                       "maximum mass count per spec for MS1, 2, ... (e.g., -max_mass_count 100 50 to specify 100 and 50 for MS1 and MS2, respectively. -1 specifies unlimited)");
    defaults_.setValue("min_intensity", .0, "intensity threshold");
    defaults_.setValue("num_overlapped_scans", 15, "number of overlapped scans for MS1 deconvolution");
    defaultsToParam_();
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

  ///Calcualte the nominla mass from double mass. Mutiply 0.999497 reduces the rounding error.
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
                                           int &massIndex)
  {
    auto *spec = &(dspec.getOriginalSpectrum());
    int msLevel = spec->getMSLevel();

    //For MS2,3,.. max mass and charge ranges should be determined by precursors
    auto currentChargeRange = dspec.getCurrentMaxCharge(maxCharge) - minCharge;
    auto currentMaxMass = dspec.getCurrentMaxMass(maxMass);

    //Prepare spectrum deconvolution
    auto sd = SpectrumDeconvolution(spec, minCharge, currentMaxMass, minMass, currentMaxMass, intensityThreshold);

    //Perform deconvolution and fill in deconvolutedSpectrum
    dspec.setPeakGroups(sd.getPeakGroupsFromSpectrum(prevMassBinVector,
                                                     prevMinBinLogMassVector,
                                                     tolerance, binWidth,
                                                     minContinuousChargePeakCount,
                                                     currentChargeRange,
                                                     currentMaxMass,
                                                     numOverlappedScans,
                                                     minChargeCosine, maxMassCount, minIsotopeCosine,
                                                     avg, msLevel));

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

    intensityThreshold = param_.getValue("min_intensity");
    minContinuousChargePeakCount = param_.getValue("min_peaks");

    tolerance = param_.getValue("tol");

    for (auto j = 0; j < (int) tolerance.size(); j++)
    {
      tolerance[j] *= 1e-6;
      binWidth.push_back(.5 / tolerance[j]);
    }

    minIsotopeCosine = param_.getValue("min_isotope_cosine");
    minChargeCosine = param_.getValue("min_charge_cosine");
    maxMassCount = param_.getValue("max_mass_count");
    numOverlappedScans = param_.getValue("num_overlapped_scans");
    //    std::cout << minCharge << " " << maxCharge << " " << minMass << " " << maxMass << std::endl;

  }

  FLASHDeconvHelperStructs::PrecalculatedAveragine FLASHDeconvAlgorithm::getAveragine()
  {
    return avg;
  }

  void FLASHDeconvAlgorithm::calculateAveragine(bool useRNAavg)
  {
    avg = FLASHDeconvHelperStructs::calculateAveragines(maxMass, useRNAavg);
  }
}


