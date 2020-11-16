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
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>
#include <Eigen/Dense>

namespace OpenMS
{
  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm() :
      DefaultParamHandler("FLASHDeconvAlgorithm")
  {
    prevMassBinVector = std::vector<std::vector<Size>>();
    prevMinBinLogMassVector = std::vector<double>();
    defaults_.setValue("min_mz", -1.0, "if set to positive value, minimum m/z to deconvolute.");
    defaults_.setValue("max_mz", -1.0, "if set to positive value, maximum m/z to deconvolute.");
    defaults_.setValue("min_RT", -1.0, "if set to positive value, minimum RT to deconvolute.");
    defaults_.setValue("max_RT", -1.0, "if set to positive value, maximum RT to deconvolute.");

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
    //defaults_.setValue("min_charge_cosine",
    //                   .5,
    //                   "cosine threshold between per-charge-intensity and fitted gaussian distribution (applies only to MS1)");
    defaults_.setValue("max_mass_count",
                       IntList{-1, -1},
                       "maximum mass count per spec for MS1, 2, ... (e.g., -max_mass_count 100 50 to specify 100 and 50 for MS1 and MS2, respectively. -1 specifies unlimited)");
    defaults_.setValue("min_intensity", .0, "intensity threshold");
    defaults_.setValue("num_overlapped_scans", 15, "number of overlapped scans for MS1 deconvolution");
    defaultsToParam_();
  }

  FLASHDeconvAlgorithm::~FLASHDeconvAlgorithm()
  {
    if (logMzPeaks.empty())
    {
      return;
    }
    std::vector<LogMzPeak>().swap(logMzPeaks);
    std::vector<PeakGroup>().swap(peakGroups);
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
    if (minRT > 0 && spec->getRT() < minRT)
    {
      return;
    }
    if (maxRT > 0 && spec->getRT() > maxRT)
    {
      return;
    }

    updateLogMzPeaks(spec);
    msLevel = spec->getMSLevel();
    //For MS2,3,.. max mass and charge ranges should be determined by precursors
    currentMaxCharge = dspec.getCurrentMaxCharge(maxCharge);
    currentMaxMass = dspec.getCurrentMaxMass(maxMass);
    //Perform deconvolution and fill in deconvolutedSpectrum
    generatePeakGroupsFromSpectrum();
    dspec.setPeakGroups(peakGroups);
    if (dspec.empty())
    {
      return;
    }
    //Update peakGroup information after deconvolution
    for (auto &pg : dspec)
    {
      sort(pg.begin(), pg.end());
      pg.spec = spec;
      pg.specIndex = specIndex;
      pg.scanNumber = scanNumber;
      pg.massIndex = massIndex++;
    }
    specIndex++;
  }

  void FLASHDeconvAlgorithm::updateMembers_()
  {
    minMz = param_.getValue("min_mz");
    maxMz = param_.getValue("max_mz");

    minRT = param_.getValue("min_RT");
    maxRT = param_.getValue("max_RT");

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
    minSupportPeakCount = param_.getValue("min_peaks");

    binWidth.clear();
    tolerance = param_.getValue("tol");

    for (auto j = 0; j < (int) tolerance.size(); j++)
    {
      tolerance[j] *= 1e-6;
      binWidth.push_back(.5 / tolerance[j]);
    }

    minIsotopeCosine = param_.getValue("min_isotope_cosine");
    //minChargeCosine = param_.getValue("min_charge_cosine");
    maxMassCount = param_.getValue("max_mass_count");
    numOverlappedScans = param_.getValue("num_overlapped_scans");
    setFilters();
  }

  FLASHDeconvHelperStructs::PrecalculatedAveragine FLASHDeconvAlgorithm::getAveragine()
  {
    return avg;
  }

  void FLASHDeconvAlgorithm::calculateAveragine(bool useRNAavg)
  {
    avg = FLASHDeconvHelperStructs::calculateAveragines(maxMass, useRNAavg);
  }


  // generate filters
  void FLASHDeconvAlgorithm::setFilters()
  {
    filter.clear();
    harmonicFilter.clear();
    auto chargeRange = maxCharge - minCharge + 1;
    for (int i = 0; i < chargeRange; i++)
    {
      filter.push_back(log(1.0 / abs(i + minCharge)));
    }

    harmonicFilter.resize(hCharges.size(), chargeRange);

    for (Size k = 0; k < hCharges.size(); k++)
    {
      auto &hc = hCharges[k];
      auto n = hc / 2;

      for (int i = 0; i < chargeRange; i++)
      {
        double factor = -1;
        //if (abs(i + minCharge) == 1)
        // {
        //  factor = 1;
        //}
        //auto harmonicFilter = log(1.0 / (i + n / hc + param.minCharge));
        //        hBinOffsets[i][k] = (long) round((filter[i] - harmonicFilter) * param.binWidth);

        harmonicFilter.setValue(k, i, log(1.0 / (factor * n / hc + abs(i + minCharge))));
      }
    }
  }

  //Generate uncharged log mz transformated peaks
  void FLASHDeconvAlgorithm::updateLogMzPeaks(MSSpectrum *spec)
  {
    std::vector<LogMzPeak>().swap(logMzPeaks);
    logMzPeaks.reserve(spec->size());
    int index = 0;
    for (auto &peak : *spec)
    {
      if (minMz > 0 && peak.getMZ() < minMz)
      {
        continue;
      }
      if (maxMz > 0 && peak.getMZ() > maxMz)
      {
        break;
      }
      if (peak.getIntensity() <= intensityThreshold)//
      {
        continue;
      }
      LogMzPeak logMzPeak(peak, index++, minCharge > 0);
      logMzPeaks.push_back(logMzPeak);
    }
  }


  double FLASHDeconvAlgorithm::getBinValue(Size bin, double minV, double binWidth)
  {
    return minV + bin / binWidth;
  }

  Size FLASHDeconvAlgorithm::getBinNumber(double v, double minV, double binWidth)
  {
    if (v < minV)
    {
      return 0;
    }
    return (Size) (((v - minV) * binWidth) + .5);

  }

  // From log mz to mz bins. To reduce edge effect, its adjacent bin is set when a bin is set.
  void FLASHDeconvAlgorithm::updateMzBins(Size &binNumber,
                                          std::vector<float> &mzBinIntensities)
  {
    mzBins = boost::dynamic_bitset<>(binNumber);
    //std::fill(mzBinIntensities.begin(), mzBinIntensities.end(), .0);

    for (auto &p : logMzPeaks)
    {
      Size bi = getBinNumber(p.logMz, mzBinMinValue, binWidth[msLevel - 1]);
      if (bi >= binNumber)
      {
        continue;
      }
      mzBins.set(bi);
      mzBinIntensities[bi] += p.intensity;
      auto delta = (p.logMz - getBinValue(bi, mzBinMinValue, binWidth[msLevel - 1]));

      if (delta > 0)
      {
        if (bi < binNumber - 1
            && !mzBins[bi + 1]
            )
        {
          mzBins.set(bi + 1);
          mzBinIntensities[bi + 1] += p.intensity;
        }
      }
      else if (delta < 0)
      {
        if (bi > 0
            && !mzBins[bi - 1]
            )
        {
          mzBins.set(bi - 1);
          mzBinIntensities[bi - 1] += p.intensity;
        }
      }
    }
  }

  //take the mass bins from previous overlapping spectra and put them in the candidate mass bins.
  void FLASHDeconvAlgorithm::unionPrevMassBins()
  {
    if (massBins.empty())
    {
      return;
    }
    for (Size i = 0; i < prevMassBinVector.size(); i++)
    {
      auto &pmb = prevMassBinVector[i];
      if (pmb.empty())
      {
        continue;
      }
      long shift = (long) (round((massBinMinValue - prevMinBinLogMassVector[i]) * binWidth[msLevel - 1]));

      for (auto &index : pmb)
      {
        long j = (long) index - shift;
        if (j < 0)
        {
          continue;
        }
        if ((Size) j >= massBins.size())
        {
          break;
        }
        massBins[j] = true;
      }
    }
  }


  //Find candidate mass bins from the current spectrum. The runtime of FLASHDeconv is deteremined by this function..
  boost::dynamic_bitset<> FLASHDeconvAlgorithm::getCandidateMassBinsForThisSpectrum(std::vector<float> &massIntensitites,
                                                                                    std::vector<float> &mzIntensities)
  {
    int chargeRange = currentMaxCharge - minCharge + 1;
    int hChargeSize = (int) hCharges.size();
    int minPeakCntr = minSupportPeakCount[msLevel - 1];
    long binEnd = (long) massBins.size();
    auto candidateMassBinsForThisSpectrum = boost::dynamic_bitset<>(massBins.size());
    // how many peaks of continuous charges per mass
    auto supportPeakPairCount = std::vector<int>(massBins.size(), 0);

    auto mzBinIndex = mzBins.find_first();
    //std::fill(massIntensitites.begin(), massIntensitites.end(), 0);

    // to calculate continuous charges, the previous charge value per mass should be stored
    auto prevCharges = std::vector<int>(massBins.size(), chargeRange + 2);

    // not just charges but intensities are stored to see the intensity fold change
    auto prevIntensities = std::vector<float>(massBins.size(), 1.0f);

    // harmonic peak intensities contribute to noise
    //auto noise = new float *[hChargeSize + 1];
    //for (auto k = 0; k < hChargeSize + 1; k++)
    //{
    //  noise[k] = new float[massBins.size()];
    //  std::fill_n(noise[k], massBins.size(), 0);
    // }
    // bin width per ms level
    auto bw = binWidth[msLevel - 1];

    // intensity change ratio should not exceed the factor.
    float factor = 4.0;
    //float factor2 = (1 - 1.0f / factor);
    float hfactor = 8.0;
    //float hfactor2 = (1 - 1.0f / hfactor);
    // scan through log mz peaks

    while (mzBinIndex != mzBins.npos)
    {
      auto &intensity = mzIntensities[mzBinIndex];
      double mz = -1.0, logMz = 0;

      // scan through charges
      for (int j = 0; j < chargeRange; j++)
      {
        // mass is given by shifting by binOffsets[j]
        long massBinIndex = mzBinIndex + binOffsets[j];

        if (massBinIndex < 0)
        {
          continue;
        }
        if (massBinIndex >= binEnd)
        {
          break;
        }

        // intensity of previous charge
        auto &prevIntensity = prevIntensities[massBinIndex];
        auto isoIntensity = .0;

        // intensity ratio between current and previous charges
        float intensityRatio = intensity / prevIntensity;
        intensityRatio = intensityRatio < 1 ? 1.0f / intensityRatio : intensityRatio;

        // check if peaks of continuous charges are present
        auto &prevCharge = prevCharges[massBinIndex];
        bool chargeNotContinous = prevCharge - j != 1;
        // if passFirstStep is true, harmonic peak filtration is done
        bool passFirstStep = false;

        // for MS2,3,... iostopic peaks or water nh3 loss peaks are considered
        if (msLevel > 1)
        {
          auto acharge = abs(j + minCharge); // abs charge

          if (mz <= 0)
          {
            logMz = getBinValue(mzBinIndex, mzBinMinValue, bw);
            mz = exp(logMz);
          }
          double diff = Constants::ISOTOPE_MASSDIFF_55K_U / acharge / mz;
          auto nextIsoMz = logMz + diff;//log(mz + Constants::C13C12_MASSDIFF_U / charge);
          auto nextIsoBin = getBinNumber(nextIsoMz, mzBinMinValue, bw);

          if (nextIsoBin < mzBins.size() && mzBins[nextIsoBin] && intensity > mzIntensities[nextIsoBin])
          {
            isoIntensity = mzIntensities[nextIsoBin];
            passFirstStep = true;
          }

          if (passFirstStep) // if isotopic peaks are present
          {
            auto waterAddMz = log(mz + 18.010565 / acharge); // 17.026549
            auto waterAddBin = getBinNumber(waterAddMz, mzBinMinValue, bw);
            float waterAddIntensity = .0;

            if (waterAddBin < mzBins.size() && mzBins[waterAddBin])
            {
              waterAddIntensity = mzIntensities[waterAddBin];
            }

            if (waterAddIntensity > 0 && waterAddIntensity < intensity)
            {
              isoIntensity += waterAddIntensity;
            }

            auto amoniaLossMz = log(mz - 17.026549 / acharge); // 17.026549
            auto amoniaLossBin = getBinNumber(amoniaLossMz, mzBinMinValue, bw);
            float amoniaLossntensity = .0;

            if (amoniaLossBin >= 0 && mzBins[amoniaLossBin])
            {
              amoniaLossntensity = mzIntensities[amoniaLossBin];
            }

            if (amoniaLossntensity > 0 && amoniaLossntensity < intensity)
            {
              isoIntensity += amoniaLossntensity;
            }
          }
          else
          {
            massIntensitites[massBinIndex] -= intensity;
          }
        }
        //for MS1, if charge not continuous and intensity change ratio is low it is noise
        //else if (prevCharge < chargeRange && chargeNotContinous && intensityRatio < factor)
        //{
        //noise[hChargeSize][massBinIndex] += intensity;
        //}

        // if charge not continous or intensity ratio is too high reset continuousChargePeakPairCount
        if (chargeNotContinous || intensityRatio > factor)
        {
          supportPeakPairCount[massBinIndex] = 0;
        }
        else
        {
          passFirstStep = true;
        }

        if (passFirstStep)
        {
          //check harmonic
          // int maxHcharge = -1;
          float maxHint = .0;
          auto highThreshold = intensity * hfactor;
          auto lowThreshold = intensity / hfactor;// / factor;
          for (auto k = 0; k < hChargeSize; k++)//
          {
            auto hmzBinIndex = massBinIndex - hBinOffsets.getValue(k, j);
            //for (auto l = -1; l <= 1; l++)
            {//hChargeSize
              //auto hmzBinIndex = _hmzBinIndex + l;
              if (hmzBinIndex > 0 && hmzBinIndex < mzBins.size() && mzBins[hmzBinIndex])
              {
                auto &hintensity = mzIntensities[hmzBinIndex];
                if (hintensity > lowThreshold
                    &&
                    hintensity < highThreshold
                    )
                {
                  //if (hintensity < maxHint)
                  //{
                  //  continue;
                  //}
                  maxHint = hintensity;
                  break;//
                  //maxHcharge = k;
                }
              }
            }
          }

          auto &cc = supportPeakPairCount[massBinIndex];
          if (maxHint > 0) //
          {
            //massIntensitites[massBinIndex] -= intensity + isoIntensity;
            //noise[maxHcharge][massBinIndex] += maxHint;
            //cc = 0;
            //candidateMassBinsForThisSpectrum[massBinIndex] = false;
          }
          else
          {
            if (cc == 0)
            {
              massIntensitites[massBinIndex] += prevIntensity;
            }

            massIntensitites[massBinIndex] += intensity + isoIntensity;
            candidateMassBinsForThisSpectrum[massBinIndex] = (++cc >= minPeakCntr);
          }
        }
        prevIntensity = intensity;
        prevCharge = j;
      }
      mzBinIndex = mzBins.find_next(mzBinIndex);
    }

    /*
    // Subtract noise intensity from each mass
    auto mindex = candidateMassBinsForThisSpectrum.find_first();
    while (mindex != candidateMassBinsForThisSpectrum.npos)
    {
      auto &s = massIntensitites[mindex];
      float maxNoise = .0;
      for (auto k = 0; k <= hChargeSize; k++)
      {
        maxNoise = std::max(maxNoise, noise[k][mindex]);
      }
      s -= maxNoise;
      mindex = candidateMassBinsForThisSpectrum.find_next(mindex);
    }
*/
    //for (auto k = 0; k <= hChargeSize; k++)
    //{
    //  delete[] noise[k];
    //}
    //delete[] noise;
    return candidateMassBinsForThisSpectrum;
  }

  // Subfunction of updateMassBins. If a peak corresponds to multiple masses, only one mass is selected baed on intensities..
  Matrix<int> FLASHDeconvAlgorithm::updateMassBins_(boost::dynamic_bitset<> &candidateMassBinsForThisSpectrum,
                                                    std::vector<float> &massIntensities,
                                                    long &binStart, long &binEnd)
  {
    //int chargeRange = param.currentChargeRange;
    auto chargeRange = currentMaxCharge - minCharge + 1;
    Matrix<int> chargeRanges(3, massBins.size(), INT_MAX);
    for (auto i = 0; i < massBins.size(); i++)
    {
      chargeRanges.setValue(1, i, INT_MIN);
    }

    auto mzBinIndex = mzBins.find_first();
    long binSize = (long) massBins.size();

    massBinsForThisSpectrum = boost::dynamic_bitset<>(massBins.size());
    auto toSkip = (candidateMassBinsForThisSpectrum | massBins).flip();
    massBins.reset();

    while (mzBinIndex != mzBins.npos)
    {
      long maxIndex = -1;
      float maxCount = -1e11;
      int charge = 0;

      for (int j = 0; j < chargeRange; j++)
      {
        long massBinIndex = mzBinIndex + binOffsets[j];

        if (massBinIndex < 0)
        {
          continue;
        }
        if (massBinIndex >= binSize)
        {
          break;
        }

        if (toSkip[massBinIndex])
        {
          continue;
        }

        auto &t = massIntensities[massBinIndex];// + noneContinuousChargePeakPairCount[massBinIndex];//

        if (t == 0)
        { // no signal
          continue;
        }

        if (maxCount < t)
        {
          maxCount = t;
          maxIndex = massBinIndex;
          charge = j;
        }

      }

      if (maxIndex > binStart && maxIndex < binEnd)
      {
        {
          chargeRanges.setValue(0, maxIndex, std::min(chargeRanges.getValue(0, maxIndex), charge));
          chargeRanges.setValue(1, maxIndex, std::max(chargeRanges.getValue(1, maxIndex), charge));
          massBinsForThisSpectrum[maxIndex] = candidateMassBinsForThisSpectrum[maxIndex];
          if (msLevel == 1)
          {
            chargeRanges.setValue(2, mzBinIndex, charge);
          }
          massBins[maxIndex] = true;
        }
      }
      mzBinIndex = mzBins.find_next(mzBinIndex);
    }
    return chargeRanges;
  }

  //update mass bins which will be used to select peaks in the input spectrum...
  Matrix<int> FLASHDeconvAlgorithm::updateMassBins(//float *massIntensities,
      std::vector<float> &mzIntensities)
  {
    auto bw = binWidth[msLevel - 1];
    long binThresholdMinMass = (long) getBinNumber(log(minMass), massBinMinValue, bw);
    long binThresholdMaxMass = (long) std::min(massBins.size(),
                                               1 + getBinNumber(log(currentMaxMass),
                                                                massBinMinValue,
                                                                bw));

    auto massIntensities = std::vector<float>(massBins.size(), 0);
    auto candidateMassBins = getCandidateMassBinsForThisSpectrum(massIntensities,
                                                                 mzIntensities);

    auto perMassChargeRanges = updateMassBins_(candidateMassBins,
                                               massIntensities,
                                               binThresholdMinMass,
                                               binThresholdMaxMass);
    return perMassChargeRanges;
  }

  //With massBins, select peaks from the same mass in the original input spectrum
  void FLASHDeconvAlgorithm::getCandidatePeakGroups(//float *massIntensities,
      Matrix<int> &chargeRanges)
  {
    const int maxMissingIsotope = 2;//msLevel == 1 ? 2 : 1;
    double bw = binWidth[msLevel - 1];
    double tol = tolerance[msLevel - 1];
    int chargeRange = currentMaxCharge - minCharge + 1;
    auto mzBinSize = mzBins.size();
    auto massBinSize = massBins.size();
    int logMzPeakSize = (int) logMzPeaks.size();
    auto currentPeakIndex = std::vector<int>(chargeRange, 0);
    peakGroups.reserve(massBins.count());
    auto massBinIndex = massBins.find_first();
    auto peakBinNumbers = std::vector<Size>(logMzPeakSize);

    for (int i = 0; i < logMzPeakSize; i++)
    {
      peakBinNumbers[i] = getBinNumber(logMzPeaks[i].logMz, mzBinMinValue, bw);
    }

    while (massBinIndex != massBins.npos)
    {
      double logM = getBinValue(massBinIndex, massBinMinValue, bw);
      double mass = exp(logM);
      /*if (false)
      {
        double diff = Constants::ISOTOPE_MASSDIFF_55K_U / mass;
        double isoLogM1 = logM - diff;
        double isoLogM2 = logM + diff;

        auto b1 = getBinNumber(isoLogM1, massBinMinValue, bw);
        auto b2 = getBinNumber(isoLogM2, massBinMinValue, bw);
        if (b1 > 0 && b1 < massBinIndex && b2 < massBinSize && b2 > massBinIndex)
        {
          if (massIntensities[b1] == 0 && massIntensities[b2] == 0) //
          {
            massBinIndex = massBins.find_next(massBinIndex);
            continue;
          }
          if (massIntensities[massBinIndex] < massIntensities[b1]
          || massIntensities[massBinIndex] < massIntensities[b2])
          {
            massBinIndex = massBins.find_next(massBinIndex);
            continue;
          }
        }*/

      /*
              if (b2 < massBinSize && b2 > massBinIndex)
              {
                if (massIntensities[massBinIndex] < massIntensities[b2])
                {
                  massBinIndex = massBins.find_next(massBinIndex);
                  continue;
                }
              }*/
      PeakGroup pg;

      pg.reserve(chargeRange * 30);
      Size rightIndex = avg.getIsotopeEndIndex(mass);
      Size leftIndex = avg.getIsotopeStartIndex(mass);

      for (int j = chargeRanges.getValue(0, massBinIndex); j <= chargeRanges.getValue(1, massBinIndex); j++)
      {
        auto &binOffset = binOffsets[j];
        auto bi = massBinIndex - binOffset;

        if (bi >= mzBinSize || (chargeRanges.getValue(2, bi) < chargeRange && chargeRanges.getValue(2, bi) != j))
        {
          continue;
        }

        double maxIntensity = -1.0;
        int charge = j + minCharge;
        auto &cpi = currentPeakIndex[j];
        int maxPeakIndex = -1;

        while (cpi < logMzPeakSize - 1)
        {
          if (peakBinNumbers[cpi] == bi)
          {
            auto intensity = logMzPeaks[cpi].intensity;
            if (intensity > maxIntensity)
            {
              maxIntensity = intensity;
              maxPeakIndex = cpi;
            }

          }
          else if (peakBinNumbers[cpi] > bi)
          {
            break;
          }
          cpi++;
        }

        if (maxPeakIndex < 0)
        {
          continue;
        }

        const double mz = logMzPeaks[maxPeakIndex].mz;
        const double isof = Constants::ISOTOPE_MASSDIFF_55K_U / abs(charge);
        double mzDelta = tol * mz * 2; //

        if (pg.perChargeInfo.find(charge) == pg.perChargeInfo.end())
        {
          pg.perChargeInfo[charge] = std::vector<float>(3);
          pg.perChargeInfo[charge].push_back(0);
        }
        auto &np = pg.perChargeInfo[charge][0];

        int pi = 0;
        // int peakcntr = 0;
        for (int peakIndex = maxPeakIndex; peakIndex < logMzPeakSize; peakIndex++)
        {
          const double observedMz = logMzPeaks[peakIndex].mz;
          const double intensity = logMzPeaks[peakIndex].intensity;
          //observedMz = mz + isof * i * d - d * mzDelta;
          double di = observedMz - mz;

          int i = (int) (.5 + di / isof);

          if (i > (int) rightIndex)
          {
            break;
          }

          if (i - pi > maxMissingIsotope)
          {
            break;
          }

          if (abs(di - i * isof) >= mzDelta) // noise
          {
            np += intensity * intensity;
            //peakcntr++;
          }
          else
          {
            const auto bin = peakBinNumbers[peakIndex] + binOffset;
            if (bin < massBinSize)
            {
              LogMzPeak p(logMzPeaks[peakIndex]);
              p.charge = charge;
              pg.push_back(p);
            }
            pi = i;
          }
        }

        pi = 0;

        for (int peakIndex = maxPeakIndex - 1; peakIndex >= 0; peakIndex--)
        {
          const double observedMz = logMzPeaks[peakIndex].mz;
          const double intensity = logMzPeaks[peakIndex].intensity;

          //observedMz = mz + isof * i * d - d * mzDelta;
          double di = mz - observedMz;
          int i = (int) (.5 + di / isof);

          if (i > (int) leftIndex)
          {
            break;
          }

          if (i - pi > maxMissingIsotope)
          {
            break;
          }

          if (abs(di - i * isof) >= mzDelta)
          {
            np += intensity * intensity;
            //continue;
          }
          else
          {
            const auto bin = peakBinNumbers[peakIndex] + binOffset;

            if (bin < massBinSize)
            {
              LogMzPeak p(logMzPeaks[peakIndex]);
              p.charge = charge;
              pg.push_back(p);
            }

            pi = i;
          }
        }
      }

      if (!pg.empty())
      {
        double maxIntensity = -1.0;
        double tmaxMass = .0;
        auto newPeaks = std::vector<LogMzPeak>();
        newPeaks.reserve(pg.size());
        for (auto &p : pg)
        {
          if (maxIntensity < p.intensity)
          {
            maxIntensity = p.intensity;
            tmaxMass = p.getUnchargedMass();
          }
        }
        double isoDelta = tol * tmaxMass;
        int minOff = 10000;
        for (auto &p : pg)
        {
          p.isotopeIndex = round((p.getUnchargedMass() - tmaxMass) / Constants::ISOTOPE_MASSDIFF_55K_U);
          if (abs(tmaxMass - p.getUnchargedMass() + Constants::ISOTOPE_MASSDIFF_55K_U * p.isotopeIndex) >
              isoDelta)
          {
            continue;
          }
          newPeaks.push_back(p);
          minOff = minOff > p.isotopeIndex ? p.isotopeIndex : minOff;
        }

        pg.swap(newPeaks);
        std::vector<LogMzPeak>().swap(newPeaks);

        for (auto &p : pg)
        {
          p.isotopeIndex -= minOff;
        }
        pg.updateMassesAndIntensity(avg);
        peakGroups.push_back(pg); //
      }
      massBinIndex = massBins.find_next(massBinIndex);
    }
  }

  bool FLASHDeconvAlgorithm::empty()
  {
    return logMzPeaks.empty();
  }

  //spectral deconvolution main function
  void FLASHDeconvAlgorithm::generatePeakGroupsFromSpectrum()
  {
    std::vector<PeakGroup>().swap(peakGroups);
    auto minPeakCntr =
        minSupportPeakCount[msLevel - 1];
    auto currentChargeRange = currentMaxCharge - minCharge + 1;
    auto tmp = currentChargeRange - minPeakCntr;
    tmp = tmp < 0 ? 0 : tmp;
    double massBinMaxValue = std::min(
        logMzPeaks[logMzPeaks.size() - 1].logMz -
        filter[tmp],
        log(currentMaxMass));

    auto bw = binWidth[msLevel - 1];
    tmp = minPeakCntr - 1;
    tmp = tmp < 0 ? 0 : tmp;
    massBinMinValue = logMzPeaks[0].logMz - filter[tmp];
    mzBinMinValue = logMzPeaks[0].logMz;

    double mzBinMaxValue = logMzPeaks[logMzPeaks.size() - 1].logMz;
    Size massBinNumber = getBinNumber(massBinMaxValue, massBinMinValue, bw) + 1;

    for (int i = 0; i < currentChargeRange; i++)
    {
      binOffsets.push_back((int) round((mzBinMinValue - filter[i] - massBinMinValue) * bw));
    }

    hBinOffsets.resize(hCharges.size(), currentChargeRange);
    for (Size k = 0; k < hCharges.size(); k++)
    {
      std::vector<int> _hBinOffsets;
      for (int i = 0; i < currentChargeRange; i++)
      {
        hBinOffsets
            .setValue(k, i, (int) round((mzBinMinValue - harmonicFilter.getValue(k, i) - massBinMinValue) * bw));
      }
    }

    Size mzBinNumber = getBinNumber(mzBinMaxValue, mzBinMinValue, bw) + 1;
    auto mzBinIntensities = std::vector<float>(mzBinNumber, .0f);

    updateMzBins(mzBinNumber, mzBinIntensities);
    //auto *massIntensities = new float[massBinNumber];

    massBins = boost::dynamic_bitset<>(massBinNumber);
    massBinsForThisSpectrum = boost::dynamic_bitset<>(massBinNumber);

    if (msLevel == 1)
    {
      unionPrevMassBins();
    }
    auto perMassChargeRanges = updateMassBins(//massIntensities,
        mzBinIntensities);

    getCandidatePeakGroups(//massIntensities,
        perMassChargeRanges);

    scoreAndFilterPeakGroups();
    //reassignPeaksinPeakGroups();
    //scoreAndFilterPeakGroups();

    removeOverlappingPeakGroups(tolerance[msLevel - 1]);
    removeHarmonicPeakGroups(tolerance[msLevel - 1]); //



    if (msLevel == 1)
    {
      if (!prevMassBinVector.empty() && prevMassBinVector.size() >= (Size) numOverlappedScans)//
      {
        prevMassBinVector.erase(prevMassBinVector.begin());
        prevMinBinLogMassVector.erase(prevMinBinLogMassVector.begin());
      }

      std::vector<Size> mb;
      mb.reserve(peakGroups.size());
      for (auto &pg : peakGroups)//filteredPeakGroups
      {
        pg.shrink_to_fit();
        //if (massBinsForThisSpectrum[pg.massBinIndex])
        //{
        auto pgBin = getBinNumber(pg.avgMass, massBinMinValue, binWidth[msLevel - 1]);
        mb.push_back(pgBin);
        //  }
      }

      prevMassBinVector.push_back(mb); //
      prevMinBinLogMassVector.push_back(massBinMinValue);
      prevMassBinVector.shrink_to_fit();
      prevMinBinLogMassVector.shrink_to_fit();
    }
  }


  double
  FLASHDeconvAlgorithm::getCosine(std::vector<double> &a,
                                  int &aStart,
                                  int &aEnd,
                                  IsotopeDistribution &b,
                                  int &bSize,
                                  double &bNorm,
                                  int offset)
  {
    double n = .0, d1 = .0;
    //int c = 0;
    for (int j = aStart; j <= aEnd; j++)
    {
      d1 += a[j] * a[j];
      int i = j - offset;
      if (i < 0 || i >= bSize)
      {
        continue;
      }
      n += a[j] * b[i].getIntensity(); //
    }
    double d = (d1 * bNorm);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }

  double FLASHDeconvAlgorithm::getCosine(const double *a, double *b, Size size)
  {
    double n = .0, d1 = .0, d2 = .0;
    //int overlapCntr = 0;
    for (Size j = 0; j < size; j++)
    {
      d1 += a[j] * a[j];
      d2 += b[j] * b[j];
      n += a[j] * b[j];
      //  if(a[j] > 0 && b[j] > 0) overlapCntr++;
    }

    //if(overlapCntr < 2) return 0; //
    double d = (d1 * d2);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }


  double FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                                        std::vector<double> &perIsotopeIntensities,
                                                                        int &offset,
                                                                        PrecalculatedAveragine &avg)
  {
    auto iso = avg.get(mass);
    auto isoNorm = avg.getNorm(mass);

    int isoSize = (int) iso.size();

    offset = 0;
    double maxCosine = -1;
    int isotopeLength = 0;
    int maxIsotopeIndex = 0, minIsotopeIndex = -1;

    for (int i = 0; i < avg.maxIsotopeIndex; i++)
    {
      if (perIsotopeIntensities[i] <= 0)
      {
        continue;
      }
      isotopeLength++;
      maxIsotopeIndex = i;
      if (minIsotopeIndex < 0)
      {
        minIsotopeIndex = i;
      }
    }


    auto maxCntr = 0;
    for (int f = -isoSize - minIsotopeIndex; f <= maxIsotopeIndex; f++)
    {
      auto cos = getCosine(perIsotopeIntensities,
                           minIsotopeIndex,
                           maxIsotopeIndex,
                           iso,
                           isoSize,
                           isoNorm,
                           f);

      if (maxCosine <= cos)
      {
        if (maxCosine == cos)
        {
          maxCntr++;
          offset += f;
        }
        else
        {
          maxCosine = cos;
          maxCntr = 1;
          offset = f;
        }
      }
    }
    offset /= maxCntr;
    return maxCosine;
  }


  bool FLASHDeconvAlgorithm::checkChargeDistribution(std::vector<double> &perChargeIntensity)
  {
    double maxPerChargeIntensity = .0;
    int nonZeroStart = -1, nonZeroEnd = 0;
    auto chargeRange = currentMaxCharge - minCharge + 1;
    for (int i = 0; i < chargeRange; i++)
    {
      if (perChargeIntensity[i] > 0)
      {
        maxPerChargeIntensity = std::max(maxPerChargeIntensity, perChargeIntensity[i]);
        if (nonZeroStart < 0)
        {
          nonZeroStart = i;
        }
        nonZeroEnd = i;
      }
    }

    int prevCharge = nonZeroStart;

    int n_r = 0;
    double factor = 4.0;
    double intThreshold = maxPerChargeIntensity / factor;//intensities[intensities.size()*95/100] / 5.0;
    for (int k = prevCharge + 1; k <= nonZeroEnd; k++)
    {
      if (perChargeIntensity[k] <= intThreshold)
      {
        continue;
      }

      if (k - prevCharge == 1)
      {
        n_r++;
      }
      if (n_r >= minSupportPeakCount[0])
      {
        return true;
      }
      prevCharge = k;
    }

    return false;
  }

  void FLASHDeconvAlgorithm::scoreAndFilterPeakGroups()
  {
    std::vector<PeakGroup> filteredPeakGroups;
    filteredPeakGroups.reserve(peakGroups.size());
    double threshold = .0;
    auto chargeRange = currentMaxCharge - minCharge + 1;

    auto mc = maxMassCount.size() > msLevel - 1 ? maxMassCount[msLevel - 1] : -1;
    if (mc > 0)
    {
      std::vector<double> intensities;
      intensities.reserve(peakGroups.size());

      for (auto &pg : peakGroups)
      {
        //pg.updateMassesAndIntensity(avg);
        intensities.push_back(pg.intensity);
      }

      if (intensities.size() > (Size) mc)
      {
        sort(intensities.begin(), intensities.end());
        threshold = intensities[intensities.size() - mc];
      }
      //std::vector<double>().swap(intensities);
    }

    auto perIsotopeIntensity = std::vector<double>(avg.maxIsotopeIndex);
    auto perChargeIntensity = std::vector<double>(chargeRange);

    for (auto &pg : peakGroups)
    {
      if (pg.intensity < threshold)
      {
        continue; //
      }

      auto indices = updatePerChargeIsotopeIntensity(
          perIsotopeIntensity, perChargeIntensity,
          avg.maxIsotopeIndex, pg);

      //pg.chargeCosineScore = getChargeFitScore(perChargeIntensity, chargeRange);

      if (msLevel == 1)
      {
        //if (pg.empty() ||
        //    pg.chargeCosineScore <= minChargeCosine)
        //{
        //continue;
        //}

        bool isChargeWellDistributed = checkChargeDistribution(perChargeIntensity);

        if (!isChargeWellDistributed)
        {
          continue;
        }

      }


      int offset = 0;
      pg.isotopeCosineScore = getIsotopeCosineAndDetermineIsotopeIndex(pg[0].getUnchargedMass(),
                                                                       perIsotopeIntensity,
                                                                       offset, avg);

      if (pg.empty() ||
          (pg.isotopeCosineScore <=
           minIsotopeCosine[msLevel - 1]))// (msLevel <= 1 ? param.minIsotopeCosineSpec : param.minIsotopeCosineSpec2)))
      {
        continue;
      }

      pg.updateMassesAndIntensity(avg, offset, avg.maxIsotopeIndex);

      auto iso = avg.get(pg.monoisotopicMass);
      auto isoNorm = avg.getNorm(pg.monoisotopicMass);
      int isoSize = (int) iso.size();
      float totalNoise = .0;
      float totalSignal = .0;
      auto perChargeMaxIntensity = std::vector<double>(chargeRange);

      for (auto charge = pg.minCharge; charge <= pg.maxCharge; charge++)
      {
        int j = charge - minCharge;
        auto perIsotopeIntensities = std::vector<double>(avg.maxIsotopeIndex, 0);

        int minIsotopeIndex = avg.maxIsotopeIndex;
        int maxIsotopeIndex = 0;

        //double minMz = pg.monoisotopicMass * 2;
        //double maxMz = 0;
        double maxIntensity = .0;
        double sumIntensity = .0;
        double sp = .0;

        for (auto &p:pg)
        {
          if (p.charge != charge)
          {
            continue;
          }

          if (p.isotopeIndex > isoSize)
          {
            continue;
          }

          perIsotopeIntensities[p.isotopeIndex] += p.intensity;
          sumIntensity += p.intensity;
          minIsotopeIndex = minIsotopeIndex < p.isotopeIndex ? minIsotopeIndex : p.isotopeIndex;
          maxIsotopeIndex = maxIsotopeIndex < p.isotopeIndex ? p.isotopeIndex : maxIsotopeIndex;

          //minMz = minMz < p.mz ? minMz : p.mz;
          // maxMz = maxMz > p.mz ? maxMz : p.mz;
          if (maxIntensity < p.intensity)
          {
            maxIntensity = p.intensity;
            perChargeMaxIntensity[j] = maxIntensity;
          }
          // sp += p.intensity * p.intensity;
        }
        if (maxIntensity <= 0)
        {
          continue;
        }

        for (int k = minIsotopeIndex; k <= maxIsotopeIndex; ++k)
        {
          if (k > isoSize)
          {
            break;
          }
          sp += perIsotopeIntensities[k] * perIsotopeIntensities[k];
        }

        auto cos = getCosine(perIsotopeIntensities,
                             minIsotopeIndex,
                             maxIsotopeIndex,
                             iso,
                             isoSize,
                             isoNorm,
                             0);

        double cos2 = cos * cos;

        if (pg.perChargeInfo.find(charge) == pg.perChargeInfo.end())
        {
          pg.perChargeInfo[charge] = std::vector<float>(3);
          pg.perChargeInfo[charge][0] = 0;
        }
        pg.perChargeInfo[charge][1] = cos;
        pg.perChargeInfo[charge][2] = sumIntensity;

        auto dno = (1 - cos2) * sp + pg.perChargeInfo[charge][0] + 1;
        auto no = cos2 * sp + 1;

        pg.perChargeInfo[charge][0] = no / dno;

        totalNoise += dno;
        totalSignal += no;
      }

      pg.totalSNR = totalSignal / totalNoise;
      pg.qScore = -10000;

      for (auto &item : pg.perChargeInfo)
      {
        auto charge = item.first;
        int j = charge - minCharge;

        auto score = QScore::getQScore(&pg, perChargeMaxIntensity[j], charge);

        if (score < pg.qScore)
        {
          continue;
        }
        //pg.maxScorePeakIntensity = pg.perChargeSNR[charge];
        pg.maxQScoreCharge = charge;
        pg.qScore = score;
      }

      pg.maxQScoreMzStart = pg.monoisotopicMass * 2;
      pg.maxQScoreMzEnd = 0;
      for (auto &p:pg)
      {
        if (p.charge != pg.maxQScoreCharge)
        {
          continue;
        }
        if (p.isotopeIndex > isoSize)
        {
          continue;
        }

        pg.maxQScoreMzStart = pg.maxQScoreMzStart < p.mz ? pg.maxQScoreMzStart : p.mz;
        pg.maxQScoreMzEnd = pg.maxQScoreMzEnd > p.mz ? pg.maxQScoreMzEnd : p.mz;
      }
      if (pg.maxQScoreMzStart > pg.maxQScoreMzEnd)
      {
        continue;
      }
      if (msLevel == 1 || pg.totalSNR > 0.1) //
        //if(msLevel == 1 || getPeakGroupScore(pg)>.5) //
      {
        filteredPeakGroups.push_back(pg);
      }
    }

    peakGroups.swap(filteredPeakGroups);
    std::vector<PeakGroup>().swap(filteredPeakGroups);

    if (msLevel > 1)
    {
      filterPeakGroupsByIsotopeCosine(mc);
    }
    else
    {
      filterPeakGroupsByQScore(mc);
    }
  }


  void FLASHDeconvAlgorithm::filterPeakGroupsByIsotopeCosine(int currentMaxMassCount)
  {
    if (currentMaxMassCount <= 0 || peakGroups.size() <= (Size) currentMaxMassCount)
    {
      return;
    }

    Size mc = (Size) currentMaxMassCount;
    std::vector<double> scores;
    scores.reserve(peakGroups.size());
    for (auto &pg : peakGroups)
    {
      scores.push_back(pg.isotopeCosineScore);
    }

    sort(scores.begin(), scores.end());

    auto newPeakGroups = std::vector<PeakGroup>();
    newPeakGroups.reserve(peakGroups.size());
    auto threshold = scores[scores.size() - mc];
    for (auto &pg : peakGroups)
    {
      if (newPeakGroups.size() > mc)
      {
        break;
      }

      if (pg.isotopeCosineScore >= threshold)
      {
        newPeakGroups.push_back(pg);
      }
    }
    std::vector<PeakGroup>().swap(peakGroups);
    newPeakGroups.swap(peakGroups);
  }

  void FLASHDeconvAlgorithm::filterPeakGroupsByQScore(int currentMaxMassCount)
  {
    if (currentMaxMassCount <= 0 || peakGroups.size() <= (Size) currentMaxMassCount)
    {
      return;
    }

    Size mc = (Size) currentMaxMassCount;
    std::vector<double> scores;
    scores.reserve(peakGroups.size());
    for (auto &pg : peakGroups)
    {
      scores.push_back(pg.qScore);
    }

    sort(scores.begin(), scores.end());

    auto newPeakGroups = std::vector<PeakGroup>();
    newPeakGroups.reserve(peakGroups.size());
    auto threshold = scores[scores.size() - mc];
    for (auto &pg : peakGroups)
    {
      if (newPeakGroups.size() > mc)
      {
        break;
      }

      if (pg.qScore >= threshold)
      {
        newPeakGroups.push_back(pg);
      }
    }
    std::vector<PeakGroup>().swap(peakGroups);
    newPeakGroups.swap(peakGroups);
  }


  void FLASHDeconvAlgorithm::removeHarmonicPeakGroups(double tol)
  {
    sort(peakGroups.begin(), peakGroups.end());
    std::vector<PeakGroup> merged;
    merged.reserve(peakGroups.size());

    std::vector<double> masses;
    masses.reserve(peakGroups.size());
    for (auto &peakGroup : peakGroups)
    {
      masses.push_back(peakGroup.avgMass);
    }
    for (auto &pg : peakGroups)
    {
      bool select = true;
      for (int h = 2; h <= 3; h++)
      {
        for (int k = 0; k < 2; k++)
        {
          for (int i = -2; i <= 2; ++i)
          {
            auto omass = pg.avgMass + i * Constants::ISOTOPE_MASSDIFF_55K_U;
            auto hmass = k == 0 ? omass * h : omass / h;
            double massTol = 2 * hmass * tol;
            auto iter = std::lower_bound(masses.begin(), masses.end(), hmass - massTol);
            Size j = iter - masses.begin();

            if (j >= 0 && j < peakGroups.size())
            {
              for (; j < peakGroups.size(); j++)
              {
                auto &pgo = peakGroups[j];
                if (hmass - pgo.avgMass > massTol)
                {
                  continue;
                }

                if (!select || pgo.avgMass - hmass > massTol)
                {
                  break;
                }
                select &= pg.totalSNR >= pgo.totalSNR;
                if (!select)
                {
                  break;
                }
              }
              if (!select)
              {
                break;
              }
            }
          }
          if (!select)
          {
            break;
          }
        }
      }
      if (!select)
      {
        continue;
      }
      merged.push_back(pg);
    }
    std::vector<PeakGroup>().swap(peakGroups);
    merged.swap(peakGroups);
  }

  void FLASHDeconvAlgorithm::removeOverlappingPeakGroups(double tol)
  { // pgs are sorted
    int isoLength = 1; // inclusive
    sort(peakGroups.begin(), peakGroups.end());
    std::vector<PeakGroup> filtered;
    filtered.reserve(peakGroups.size());

    for (Size i = 0; i < peakGroups.size(); i++)
    {
      bool select = true;
      auto &pg = peakGroups[i];

      if (pg.avgMass <= 0)
      {
        continue;
      }
      if (i > 0 && abs(pg.avgMass - peakGroups[i - 1].avgMass) < 1e-3)
      {
        continue;
      }

      double massTol = pg.avgMass * tol;

      int j = i + 1;
      for (int l = 0; l <= isoLength; l++)
      {
        auto off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j < peakGroups.size(); j++)
        {
          auto &pgo = peakGroups[j];

          if (l != 0 && pgo.avgMass - pg.avgMass < off - massTol)
          {
            continue;
          }

          if (!select || pgo.avgMass - pg.avgMass > off + massTol)
          {
            break;
          }
          select &= pg.intensity > pgo.intensity;
        }
      }

      if (!select)
      {
        continue;
      }

      j = i - 1;
      for (int l = 0; l <= isoLength; l++)
      {
        auto off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j >= 0; j--)
        {
          auto &pgo = peakGroups[j];

          if (l != 0 && pg.avgMass - pgo.avgMass < off - massTol)
          {
            continue;
          }

          if (!select || pg.avgMass - pgo.avgMass > off + massTol)
          {
            break;
          }
          select &= pg.intensity > pgo.intensity;
        }
      }
      if (!select)
      {
        continue;
      }
      filtered.push_back(pg);
    }
    std::vector<PeakGroup>().swap(peakGroups);
    filtered.swap(peakGroups);
  }


  void FLASHDeconvAlgorithm::reassignPeaksinPeakGroups()
  {
    //auto maxSNR = std::vector<double>(logMzPeaks.size() , 0);
    auto maxIntensity = std::vector<double>(logMzPeaks.size(), 0);
    //auto maxIntensityCharge = std::vector<int>(logMzPeaks.size() , 0);

    for (auto &pg : peakGroups)
    {
      //auto snr = pg.totalSNR;
      auto intensity = pg.intensity / (pg.maxCharge - pg.minCharge + 1);

      for (auto &p: pg)
      {
        auto idx = p.index;
        //if (maxSNR[idx] < snr)
        //{
        //  maxSNR[idx] = snr;
        //}
        if (maxIntensity[idx] < intensity)
        {
          maxIntensity[idx] = intensity;
          //maxIntensityCharge[idx] = p.charge;
        }
      }
    }

    for (auto &pg : peakGroups)
    {
      //auto snr = pg.totalSNR;
      auto intensity = pg.intensity / (pg.maxCharge - pg.minCharge + 1);

      std::vector<LogMzPeak> tmp;
      tmp.swap(pg);
      pg.reserve(tmp.size());
      for (auto &p: tmp)
      {
        auto idx = p.index;
        if (//maxSNR[idx] * .5 > snr
          //||
            maxIntensity[idx] > intensity
          //maxIntensityCharge[idx] != p.charge
            )
        {
          continue;
        }
        pg.push_back(p);
      }
    }
  }


  std::vector<int> FLASHDeconvAlgorithm::updatePerChargeIsotopeIntensity(
      std::vector<double> &perIsotopeIntensity,
      std::vector<double> &perChargeIntensity,
      int maxIsotopeCount,
      PeakGroup &pg)
  {
    auto chargeRange = currentMaxCharge - minCharge + 1;
    std::fill(perIsotopeIntensity.begin(), perIsotopeIntensity.end(), 0);
    std::fill(perChargeIntensity.begin(), perChargeIntensity.end(), 0);

    int minPgCharge = INT_MAX;
    int maxPgCharge = INT_MIN;
    //    double maxIntensity = -1;
    int maxIntChargeIndex = -1;
    //    double maxIntensity2 = -1;
    int maxIntIsoIndex = -1;

    for (auto &p : pg)
    {
      if (p.isotopeIndex < 0 || p.isotopeIndex >= maxIsotopeCount)
      {
        continue;
      }
      minPgCharge = std::min(minPgCharge, p.charge);
      maxPgCharge = std::max(maxPgCharge, p.charge);

      int index = p.charge - minCharge;
      perIsotopeIntensity[p.isotopeIndex] += p.intensity;
      perChargeIntensity[index] += p.intensity;
    }
    pg.maxCharge = maxPgCharge;
    pg.minCharge = minPgCharge;

    return std::vector<int>{maxIntChargeIndex, maxIntIsoIndex};
  }

  double FLASHDeconvAlgorithm::getCosine(std::vector<double> &a, std::vector<double> &b, int off)
  {
    double n = .0, d1 = .0, d2 = .0;
    Size size = a.size();
    //int overlapCntr = 0;
    for (Size j = off; j < size - off; j++)
    {
      d1 += a[j] * a[j];
      d2 += b[j] * b[j];
      n += a[j] * b[j];
      //  if(a[j] > 0 && b[j] > 0) overlapCntr++;
    }

    //if(overlapCntr < 2) return 0; //
    double d = (d1 * d2);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }


  /*double FLASHDeconvAlgorithm::getChargeFitScore(double *perChargeIntensity, int chargeRange)
  {
    double maxPerChargeIntensity = .0;
    std::vector<double> xs;
    std::vector<double> ys;

    xs.reserve(+2);
    ys.reserve(chargeRange + 2);

    for (int i = 0; i < chargeRange; i++)
    {
      maxPerChargeIntensity = std::max(maxPerChargeIntensity, perChargeIntensity[i]);
    }

    double th = maxPerChargeIntensity * .02;// as recommended in the original paper...
    int first = -1, last = 0;
    for (int i = 0; i < chargeRange; i++)
    {
      if (perChargeIntensity[i] <= th)
      {
        continue;
      }
      if (first < 0)
      {
        first = i;
      }

      last = i;
    }

    for (int i = first; i <= last; i++)
    {
      xs.push_back(i);
      ys.push_back((1 + perChargeIntensity[i]));
    }

    if (xs.size() <= 3)
    {
      return 0.5;
    }

    Eigen::Matrix3d m;
    Eigen::Vector3d v;

    double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
    double t0 = 0, t1 = 0, t2 = 0;

    for (Size i = 0; i < xs.size(); i++)
    {
      auto &x = xs[i];
      auto y = log(ys[i]);
      s0++;
      s1 += x;
      s2 += x * x;
      s3 += x * x * x;
      s4 += x * x * x * x;
      t0 += y;
      t1 += y * x;
      t2 += y * x * x;
    }
    m(0, 0) = s0;
    m(1, 0) = m(0, 1) = s1;
    m(2, 0) = m(1, 1) = m(0, 2) = s2;
    m(2, 1) = m(1, 2) = s3;
    m(2, 2) = s4;

    auto im = m.inverse();
    v(0) = t0;
    v(1) = t1;
    v(2) = t2;
    //cout<<v<<endl;
    auto abc = im * v;
    //cout<<abc<<endl;
    double mu = -abc(1) / abc(2) / 2;
    double omega = -1 / abc(2) / 2;

    if (omega <= 0)
    {
      return 0;
    }
    std::vector<double> tys;

    for (Size i = 0; i < ys.size(); i++)
    {
      double ty = exp(-(xs[i] - mu) * (xs[i] - mu) / 2 / omega);
      tys.push_back(ty);
    }
    return getCosine(ys, tys);
  }*/

  double FLASHDeconvAlgorithm::getChargeFitScore(double *perChargeIntensity, int chargeRange)
  {
    double maxPerChargeIntensity = .0;
    std::vector<double> xs;
    std::vector<double> ys;

    xs.reserve(chargeRange + 2);
    ys.reserve(chargeRange + 2);

    for (int i = 0; i < chargeRange; i++)
    {
      maxPerChargeIntensity = std::max(maxPerChargeIntensity, perChargeIntensity[i]);
    }

    /*double th = maxPerChargeIntensity * .02;// as recommended in the original paper...
    for (int i = 0; i < chargeRange; i++)
    {
      if (perChargeIntensity[i] <= th)
      {
        continue;
      }

      xs.push_back(i);
      ys.push_back((1 + perChargeIntensity[i]));
    }*/
    double th = maxPerChargeIntensity * .02;// as recommended in the original paper...
    int first = -1, last = 0;
    for (int i = 0; i < chargeRange; i++)
    {
      if (perChargeIntensity[i] <= th)
      {
        continue;
      }
      if (first < 0)
      {
        first = i;
      }

      last = i;
    }

    for (int i = first; i <= last; i++)
    {
      if (perChargeIntensity[i] <= 0)
      {
        continue;
      }
      xs.push_back(i);
      ys.push_back((perChargeIntensity[i]));
    }

    if (xs.size() <= 3)
    {
      return 0.5;
    }

    Eigen::Matrix3d m;
    Eigen::Vector3d v;

    double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
    double t0 = 0, t1 = 0, t2 = 0;

    for (Size i = 0; i < xs.size(); i++)
    {
      auto &x = xs[i];
      auto y = log(ys[i]);
      s0++;
      s1 += x;
      s2 += x * x;
      s3 += x * x * x;
      s4 += x * x * x * x;
      t0 += y;
      t1 += y * x;
      t2 += y * x * x;
    }
    m(0, 0) = s0;
    m(1, 0) = m(0, 1) = s1;
    m(2, 0) = m(1, 1) = m(0, 2) = s2;
    m(2, 1) = m(1, 2) = s3;
    m(2, 2) = s4;

    auto im = m.inverse();
    v(0) = t0;
    v(1) = t1;
    v(2) = t2;
    //cout<<v<<endl;
    auto abc = im * v;
    //cout<<abc<<endl;
    double mu = -abc(1) / abc(2) / 2;
    double omega = -1 / abc(2) / 2;

    if (omega <= 0)
    {
      return -.1;
    }
    std::vector<double> tys;

    for (Size i = 0; i < ys.size(); i++)
    {
      double ty = exp(-(xs[i] - mu) * (xs[i] - mu) / 2 / omega);
      tys.push_back(ty);
    }
    return getCosine(ys, tys);
  }

}