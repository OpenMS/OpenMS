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
    auto chargeRange = maxCharge - minCharge;
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
        if (abs(i + minCharge) == 1)
        {
          factor = 1;
        }
        harmonicFilter.setValue(k, i, log(1.0 / (factor * n / hc + abs(i + minCharge))));
      }
    }
  }

  //Generate uncharged log mz transformated peaks
  void FLASHDeconvAlgorithm::updateLogMzPeaks(MSSpectrum *spec)
  {
    std::vector<LogMzPeak>().swap(logMzPeaks);
    logMzPeaks.reserve(spec->size());
    for (auto &peak: *spec)
    {
      if (peak.getIntensity() <= intensityThreshold)//
      {
        continue;
      }
      LogMzPeak logMzPeak(peak, minCharge > 0);
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
                                          float *mzBinIntensities)
  {
    mzBins = boost::dynamic_bitset<>(binNumber);
    std::fill_n(mzBinIntensities, binNumber, .0);

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
        if (bi < binNumber - 1 && !mzBins[bi + 1])
        {
          mzBins.set(bi + 1);
          mzBinIntensities[bi + 1] += p.intensity;
        }
      }
      else if (delta < 0)
      {
        if (bi > 0 && !mzBins[bi - 1])
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
  boost::dynamic_bitset<> FLASHDeconvAlgorithm::getCandidateMassBinsForThisSpectrum(float *massIntensitites,
                                                                                    float *mzIntensities)
  {
    int chargeRange = currentMaxCharge - minCharge;
    int hChargeSize = (int) hCharges.size();
    int minPeakCntr = minSupportPeakCount[msLevel - 1];
    long binEnd = (long) massBins.size();
    auto candidateMassBinsForThisSpectrum = boost::dynamic_bitset<>(massBins.size());
    // how many peaks of continuous charges per mass
    int *continuousChargePeakPairCount = new int[massBins.size()];
    std::fill_n(continuousChargePeakPairCount, massBins.size(), 0);

    auto mzBinIndex = mzBins.find_first();
    std::fill_n(massIntensitites, massBins.size(), 0);

    // to calculate continuous charges, the previous charge value per mass should be stored
    int *prevCharges = new int[massBins.size()];
    std::fill_n(prevCharges, massBins.size(), (chargeRange + 2));

    // not just charges but intensities are stored to see the intensity fold change
    auto *prevIntensities = new float[massBins.size()];
    std::fill_n(prevIntensities, massBins.size(), 1);

    // harmonic peak intensities contribute to noise
    auto noise = new float *[hChargeSize + 1];
    for (auto k = 0; k < hChargeSize + 1; k++)
    {
      noise[k] = new float[massBins.size()];
      std::fill_n(noise[k], massBins.size(), 0);
    }
    // bin width per ms level
    auto bw = binWidth[msLevel - 1];

    // intensity change ratio should not exceed the factor. factor2 just to reduce run time
    float factor = 4.0;
    float factor2 = (1 - 1.0f / factor);
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
        else if (prevCharge < chargeRange && chargeNotContinous && intensityRatio < factor)
        {
          noise[hChargeSize][massBinIndex] += intensity;
        }

        // if charge not continous or intensity ratio is too high reset continuousChargePeakPairCount
        if (chargeNotContinous || intensityRatio > factor)
        {
          continuousChargePeakPairCount[massBinIndex] = 0;
        }
        else
        {
          passFirstStep = true;
        }

        if (passFirstStep)
        {
          //check harmonic
          int maxHcharge = -1;
          float maxHint = .0;
          auto highThreshold = intensity * factor;
          auto lowThreshold = intensity * factor2;// / factor;
          for (auto k = 0; k < hChargeSize; k++)
          {
            auto hmzBinIndex = massBinIndex - hBinOffsets.getValue(k, j);
            if (hmzBinIndex > 0 && hmzBinIndex < (long) mzBins.size() && mzBins[hmzBinIndex])
            {
              auto &hintensity = mzIntensities[hmzBinIndex];
              if (hintensity > lowThreshold
                  &&
                  hintensity < highThreshold)
              {
                if (hintensity < maxHint)
                {
                  continue;
                }
                maxHint = hintensity;
                maxHcharge = k;
              }
            }
          }

          auto &cc = continuousChargePeakPairCount[massBinIndex];
          if (maxHcharge >= 0) //
          {
            noise[maxHcharge][massBinIndex] += maxHint;
            cc = 0;
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

    delete[] prevIntensities;
    delete[] prevCharges;
    for (auto k = 0; k <= hChargeSize; k++)
    {
      delete[] noise[k];
    }
    delete[] noise;
    delete[] continuousChargePeakPairCount;
    return candidateMassBinsForThisSpectrum;
  }


  // Subfunction of updateMassBins. If a peak corresponds to multiple masses, only one mass is selected baed on intensities..
  Matrix<int> FLASHDeconvAlgorithm::updateMassBins_(boost::dynamic_bitset<> &candidateMassBinsForThisSpectrum,
                                                    float *massIntensities,
                                                    long &binStart, long &binEnd)
  {
    //int chargeRange = param.currentChargeRange;
    auto chargeRange = currentMaxCharge - minCharge;
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
  Matrix<int> FLASHDeconvAlgorithm::updateMassBins(float *massIntensities,
                                                   float *mzIntensities)
  {
    auto bw = binWidth[msLevel - 1];
    long binThresholdMinMass = (long) getBinNumber(log(minMass), massBinMinValue, bw);
    long binThresholdMaxMass = (long) std::min(massBins.size(),
                                               1 + getBinNumber(log(currentMaxMass),
                                                                massBinMinValue,
                                                                bw));

    auto candidateMassBins = getCandidateMassBinsForThisSpectrum(massIntensities,
                                                                 mzIntensities);

    auto perMassChargeRanges = updateMassBins_(candidateMassBins,
                                               massIntensities,
                                               binThresholdMinMass,
                                               binThresholdMaxMass);

    return perMassChargeRanges;
  }

  //With massBins, select peaks from the same mass in the original input spectrum
  void FLASHDeconvAlgorithm::getCandidatePeakGroups(float *massIntensities,
                                                    Matrix<int> &chargeRanges)
  {
    auto maxMissingIsotope = 2;//msLevel == 1 ? 2 : 1;
    double bw = binWidth[msLevel - 1];
    double tol = tolerance[msLevel - 1];
    int chargeRange = currentMaxCharge - minCharge;
    auto mzBinSize = mzBins.size();
    auto massBinSize = massBins.size();
    int logMzPeakSize = (int) logMzPeaks.size();
    int *currentPeakIndex = new int[chargeRange];
    std::fill_n(currentPeakIndex, chargeRange, 0);
    peakGroups.reserve(massBins.count());
    auto massBinIndex = massBins.find_first();
    Size *peakBinNumbers = new Size[logMzPeakSize];

    for (int i = 0; i < logMzPeakSize; i++)
    {
      peakBinNumbers[i] = getBinNumber(logMzPeaks[i].logMz, mzBinMinValue, bw);
    }

    while (massBinIndex != massBins.npos)
    {
      double logM = getBinValue(massBinIndex, massBinMinValue, bw);
      double mass = exp(logM);
      if (msLevel == 1)
      {
        double diff = Constants::ISOTOPE_MASSDIFF_55K_U / mass;
        double isoLogM1 = logM - diff;
        double isoLogM2 = logM + diff;

        auto b1 = getBinNumber(isoLogM1, massBinMinValue, bw);

        if (b1 > 0 && b1 < massBinIndex)
        {
          if (massIntensities[massBinIndex] < massIntensities[b1])
          {
            massBinIndex = massBins.find_next(massBinIndex);
            continue;
          }
        }
        auto b2 = getBinNumber(isoLogM2, massBinMinValue, bw);

        if (b2 < massBinSize && b2 > massBinIndex)
        {
          if (massIntensities[massBinIndex] < massIntensities[b2])
          {
            massBinIndex = massBins.find_next(massBinIndex);
            continue;
          }
        }

        if (massIntensities[b1] == 0 && massIntensities[b2] == 0) //
        {
          massBinIndex = massBins.find_next(massBinIndex);
          continue;
        }
      }

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

        pg.massBinIndex = massBinIndex;
        peakGroups.push_back(pg); //
      }
      massBinIndex = massBins.find_next(massBinIndex);
    }
    delete[] currentPeakIndex;
    delete[] peakBinNumbers;
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
    auto currentChargeRange = currentMaxCharge - minCharge;
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
    auto *mzBinIntensities = new float[mzBinNumber];

    updateMzBins(mzBinNumber, mzBinIntensities);
    auto *massIntensities = new float[massBinNumber];

    massBins = boost::dynamic_bitset<>(massBinNumber);
    massBinsForThisSpectrum = boost::dynamic_bitset<>(massBinNumber);

    if (msLevel == 1)
    {
      unionPrevMassBins();
    }
    auto perMassChargeRanges = updateMassBins(massIntensities, mzBinIntensities);

    getCandidatePeakGroups(massIntensities, perMassChargeRanges);

    PeakGroupScoring scorer = PeakGroupScoring(peakGroups,
                                               minCharge,
                                               currentMaxCharge,
                                               minChargeCosine,
                                               tolerance,
                                               maxMassCount,
                                               minSupportPeakCount,
                                               minIsotopeCosine);

    peakGroups = scorer.scoreAndFilterPeakGroups(msLevel, avg);

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
        if (massBinsForThisSpectrum[pg.massBinIndex])
        {
          mb.push_back(pg.massBinIndex);
        }
      }

      prevMassBinVector.push_back(mb); //
      prevMinBinLogMassVector.push_back(massBinMinValue);
      prevMassBinVector.shrink_to_fit();
      prevMinBinLogMassVector.shrink_to_fit();
    }

    delete[] mzBinIntensities;
    delete[] massIntensities;
  }


}