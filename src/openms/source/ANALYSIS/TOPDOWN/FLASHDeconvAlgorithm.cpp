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
    prev_mass_bin_vector = std::vector<std::vector<Size>>();
    prev_minbin_logmass_vector = std::vector<double>();
    defaults_.setValue("tol",
                       DoubleList{10.0, 5.0},
                       "ppm tolerance for MS1, 2, ... (e.g., -tol 10.0 5.0 to specify 10.0 and 5.0 ppm for MS1 and MS2, respectively)");

    defaults_.setValue("min_mass", 50.0, "minimum mass (Da)");
    defaults_.setValue("max_mass", 100000.0, "maximum mass (Da)");

    defaults_.setValue("min_charge", 1, "minimum charge state (can be negative for negative mode)");
    defaults_.setValue("max_charge", 100, "maximum charge state (can be negative for negative mode)");

    defaults_.setValue("min_mz", -1.0, "if set to positive value, minimum m/z to deconvolute.");
    defaults_.setValue("max_mz", -1.0, "if set to positive value, maximum m/z to deconvolute.");
    defaults_.setValue("min_rt", -1.0, "if set to positive value, minimum RT to deconvolute.");
    defaults_.setValue("max_rt", -1.0, "if set to positive value, maximum RT to deconvolute.");

    defaults_.setValue("min_isotope_cosine",
                       DoubleList{.75, .75},
                       "cosine threshold between avg. and observed isotope pattern for MS1, 2, ... (e.g., -min_isotope_cosine 0.8 0.6 to specify 0.8 and 0.6 for MS1 and MS2, respectively)");

    defaults_.setValue("min_peaks",
                       IntList{3, 1},
                       "minimum number of supporting peaks for MS1, 2, ...  (e.g., -min_peaks 3 2 to specify 3 and 2 for MS1 and MS2, respectively)");

    //defaults_.setValue("min_charge_score",
    //                   DoubleList{.0, .0},
    //                   "charge score threshold for MS1, 2, ... (e.g., -min_charge_score 0.7 0.3 to specify 0.7 and 0.3 for MS1 and MS2, respectively)");

    //defaults_.setValue("min_charge_cosine",
    //                   .5,
    //                   "cosine threshold between per-charge-intensity and fitted gaussian distribution (applies only to MS1)");
    defaults_.setValue("max_mass_count",
                       IntList{-1, -1},
                       "maximum mass count per spec for MS1, 2, ... (e.g., -max_mass_count 100 50 to specify 100 and 50 for MS1 and MS2, respectively. -1 specifies unlimited)");

    defaults_.setValue("min_mass_count",
                       IntList{-1, -1},
                       "minimum mass count per spec for MS1, 2, ... "
                       "this parameter is only for real time acquisition. "
                       "the parameter may not be satisfied in case spectrum quality is too poor. (e.g., -max_mass_count -1 2 to specify no min limit and 2 for MS1 and MS2, respectively. -1 specifies unlimited)");
    defaults_.addTag("min_mass_count", "advanced");

    defaults_.setValue("min_intensity", .0, "intensity threshold");
    defaults_.setValue("RT_window", 20.0, "RT window for MS1 deconvolution");
    defaultsToParam_();
  }

  FLASHDeconvAlgorithm::~FLASHDeconvAlgorithm()
  {
    if (log_mz_peaks.empty())
    {
      return;
    }
    std::vector<LogMzPeak>().swap(log_mz_peaks);
  }


  FLASHDeconvAlgorithm& FLASHDeconvAlgorithm::operator=(const FLASHDeconvAlgorithm& fd)
  {
    if (this ==& fd)
    {
      return *this;
    }
    //...
    return *this;
  }

  ///Calcualte the nominla mass from double mass. Mutiply 0.999497 reduces the rounding error.
  int FLASHDeconvAlgorithm::getNominalMass(const double mass)
  {
    return (int) (mass * 0.999497 + .5);
  }


  //This function is the main function for the deconvolution. Takes empty DeconvolutedSpectrum and fill it up with peakGroups.
  // DeconvolutedSpectrum contains the recursor peak group for MSn.
  //A peakGroup is the collection of peaks from a single mass (monoisotopic mass). Thus it contains peaks from different charges and iostope indices.
  void FLASHDeconvAlgorithm::fillPeakGroupsInDeconvolutedSpectrum(DeconvolutedSpectrum& deconvoluted_spec, const int scan_number)
  {
    auto *spec =& (deconvoluted_spec.getOriginalSpectrum());
    deconvoluted_spectrum =& deconvoluted_spec;
    //std::vector<PeakGroup> _peakGroups;
    if (min_rt > 0 && spec->getRT() < min_rt)
    {
      return;
    }
    if (max_rt > 0 && spec->getRT() > max_rt)
    {
      return;
    }

    updateLogMzPeaks(spec);
    ms_level = spec->getMSLevel();
    //For MS2,3,.. max mass and charge ranges should be determined by precursors
    current_max_charge = deconvoluted_spectrum->getCurrentMaxCharge(max_charge); // TODO negative mode!!
    current_max_mass = deconvoluted_spectrum->getCurrentMaxMass(max_mass);
    //if(msLevel>1) { std::cout<<std::to_string(current_max_mass)<<" " << current_max_charge << std::endl; }
    //Perform deconvolution and fill in deconvoluted_spectrum
    generatePeakGroupsFromSpectrum();
    if (deconvoluted_spectrum->empty())
    {
      return;
    }
    //Update peakGroup information after deconvolution
    for (auto& pg : *deconvoluted_spectrum)
    {
      sort(pg.begin(), pg.end());
      pg.setScanNumber(scan_number);
    }
  }

  void FLASHDeconvAlgorithm::updateMembers_()
  {
    min_mz = param_.getValue("min_mz");
    max_mz = param_.getValue("max_mz");

    min_rt = param_.getValue("min_rt");
    max_rt = param_.getValue("max_rt");

    min_charge = param_.getValue("min_charge");
    max_charge = param_.getValue("max_charge");

    if (min_charge > max_charge)
    {
      int tmp = min_charge;
      min_charge = max_charge;
      max_charge = tmp;
    }

    max_mass = param_.getValue("max_mass");
    min_mass = param_.getValue("min_mass");

    intensity_threshold = param_.getValue("min_intensity");
    min_support_peak_count = param_.getValue("min_peaks");

    bin_width.clear();
    tolerance = param_.getValue("tol");

    for (int j = 0; j < (int) tolerance.size(); j++)
    {
      tolerance[j] *= 1e-6;
      bin_width.push_back(.5 / tolerance[j]);
    }

    min_isotope_cosine = param_.getValue("min_isotope_cosine");
    //minChargeScore = param_.getValue("min_charge_score");
    //minChargeCosine = param_.getValue("min_charge_cosine");
    max_mass_count = param_.getValue("max_mass_count");
    min_mass_count = param_.getValue("min_mass_count");
    rt_window = param_.getValue("RT_window");
    setFilters();
  }

  FLASHDeconvHelperStructs::PrecalculatedAveragine FLASHDeconvAlgorithm::getAveragine()
  {
    return avg;
  }

  void FLASHDeconvAlgorithm::calculateAveragine(const bool use_RNA_averagine)
  {
    avg = FLASHDeconvHelperStructs::calculateAveragines(max_mass, use_RNA_averagine);
  }


  // generate filters
  void FLASHDeconvAlgorithm::setFilters()
  {
    filter.clear();
    harmonic_filter_matrix.clear();
    int chargeRange = max_charge - min_charge + 1;
    for (int i = 0; i < chargeRange; i++)
    {
      filter.push_back(log(1.0 / abs(i + min_charge)));
    }

    harmonic_filter_matrix.resize(harmonic_charges.size(), chargeRange);

    for (Size k = 0; k < harmonic_charges.size(); k++)
    {
      int hc = harmonic_charges[k];
      int n = hc / 2;

      for (int i = 0; i < chargeRange; i++)
      {
        double factor = 1;
        //if (abs(i + minCharge) == 1)
        // {
        //  factor = 1;
        //}
        //auto harmonic_filter_matrix = log(1.0 / (i + n / hc + param.minCharge));
        //        harmonic_bin_offset_matrix[i][k] = (long) round((filter[i] - harmonic_filter_matrix) * param.bin_width);

        harmonic_filter_matrix.setValue(k, i, log(1.0 / (factor * n / hc + abs(i + min_charge))));
      }
    }
  }

  //Generate uncharged log mz transformated peaks
  void FLASHDeconvAlgorithm::updateLogMzPeaks(const MSSpectrum *spec)
  {
    std::vector<LogMzPeak>().swap(log_mz_peaks);
    log_mz_peaks.reserve(spec->size());
    //int index = 0;
    for (auto& peak : *spec)
    {
      if (min_mz > 0 && peak.getMZ() < min_mz)
      {
        continue;
      }
      if (max_mz > 0 && peak.getMZ() > max_mz)
      {
        break;
      }
      if (peak.getIntensity() <= intensity_threshold)//
      {
        continue;
      }
      LogMzPeak logMzPeak(peak, min_charge > 0);
      log_mz_peaks.push_back(logMzPeak);
    }
  }


  double FLASHDeconvAlgorithm::getBinValue(const Size bin, const double min_value, const double bin_width)
  {
    return min_value + bin / bin_width;
  }

  Size FLASHDeconvAlgorithm::getBinNumber(const double value, const double min_value, const double bin_width)
  {
    if (value < min_value)
    {
      return 0;
    }
    return (Size) (((value - min_value) * bin_width) + .5);

  }

  // From log mz to mz bins. To reduce edge effect, its adjacent bin is set when a bin is set.
  void FLASHDeconvAlgorithm::updateMzBins(const Size& bin_number,
                                          std::vector<float>& mz_bin_intensities)
  {
    mz_bins_for_edge_effect = boost::dynamic_bitset<>(bin_number);
    mz_bins = boost::dynamic_bitset<>(bin_number);
    //std::fill(mz_bin_intensities.begin(), mz_bin_intensities.end(), .0);

    for (auto& p : log_mz_peaks)
    {
      Size bi = getBinNumber(p.logMz, mz_bin_min_value, bin_width[ms_level - 1]);
      if (bi >= bin_number)
      {
        continue;
      }
      mz_bins.set(bi);
      mz_bins_for_edge_effect.set(bi);
      mz_bin_intensities[bi] += p.intensity;
    }
    for (auto& p : log_mz_peaks)
    {
      Size bi = getBinNumber(p.logMz, mz_bin_min_value, bin_width[ms_level - 1]);
      double delta = (p.logMz - getBinValue(bi, mz_bin_min_value, bin_width[ms_level - 1]));

      if (delta > 0)
      {
        if (bi < bin_number - 1
            && !mz_bins_for_edge_effect[bi + 1]
            )
        {
          mz_bins_for_edge_effect.set(bi + 1);
          mz_bin_intensities[bi + 1] += p.intensity;
        }
      }
      else if (delta < 0)
      {
        if (bi > 0
            && !mz_bins_for_edge_effect[bi - 1]
            )
        {
          mz_bins_for_edge_effect.set(bi - 1);
          mz_bin_intensities[bi - 1] += p.intensity;
        }
      }
    }
  }

  //take the mass bins from previous overlapping spectra and put them in the candidate mass bins.
  void FLASHDeconvAlgorithm::unionPrevMassBins()
  {
    if (mass_bins.empty())
    {
      return;
    }
    for (Size i = 0; i < prev_mass_bin_vector.size(); i++)
    {
      auto& pmb = prev_mass_bin_vector[i];
      if (pmb.empty())
      {
        continue;
      }
      long shift = (long) (round((mass_bin_min_value - prev_minbin_logmass_vector[i]) * bin_width[ms_level - 1]));

      for (Size& index : pmb)
      {
        long j = (long) index - shift;
        if (j < 0)
        {
          continue;
        }
        if ((Size) j >= mass_bins.size())
        {
          break;
        }
        mass_bins[j] = true;
      }
    }
  }

  //Find candidate mass bins from the current spectrum. The runtime of FLASHDeconv is deteremined by this function..
  void FLASHDeconvAlgorithm::updateCandidateMassBins(std::vector<float>& mass_intensitites,
                                                                        const std::vector<float>& mz_intensities)
  {
    int chargeRange = current_max_charge - min_charge + 1;
    int hChargeSize = (int) harmonic_charges.size();
    int minPeakCntr = min_support_peak_count[ms_level - 1];
    long binEnd = (long) mass_bins.size();
    //auto candidateMassBinsForThisSpectrum = boost::dynamic_bitset<>(mass_bins.size());
    // how many peaks of continuous charges per mass
    auto supportPeakCount = std::vector<int>(mass_bins.size(), 0);

    Size mzBinIndex = mz_bins.find_first();
    //std::fill(mass_intensitites.begin(), mass_intensitites.end(), 0);

    // to calculate continuous charges, the previous charge value per mass should be stored
    auto prevCharges = std::vector<int>(mass_bins.size(), chargeRange + 2);

    // not just charges but intensities are stored to see the intensity fold change
    auto prevIntensities = std::vector<float>(mass_bins.size(), 1.0f);

    double bw = bin_width[ms_level - 1];

    // intensity change ratio should not exceed the factor.
    float factor = 10.0;

    while (mzBinIndex != mz_bins.npos)
    {
      float intensity = mz_intensities[mzBinIndex];
      double mz = -1.0, logMz = 0;
      if (ms_level > 1)
      {
        logMz = getBinValue(mzBinIndex, mz_bin_min_value, bw);
        mz = exp(logMz);
      }

      // scan through charges
      for (int j = 0; j < chargeRange; j++)
      {
        // mass is given by shifting by bin_offsets[j]
        long massBinIndex = mzBinIndex + bin_offsets[j];

        if (massBinIndex < 0)
        {
          continue;
        }


        if (massBinIndex >= binEnd)
        {
          break;
        }

        auto& spc = supportPeakCount[massBinIndex];

        if (ms_level == 1)
        {
          // intensity of previous charge
          float& prevIntensity = prevIntensities[massBinIndex];
          // intensity ratio between current and previous charges
          float intensityRatio = intensity / prevIntensity;
          intensityRatio = intensityRatio < 1 ? 1.0f / intensityRatio : intensityRatio;

          // check if peaks of continuous charges are present
          int& prevCharge = prevCharges[massBinIndex];
          bool chargeNotContinous = prevCharge - j != 1;
          // if charge not continous or intensity ratio is too high reset continuousChargePeakPairCount
          if (chargeNotContinous || intensityRatio > factor)
          {
            spc = 0;
          }
          else
          { // check harmonic artifacts
            float highThreshold = intensity * factor;
            float lowThreshold = intensity / factor;// / factor;
            bool isHarmonic = false;
            for (int k = 0; k < hChargeSize; k++)//
            {
              int hmzBinIndex = massBinIndex - harmonic_bin_offset_matrix.getValue(k, j);
              if (hmzBinIndex > 0 && hmzBinIndex < mz_bins_for_edge_effect.size() && mz_bins_for_edge_effect[hmzBinIndex])
              {
                float hintensity = mz_intensities[hmzBinIndex];
                if (hintensity > lowThreshold
                    &&
                    hintensity < highThreshold
                    )
                {
                  isHarmonic = true;
                  break;//
                  //maxHcharge = k;
                }
              }
            }
            if (!isHarmonic)
            {
              if (spc == 0)
              {
                mass_intensitites[massBinIndex] += prevIntensity;
              }

              mass_intensitites[massBinIndex] += intensity;
              mass_bins[massBinIndex] = (++spc >= minPeakCntr);
            }
          }
          prevIntensity = intensity;
          prevCharge = j;
        }
        else // for MS2,3,... iostopic peaks or water nh3 loss peaks are considered
        {
          int acharge = abs(j + min_charge); // abs charge
          bool supportPeakPresent = false;
          double isoIntensity = .0;
          double diff = Constants::ISOTOPE_MASSDIFF_55K_U / acharge / mz;
          Size nextIsoBin = getBinNumber(logMz + diff, mz_bin_min_value, bw);

          if (nextIsoBin < mz_bins_for_edge_effect.size() && mz_bins_for_edge_effect[nextIsoBin])
          {
            isoIntensity = mz_intensities[nextIsoBin];
            supportPeakPresent = true;
            //spc++;
          }

          if (supportPeakPresent)
          {
            double waterLossMz = log(mz - 18.010565 / acharge); // 17.026549
            Size waterLossBin = getBinNumber(waterLossMz, mz_bin_min_value, bw);

            if (waterLossBin >= 0 && mz_bins_for_edge_effect[waterLossBin])
            {
              float waterLossIntensity = mz_intensities[waterLossBin];
              if (waterLossIntensity < intensity)
              {
                isoIntensity += waterLossIntensity;
              }
            }

            double amoniaLossMz = log(mz - 17.026549 / acharge); // 17.026549
            Size amoniaLossBin = getBinNumber(amoniaLossMz, mz_bin_min_value, bw);

            if (amoniaLossBin >= 0 && mz_bins_for_edge_effect[amoniaLossBin])
            {
              float amoniaLossntensity = mz_intensities[amoniaLossBin];
              if (amoniaLossntensity < intensity)
              {
                isoIntensity += amoniaLossntensity;
              }
            }

            mass_intensitites[massBinIndex] += intensity + isoIntensity;
            mass_bins[massBinIndex] = (++spc >= minPeakCntr);
          }
          else
          {
            mass_intensitites[massBinIndex] -= intensity;
            //spc = 0;
          }
        }
      }
      mzBinIndex = mz_bins.find_next(mzBinIndex);
    }
  }

  // Subfunction of updateMassBins. If a peak corresponds to multiple masses, only one mass is selected based on intensities..
  Matrix<int> FLASHDeconvAlgorithm::filterMassBins(const std::vector<float>& massIntensities)
  {
    //int chargeRange = param.currentChargeRange;
    int chargeRange = current_max_charge - min_charge + 1;
    Matrix<int> chargeRanges(2, mass_bins.size(), INT_MAX);
    for (int i = 0; i < mass_bins.size(); i++)
    {
      chargeRanges.setValue(1, i, INT_MIN);
    }
    Size mzBinIndex = mz_bins.find_first();
    long binSize = (long) mass_bins.size();

    //massBinsForThisSpectrum = boost::dynamic_bitset<>(mass_bins.size());
    auto toSkip = mass_bins.flip();
    mass_bins.reset();

    while (mzBinIndex != mz_bins.npos)
    {
      long maxIndex = -1;
      float maxCount = -1e11;
      int charge = 0;

      for (int j = 0; j < chargeRange; j++)
      {
        long massBinIndex = mzBinIndex + bin_offsets[j];

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

        float t = massIntensities[massBinIndex];

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

      if (maxIndex >= 0 && maxIndex < binSize)
      {
        chargeRanges.setValue(0, maxIndex, std::min(chargeRanges.getValue(0, maxIndex), charge));
        chargeRanges.setValue(1, maxIndex, std::max(chargeRanges.getValue(1, maxIndex), charge));
        //massBinsForThisSpectrum[maxIndex] = candidateMassBinsForThisSpectrum[maxIndex];
        mass_bins[maxIndex] = true;
      }
      mzBinIndex = mz_bins.find_next(mzBinIndex);
    }
    return chargeRanges;
  }

  //update mass bins which will be used to select peaks in the input spectrum...
  Matrix<int> FLASHDeconvAlgorithm::updateMassBins(//float *massIntensities,
      const std::vector<float>& mz_intensities)
  {
    auto mass_intensities = std::vector<float>(mass_bins.size(), 0);
    updateCandidateMassBins(mass_intensities, mz_intensities);
    auto per_mass_charge_ranges = filterMassBins(mass_intensities);

    return per_mass_charge_ranges;
  }

  //With mass_bins, select peaks from the same mass in the original input spectrum
  void FLASHDeconvAlgorithm::getCandidatePeakGroups(const Matrix<int>& charge_ranges)
  {
    const int maxMissingIsotope = 2;
    double bw = bin_width[ms_level - 1];
    double tol = tolerance[ms_level - 1];
    int chargeRange = current_max_charge - min_charge + 1;
    Size massBinSize = mass_bins.size();
    int logMzPeakSize = (int) log_mz_peaks.size();
    auto currentPeakIndex = std::vector<int>(chargeRange, 0);
    deconvoluted_spectrum->reserve(mass_bins.count());
    Size massBinIndex = mass_bins.find_first();
    auto peakBinNumbers = std::vector<Size>(logMzPeakSize);

    for (int i = 0; i < logMzPeakSize; i++)
    {
      peakBinNumbers[i] = getBinNumber(log_mz_peaks[i].logMz, mz_bin_min_value, bw);
    }

    while (massBinIndex != mass_bins.npos)
    {
      double logM = getBinValue(massBinIndex, mass_bin_min_value, bw);
      double mass = exp(logM);
      PeakGroup pg(min_charge, charge_ranges.getValue(1, massBinIndex) + min_charge);

      pg.reserve(chargeRange * 30);
      Size rightIndex = avg.getIsotopeEndIndex(mass);
      Size leftIndex = avg.getIsotopeStartIndex(mass);

      for (int j = charge_ranges.getValue(0, massBinIndex); j <= charge_ranges.getValue(1, massBinIndex); j++)
      {
        int& binOffset = bin_offsets[j];
        int bi = massBinIndex - binOffset;

        double maxIntensity = -1.0;
        int charge = j + min_charge;
        int& cpi = currentPeakIndex[j];
        int maxPeakIndex = -1;

        while (cpi < logMzPeakSize - 1)
        {
          if (peakBinNumbers[cpi] == bi)
          {
            double intensity = log_mz_peaks[cpi].intensity;
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

        const double mz = log_mz_peaks[maxPeakIndex].mz;
        const double isof = Constants::ISOTOPE_MASSDIFF_55K_U / abs(charge);
        double mzDelta = tol * mz * 2; //

        double np = .0;

        int pi = 0;
        // int peakcntr = 0;
        for (int peakIndex = maxPeakIndex; peakIndex < logMzPeakSize; peakIndex++)
        {
          const double observedMz = log_mz_peaks[peakIndex].mz;
          const double intensity = log_mz_peaks[peakIndex].intensity;
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
            const Size bin = peakBinNumbers[peakIndex] + binOffset;
            if (bin < massBinSize)
            {
              LogMzPeak p(log_mz_peaks[peakIndex]);
              p.charge = charge;
              pg.push_back(p);
            }
            pi = i;
          }
        }

        pi = 0;

        for (int peakIndex = maxPeakIndex - 1; peakIndex >= 0; peakIndex--)
        {
          const double observedMz = log_mz_peaks[peakIndex].mz;
          const double intensity = log_mz_peaks[peakIndex].intensity;

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
            const Size bin = peakBinNumbers[peakIndex] + binOffset;

            if (bin < massBinSize)
            {
              LogMzPeak p(log_mz_peaks[peakIndex]);
              p.charge = charge;
              pg.push_back(p);
            }

            pi = i;
          }
        }
        if (np > 0)
        {
          pg.setChargeSNR(charge, np);
        }
      }

      if (!pg.empty())
      {
        double maxIntensity = -1.0;
        double tmaxMass = .0;
        auto newPeaks = std::vector<LogMzPeak>();
        newPeaks.reserve(pg.size());
        for (auto& p : pg)
        {
          if (maxIntensity < p.intensity)
          {
            maxIntensity = p.intensity;
            tmaxMass = p.getUnchargedMass();
          }
        }
        double isoDelta = tol * tmaxMass;
        int minOff = 10000;
        for (auto& p : pg)
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
        //std::vector<LogMzPeak>().swap(newPeaks);

        for (auto& p : pg)
        {
          p.isotopeIndex -= minOff;
        }
        pg.updateMassesAndIntensity();
        deconvoluted_spectrum->push_back(pg); //
      }
      massBinIndex = mass_bins.find_next(massBinIndex);
    }
  }

  bool FLASHDeconvAlgorithm::empty()
  {
    return log_mz_peaks.empty();
  }

  //spectral deconvolution main function
  void FLASHDeconvAlgorithm::generatePeakGroupsFromSpectrum()
  {
    std::vector<PeakGroup> empty;
    deconvoluted_spectrum->swap(empty);
    int minPeakCntr =
        min_support_peak_count[ms_level - 1];
    int currentChargeRange = current_max_charge - min_charge + 1;
    int tmp = currentChargeRange - minPeakCntr;

    tmp = tmp < 0 ? 0 : tmp;
    double massBinMaxValue = std::min(
        log_mz_peaks[log_mz_peaks.size() - 1].logMz -
        filter[tmp],
        log(current_max_mass + avg.getAverageMassDelta(current_max_mass) + 1));

    double bw = bin_width[ms_level - 1];
    tmp = minPeakCntr - 1;
    tmp = tmp < 0 ? 0 : tmp;
    //filter.push_back(log(1.0 / abs(i + minCharge)));

    mass_bin_min_value = log(std::max(1.0, min_mass - avg.getAverageMassDelta(min_mass)));
    mz_bin_min_value = log_mz_peaks[0].logMz;

    double mzBinMaxValue = log_mz_peaks[log_mz_peaks.size() - 1].logMz;
    Size massBinNumber = getBinNumber(massBinMaxValue, mass_bin_min_value, bw) + 1;

    for (int i = 0; i < currentChargeRange; i++)
    {
      bin_offsets.push_back((int) round((mz_bin_min_value - filter[i] - mass_bin_min_value) * bw));
    }

    //long massBinIndex = mzBinIndex + bin_offsets[j];

    harmonic_bin_offset_matrix.resize(harmonic_charges.size(), currentChargeRange);
    for (Size k = 0; k < harmonic_charges.size(); k++)
    {
      std::vector<int> _hBinOffsets;
      for (int i = 0; i < currentChargeRange; i++)
      {
        harmonic_bin_offset_matrix
            .setValue(k, i, (int) round((mz_bin_min_value - harmonic_filter_matrix.getValue(k, i) - mass_bin_min_value) * bw));
      }
    }

    Size mzBinNumber = getBinNumber(mzBinMaxValue, mz_bin_min_value, bw) + 1;
    auto mzBinIntensities = std::vector<float>(mzBinNumber, .0f);

    updateMzBins(mzBinNumber, mzBinIntensities);

    mass_bins = boost::dynamic_bitset<>(massBinNumber);
    //massBinsForThisSpectrum = boost::dynamic_bitset<>(massBinNumber);

    if (ms_level == 1)
    {
      unionPrevMassBins();
    }
    auto perMassChargeRanges = updateMassBins(//massIntensities,
        mzBinIntensities);

    getCandidatePeakGroups(//massIntensities,
        perMassChargeRanges);

    scoreAndFilterPeakGroups();

    if(ms_level == 1)
    {
      removeOverlappingPeakGroups(tolerance[ms_level - 1]);
    }
    else
    {
      removeOverlappingPeakGroupsWithNominalMass();
    }

    removeHarmonicPeakGroups(tolerance[ms_level - 1]); //

    if (ms_level == 1)
    {
      while (!prev_rt_vector.empty() &&
             deconvoluted_spectrum->getOriginalSpectrum().getRT() - prev_rt_vector[0] > rt_window)//
      {
        prev_rt_vector.erase(prev_rt_vector.begin());
        prev_mass_bin_vector.erase(prev_mass_bin_vector.begin());
        prev_minbin_logmass_vector.erase(prev_minbin_logmass_vector.begin());
      }
      std::vector<Size> mb;
      mb.reserve(deconvoluted_spectrum->size());
      for (auto& pg : *deconvoluted_spectrum)//filteredPeakGroups
      {
        pg.shrink_to_fit();
        //if (massBinsForThisSpectrum[pg.massBinIndex])
        //{
        double massDelta = avg.getAverageMassDelta(pg.getMonoMass());
        Size pgBin = getBinNumber(pg.getMonoMass() + massDelta, mass_bin_min_value, bin_width[ms_level - 1]);
        mb.push_back(pgBin);
        //  }
      }

      prev_rt_vector.push_back(deconvoluted_spectrum->getOriginalSpectrum().getRT());
      prev_mass_bin_vector.push_back(mb); //
      prev_minbin_logmass_vector.push_back(mass_bin_min_value);
      prev_rt_vector.shrink_to_fit();
      prev_mass_bin_vector.shrink_to_fit();
      prev_minbin_logmass_vector.shrink_to_fit();
    }
  }


  double
  FLASHDeconvAlgorithm::getCosine(const std::vector<double>& a,
                                  const int& aStart,
                                  const int& aEnd,
                                  const IsotopeDistribution& b,
                                  const int& bSize,
                                  const double& bNorm,
                                  const int offset)
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

  /*
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
  }*/


  double FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(const double mass,
                                                                        const std::vector<double>& per_isotope_intensities,
                                                                        int& offset,
                                                                        const PrecalculatedAveragine& avg)
  {
    auto iso = avg.get(mass);
    double isoNorm = avg.getNorm(mass);

    int isoSize = (int) iso.size();

    offset = 0;
    double maxCosine = -1;
    int isotopeLength = 0;
    int maxIsotopeIndex = 0, minIsotopeIndex = -1;

    for (int i = 0; i < avg.getMaxIsotopeIndex(); i++)
    {
      if (per_isotope_intensities[i] <= 0)
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


    int maxCntr = 0;
    for (int f = -isoSize - minIsotopeIndex; f <= maxIsotopeIndex; f++)
    {
      double cos = getCosine(per_isotope_intensities,
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


  bool FLASHDeconvAlgorithm::checkChargeDistribution(const std::vector<double>& per_charge_intensity)
  {
    double maxPerChargeIntensity = .0;
    int nonZeroStart = -1, nonZeroEnd = 0;
    int chargeRange = current_max_charge - min_charge + 1;
    for (int i = 0; i < chargeRange; i++)
    {
      if (per_charge_intensity[i] > 0)
      {
        maxPerChargeIntensity = std::max(maxPerChargeIntensity, per_charge_intensity[i]);
        if (nonZeroStart < 0)
        {
          nonZeroStart = i;
        }
        nonZeroEnd = i;
      }
    }

    int prevCharge = nonZeroStart;

    int n_r = 0;
    double factor = 5.0;
    double intThreshold = maxPerChargeIntensity / factor;//intensities[intensities.size()*95/100] / 5.0;
    for (int k = prevCharge + 1; k <= nonZeroEnd; k++)
    {
      if (per_charge_intensity[k] <= intThreshold)
      {
        continue;
      }

      if (k - prevCharge == 1)
      {
        n_r++;
      }
      if (n_r >= min_support_peak_count[ms_level - 1])//
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
    filteredPeakGroups.reserve(deconvoluted_spectrum->size());
    double minThreshold = std::numeric_limits<double>::max();
    int chargeRange = current_max_charge - min_charge + 1;

    Size maxc = max_mass_count.size() > ms_level - 1 ? max_mass_count[ms_level - 1] : -1;
    Size minc = min_mass_count.size() > ms_level - 1 ? min_mass_count[ms_level - 1] : -1;

    if (maxc > 0 || minc > 0)
    {
      std::vector<double> intensities;
      intensities.reserve(deconvoluted_spectrum->size());

      for (auto& pg : *deconvoluted_spectrum)
      {
        if (pg.getMonoMass() < min_mass || pg.getMonoMass() > max_mass)
        {
          continue;
        }
        intensities.push_back(pg.getIntensity());
      }
      sort(intensities.begin(), intensities.end());

      if (intensities.size() > (Size) minc)
      {
        minThreshold = intensities[intensities.size() - minc];
      }
    }


    for (auto& pg : *deconvoluted_spectrum)
    {
      bool pass = false;
      if (pg.getIntensity() >= minThreshold)
      {
        pass = true; //
      }

      auto perIsotopeIntensity = std::vector<double>(avg.getMaxIsotopeIndex(), 0);
      auto perChargeIntensity = std::vector<double>(chargeRange, 0);

      auto indices = calculatePerChargeIsotopeIntensity(
          perIsotopeIntensity, perChargeIntensity,
          avg.getMaxIsotopeIndex(), pg);

      double cs = getChargeFitScore(perChargeIntensity, chargeRange);
      pg.setChargeScore(cs);

      if (ms_level == 1)
      {
        bool isChargeWellDistributed = checkChargeDistribution(perChargeIntensity);

        if (!isChargeWellDistributed)
        {
          if (!pass)
          {
            continue;
          }
        }
      }

      int offset = 0;
      double cos = getIsotopeCosineAndDetermineIsotopeIndex(pg[0].getUnchargedMass(),
                                                          perIsotopeIntensity,
                                                          offset, avg);
      pg.setIsotopeCosine(cos);

      if (pg.empty() ||
          (pg.getIsotopeCosine() <=
           min_isotope_cosine[ms_level - 1]))// (msLevel <= 1 ? param.minIsotopeCosineSpec : param.minIsotopeCosineSpec2)))
      {
        if (!pass)
        {
          continue;
        }
      }

      pg.updateMassesAndIntensity(offset, avg.getMaxIsotopeIndex());
      if (pg.getMonoMass() < min_mass || pg.getMonoMass() > max_mass)
      {
        continue;
      }
      auto iso = avg.get(pg.getMonoMass());
      double isoNorm = avg.getNorm(pg.getMonoMass());
      int isoSize = (int) iso.size();
      float totalNoise = .0;
      float totalSignal = .0;
      //auto perChargeMaxIntensity = std::vector<double>(chargeRange);

      auto crange = pg.getChargeRange();
      for (int charge = std::get<0>(crange); charge <= std::get<1>(crange); charge++)
      {
        int j = charge - min_charge;
        if (perChargeIntensity[j] <= 0)
        {
          continue;
        }
        auto perIsotopeIntensities = std::vector<double>(avg.getMaxIsotopeIndex(), 0);

        int minIsotopeIndex = avg.getMaxIsotopeIndex();
        int maxIsotopeIndex = 0;

        double maxIntensity = .0;
        //double sumIntensity = .0;
        double sp = .0;

        for (auto& p:pg)
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
          //sumIntensity += p.intensity;
          minIsotopeIndex = minIsotopeIndex < p.isotopeIndex ? minIsotopeIndex : p.isotopeIndex;
          maxIsotopeIndex = maxIsotopeIndex < p.isotopeIndex ? p.isotopeIndex : maxIsotopeIndex;

          //min_mz = min_mz < p.mz ? min_mz : p.mz;
          // max_mz = max_mz > p.mz ? max_mz : p.mz;
          if (maxIntensity < p.intensity)
          {
            maxIntensity = p.intensity;
            // perChargeMaxIntensity[j] = maxIntensity;
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

        double _cos = getCosine(perIsotopeIntensities,
                              minIsotopeIndex,
                              maxIsotopeIndex,
                              iso,
                              isoSize,
                              isoNorm,
                              0);

        double cos2 = _cos * _cos;

        pg.setChargeIsotopeCosine(charge, _cos);
        pg.setChargeIntensity(charge, perChargeIntensity[j]);

        double dno = (1 - cos2) * sp + pg.getChargeSNR(charge) + 1;
        double no = cos2 * sp + 1;

        pg.setChargeSNR(charge, no / dno);

        totalNoise += dno;
        totalSignal += no;
      }

      pg.setSNR(totalSignal / totalNoise);
      pg.setQScore(-10000);

      for (int charge = std::get<0>(crange); charge <= std::get<1>(crange); charge++)
      {
        if (pg.getChargeIntensity(charge) <= 0)
        {
          continue;
        }
        int j = charge - min_charge;

        double score = QScore::getQScore(&pg, charge);

        if (score <= pg.getQScore())
        {
          continue;
        }
        pg.setRepCharge(charge);
        pg.setQScore(score);
      }

      double maxQScoreMzStart = pg.getMonoMass() * 2;
      double maxQScoreMzEnd = .0;
      for (auto& p:pg)
      {
        if (p.charge != pg.getRepCharge())
        {
          continue;
        }
        if (p.isotopeIndex > isoSize)
        {
          continue;
        }

        maxQScoreMzStart = maxQScoreMzStart < p.mz ? maxQScoreMzStart : p.mz;
        maxQScoreMzEnd = maxQScoreMzEnd > p.mz ? maxQScoreMzEnd : p.mz;
      }
      if (maxQScoreMzStart > maxQScoreMzEnd)
      {
        continue;
      }
      pg.setMaxQScoreMzRange(maxQScoreMzStart, maxQScoreMzEnd);
      filteredPeakGroups.push_back(pg);
    }
    deconvoluted_spectrum->swap(filteredPeakGroups);

    if (ms_level > 1)
    {
      filterPeakGroupsByIsotopeCosine(maxc);
    }
    else
    {
      filterPeakGroupsByQScore(maxc);
    }
  }

  void FLASHDeconvAlgorithm::filterPeakGroupsByIsotopeCosine(const int current_max_mass_count)
  {
    if (current_max_mass_count <= 0 || deconvoluted_spectrum->size() <= (Size) current_max_mass_count)
    {
      return;
    }

    std::vector<double> scores;
    scores.reserve(deconvoluted_spectrum->size());
    for (auto& pg : *deconvoluted_spectrum)
    {
      scores.push_back(pg.getIsotopeCosine());
    }

    sort(scores.begin(), scores.end());

    auto newPeakGroups = std::vector<PeakGroup>();
    newPeakGroups.reserve(deconvoluted_spectrum->size());
    double threshold = scores[scores.size() - current_max_mass_count];
    for (auto& pg : *deconvoluted_spectrum)
    {
      if (newPeakGroups.size() > current_max_mass_count)
      {
        break;
      }

      if (pg.getIsotopeCosine() >= threshold)
      {
        newPeakGroups.push_back(pg);
      }
    }
    deconvoluted_spectrum->swap(newPeakGroups);
  }

  void FLASHDeconvAlgorithm::filterPeakGroupsByQScore(const int current_max_mass_count)
  {
    if (current_max_mass_count <= 0 || deconvoluted_spectrum->size() <= (Size) current_max_mass_count)
    {
      return;
    }

    Size mc = (Size) current_max_mass_count;
    std::vector<double> scores;
    scores.reserve(deconvoluted_spectrum->size());
    for (auto& pg : *deconvoluted_spectrum)
    {
      scores.push_back(pg.getQScore());
    }

    sort(scores.begin(), scores.end());

    auto newPeakGroups = std::vector<PeakGroup>();
    newPeakGroups.reserve(deconvoluted_spectrum->size());
    double threshold = scores[scores.size() - mc];
    for (auto& pg : *deconvoluted_spectrum)
    {
      if (newPeakGroups.size() > mc)
      {
        break;
      }

      if (pg.getQScore() >= threshold)
      {
        newPeakGroups.push_back(pg);
      }
    }
    deconvoluted_spectrum->swap(newPeakGroups);
  }


  void FLASHDeconvAlgorithm::removeHarmonicPeakGroups(const double tol)
  {
    sort(deconvoluted_spectrum->begin(), deconvoluted_spectrum->end());
    std::vector<PeakGroup> merged;
    merged.reserve(deconvoluted_spectrum->size());
    std::vector<double> masses;
    masses.reserve(deconvoluted_spectrum->size());
    for (auto& peakGroup : *deconvoluted_spectrum)
    {
      masses.push_back(peakGroup.getMonoMass());
    }
    for (auto& pg : *deconvoluted_spectrum)
    {
      bool select = true;
      for (int h = 2; h <= 3; h++)
      {
        for (int k = 0; k < 2; k++)
        {
          for (int i = -2; i <= 2; ++i)
          {
            double omass = pg.getMonoMass() + i * Constants::ISOTOPE_MASSDIFF_55K_U;
            double hmass = k == 0 ? omass * h : omass / h;
            double massTol = 2 * hmass * tol;
            auto iter = std::lower_bound(masses.begin(), masses.end(), hmass - massTol);
            Size j = iter - masses.begin();

            if (j >= 0 && j < deconvoluted_spectrum->size())
            {
              for (; j < deconvoluted_spectrum->size(); j++)
              {
                auto& pgo = (*deconvoluted_spectrum)[j];
                if (hmass - pgo.getMonoMass() > massTol)
                {
                  continue;
                }

                if (!select || pgo.getMonoMass() - hmass > massTol)
                {
                  break;
                }
                select &= pg.getIntensity() >= pgo.getIntensity();
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
    deconvoluted_spectrum->swap(merged);
  }

  void FLASHDeconvAlgorithm::removeOverlappingPeakGroups(const double tol)
  {
    int isoLength = 1; // inclusive
    std::vector<PeakGroup> filtered;
    filtered.reserve(deconvoluted_spectrum->size());
    sort(deconvoluted_spectrum->begin(), deconvoluted_spectrum->end());

    for (Size i = 0; i < deconvoluted_spectrum->size(); i++)
    {
      bool select = true;
      auto& pg = (*deconvoluted_spectrum)[i];

      if (pg.getMonoMass() <= 0)
      {
        continue;
      }
      double massTol = pg.getMonoMass() * tol * 2;

      int j = i + 1;
      for (int l = 0; l <= isoLength; l++)
      {
        double off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j < deconvoluted_spectrum->size(); j++)
        {
          auto& pgo = (*deconvoluted_spectrum)[j];
          if (l != 0 && pgo.getMonoMass() - pg.getMonoMass() < off - massTol)
          {
            continue;
          }

          if (!select || pgo.getMonoMass() - pg.getMonoMass() > off + massTol)
          {
            break;
          }
          select &= pg.getIsotopeCosine() > pgo.getIsotopeCosine();
        }
      }

      if (!select)
      {
        continue;
      }

      j = i - 1;
      for (int l = 0; l <= isoLength; l++)
      {
        double off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j >= 0; j--)
        {
          auto& pgo = (*deconvoluted_spectrum)[j];

          if (l != 0 && pg.getMonoMass() - pgo.getMonoMass() < off - massTol)
          {
            continue;
          }

          if (!select || pg.getMonoMass() - pgo.getMonoMass() > off + massTol)
          {
            break;
          }
          select &= pg.getIsotopeCosine() > pgo.getIsotopeCosine();
        }
      }
      if (!select)
      {
        continue;
      }
      filtered.push_back(pg);
    }
    deconvoluted_spectrum->swap(filtered);
  }

  void FLASHDeconvAlgorithm::removeOverlappingPeakGroupsWithNominalMass()
  {
    int isoLength = 1; // inclusive
    std::vector<PeakGroup> filtered;
    filtered.reserve(deconvoluted_spectrum->size());
    sort(deconvoluted_spectrum->begin(), deconvoluted_spectrum->end());

    for (Size i = 0; i < deconvoluted_spectrum->size(); i++)
    {
      bool select = true;
      auto& pg = (*deconvoluted_spectrum)[i];

      if (pg.getMonoMass() <= 0)
      {
        continue;
      }
      int j = i + 1;
      for (int l = 0; l <= isoLength; l++)
      {
        double off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j < deconvoluted_spectrum->size(); j++)
        {
          auto& pgo = (*deconvoluted_spectrum)[j];
          if (l != 0 && getNominalMass(pgo.getMonoMass()) < getNominalMass(pg.getMonoMass()  + off))
          {
            continue;
          }

          if (!select || getNominalMass(pgo.getMonoMass()) > getNominalMass(pg.getMonoMass()  + off))
          {
            break;
          }
          select &= pg.getIsotopeCosine() > pgo.getIsotopeCosine();
        }
      }

      if (!select)
      {
        continue;
      }

      j = i - 1;
      for (int l = 0; l <= isoLength; l++)
      {
        double off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j >= 0; j--)
        {
          auto& pgo = (*deconvoluted_spectrum)[j];

          if (l != 0 && getNominalMass(pg.getMonoMass()) < getNominalMass(pgo.getMonoMass()  + off))
          {
            continue;
          }

          if (!select || getNominalMass(pg.getMonoMass()) > getNominalMass(pgo.getMonoMass()  + off))
          {
            break;
          }
          select &= pg.getIsotopeCosine() > pgo.getIsotopeCosine();
        }
      }
      if (!select)
      {
        continue;
      }
      filtered.push_back(pg);
    }
    deconvoluted_spectrum->swap(filtered);
  }

  void FLASHDeconvAlgorithm::reassignPeaksinPeakGroups()
  {
    /*
    //auto maxSNR = std::vector<double>(log_mz_peaks.size() , 0);
    auto maxIntensity = std::vector<double>(log_mz_peaks.size(), 0);
    //auto maxIntensityCharge = std::vector<int>(log_mz_peaks.size() , 0);

    for (auto& pg : *deconvoluted_spectrum)
    {

      auto intensity = pg.intensity / (pg.maxCharge - pg.minCharge + 1);

      for (auto& p: pg)
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

    for (auto& pg : *deconvoluted_spectrum)
    {
      //auto snr = pg.total_snr;
      auto intensity = pg.intensity / (pg.maxCharge - pg.minCharge + 1);

      std::vector<LogMzPeak> tmp;
      tmp.swap(pg);
      pg.reserve(tmp.size());
      for (auto& p: tmp)
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
    }*/
  }


  std::vector<int> FLASHDeconvAlgorithm::calculatePerChargeIsotopeIntensity(
      std::vector<double>& per_isotope_intensity,
      std::vector<double>& per_charge_intensity,
      const int max_isotope_count,
      PeakGroup& pg)
  {
    int minPgCharge = INT_MAX;
    int maxPgCharge = INT_MIN;
    //    double maxIntensity = -1;
    int maxIntChargeIndex = -1;
    //    double maxIntensity2 = -1;
    int maxIntIsoIndex = -1;

    for (auto& p : pg)
    {
      if (p.isotopeIndex < 0 || p.isotopeIndex >= max_isotope_count)
      {
        continue;
      }
      minPgCharge = std::min(minPgCharge, p.charge);
      maxPgCharge = std::max(maxPgCharge, p.charge);

      int index = p.charge - min_charge;
      per_isotope_intensity[p.isotopeIndex] += p.intensity;
      per_charge_intensity[index] += p.intensity;
    }
    pg.setChargeRange(minPgCharge, maxPgCharge);

    return std::vector<int>{maxIntChargeIndex, maxIntIsoIndex};
  }

  double FLASHDeconvAlgorithm::getCosine(const std::vector<double>& a, const std::vector<double>& b, const int off)
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
    if (d <= 0 || n <= 0)
    {
      return 0;
    }

    return n / sqrt(d);
  }


  /*double FLASHDeconvAlgorithm::getChargeFitScore(double *per_charge_intensity, int charge_range)
  {
    double maxPerChargeIntensity = .0;
    std::vector<double> xs;
    std::vector<double> ys;

    xs.reserve(+2);
    ys.reserve(charge_range + 2);

    for (int i = 0; i < charge_range; i++)
    {
      maxPerChargeIntensity = std::max(maxPerChargeIntensity, per_charge_intensity[i]);
    }

    double th = maxPerChargeIntensity * .02;// as recommended in the original paper...
    int first = -1, last = 0;
    for (int i = 0; i < charge_range; i++)
    {
      if (per_charge_intensity[i] <= th)
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
      ys.push_back((1 + per_charge_intensity[i]));
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
      auto& x = xs[i];
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

  double FLASHDeconvAlgorithm::getChargeFitScore(const std::vector<double>& per_charge_intensity, const int charge_range)
  {
    double maxPerChargeIntensity = .0;
    double sumIntensity = .0;
    int maxIndex = -1;
    int firstIndex = -1;
    int lastIndex = charge_range - 1;

    for (int i = 0; i < charge_range; i++)
    {
      sumIntensity += per_charge_intensity[i];
      if (per_charge_intensity[i] <= 0){
        if(firstIndex<0){
          firstIndex = i;
        }
        lastIndex = i;
      }

      if (maxPerChargeIntensity > per_charge_intensity[i])
      {
        continue;
      }
      maxPerChargeIntensity = per_charge_intensity[i];
      maxIndex = i;
    }
    firstIndex = firstIndex < 0? 0: firstIndex;

    double p = .0;
    for (int i = maxIndex; i < lastIndex - 1; i++)
    {
      double diff = per_charge_intensity[i + 1] - per_charge_intensity[i];
      double ratio = per_charge_intensity[i] / (.1 + per_charge_intensity[i + 1]);
      if (diff <= 0 && ratio < 5.0)
      {
        continue;
      }
      p += abs(diff);
    }

    for (int i = maxIndex; i > firstIndex; i--)
    {
      double diff = per_charge_intensity[i - 1] - per_charge_intensity[i];
      double ratio = per_charge_intensity[i] / (.1 + per_charge_intensity[i - 1]);

      if (diff <= 0 && ratio < 5.0)
      {
        continue;
      }
      p += abs(diff);
    }
    return std::max(.0, 1.0 - p / sumIntensity);
  }
}