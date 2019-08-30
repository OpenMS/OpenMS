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

namespace OpenMS
{

  // constructor
  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm()
  {
  }

  FLASHDeconvAlgorithm::~FLASHDeconvAlgorithm()
  {
  }

  FLASHDeconvAlgorithm::FLASHDeconvAlgorithm(const FLASHDeconvAlgorithm&)
  {
  }

  FLASHDeconvAlgorithm& FLASHDeconvAlgorithm::operator=(const FLASHDeconvAlgorithm& fd)
  {
    //ALWAYS CHECK FOR SELF ASSIGNEMT!
    if (this == &fd) return *this;
    //...
    return *this;
  }

  int FLASHDeconvAlgorithm::getNominalMass(double &m)
  {
    return (int) (m * 0.999497 + .5);
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



  std::vector<FLASHDeconvAlgorithm::PeakGroup> FLASHDeconvAlgorithm::Deconvolution(MSExperiment &map, Parameter &param,
                                         PrecalcularedAveragine &averagines, int &specCntr, int &qspecCntr,
                                         int &massCntr)
  {
    double *filter = new double[param.chargeRange];
    //long **hBinOffsets = new long *[param.chargeRange];
    auto harmonicFilter = new double *[param.hCharges.size()];
    int *random = new int[param.chargeRange];
    std::srand(std::time(nullptr));

    for (int i = 0; i < param.chargeRange; i++)
    {
      //auto r = ((float) rand()) / (float) RAND_MAX - .5;
      filter[i] = log(
          1.0 / (i + param.minCharge)); //.124 * (1 + i%2) + should be descending, and negative!  - 0.2 * (1+(i%3)/2.0))

      //cout << i<< " " << filter[i]<<endl;
    }

    for (Size k = 0; k < param.hCharges.size(); k++)
    {
      harmonicFilter[k] = new double[param.chargeRange];
      auto &hc = param.hCharges[k];
      float n = (float) (hc / 2);

      for (int i = 0; i < param.chargeRange; i++)
      {
        //auto r = ((float) rand()) / (float) RAND_MAX - .5;log(1.0 / (i + n / hc + param.minCharge));
        harmonicFilter[k][i] = log(
            1.0 / (i - n / hc +
                   param.minCharge)); //.124 * (1 + i%2) + should be descending, and negative!  - 0.2 * (1+(i%3)/2.0))
      }
    }

    if (param.jitter != 0)
    {
      double *tfilter = new double[param.chargeRange];
      auto m = filter[0];
      auto M = filter[param.chargeRange - 1];
      for (int i = 0; i < param.chargeRange; i++)
      {
        tfilter[i] = -filter[param.chargeRange - i - 1] + M + m;
      }
      filter = tfilter;
    }

    float prevProgress = .0;
    std::vector<PeakGroup> allPeakGroups;
    allPeakGroups.reserve(100000);
    //to overlap previous mass bins.
    std::vector<std::vector<Size>> prevMassBinVector;
    std::vector<double> prevMinBinLogMassVector;
    //vector <MSSpectrum> prevSpectra;

    for (auto it = map.begin(); it != map.end(); ++it)
    {
      if ((int) it->getMSLevel() > param.maxMSLevel)
      {
        continue;
      }

      //  if (it->getRT() < 3290 || it->getRT() > 3330)
      // {
      //  continue;
      // }

      float progress = (float) (it - map.begin()) / map.size();
      if (progress > prevProgress + .01)
      {
        printProgress(progress); //
        prevProgress = progress;
      }

      specCntr++;

      auto logMzPeaks = FLASHDeconvAlgorithm::getLogMzPeaks(*it, param);
      if (logMzPeaks.empty())
      {
        continue;
      }
      auto peakGroups = FLASHDeconvAlgorithm::getPeakGroupsFromSpectrum(logMzPeaks, filter, harmonicFilter,
                                                  prevMassBinVector, prevMinBinLogMassVector,
                                                  averagines,
                                                  param, specCntr);


      if (peakGroups.empty())
      {
        continue;
      }
      //vector<PeakGroup>().swap(peakGroups);
      //prevPgs = peakGroups;
      qspecCntr++;

      //allPeakGroups.reserve(allPeakGroups.size() + peakGroups.size());

      for (auto &pg : peakGroups)
      {
        massCntr++;
        pg.spec = &(*it);
        pg.massIndex = massCntr;
        pg.specIndex = qspecCntr;
        pg.massCntr = (int) peakGroups.size();
        allPeakGroups.push_back(pg);
      }
    }

    delete[] filter;

    for (Size k = 0; k < param.hCharges.size(); k++)
    {
      delete[] harmonicFilter[k];
    }
    delete[] harmonicFilter;
    delete[] random;
    printProgress(1); //
    //allPeakGroups.shrink_to_fit();
    return allPeakGroups; //
  }

  void FLASHDeconvAlgorithm::printProgress(float progress)
  {
    //return;
    int barWidth = 70;
    std::cout << "[";
    int pos = (int) (barWidth * progress);
    for (int i = 0; i < barWidth; ++i)
    {
      if (i < pos)
      {
        std::cout << "=";
      }
      else if (i == pos)
      {
        std::cout << ">";
      }
      else
      {
        std::cout << " ";
      }
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
  }

  std::vector<FLASHDeconvAlgorithm::LogMzPeak> FLASHDeconvAlgorithm::getLogMzPeaks(MSSpectrum &spec, const Parameter &param)
  {
    std::vector<LogMzPeak> logMzPeaks;
    logMzPeaks.reserve(spec.size());
    for (auto &peak: spec)
    {
      if (peak.getIntensity() <= param.intensityThreshold)
      {
        continue;
      }
      LogMzPeak logMzPeak(peak);
      logMzPeaks.push_back(logMzPeak);
    }

    //logMzPeaks.shrink_to_fit();

    return logMzPeaks;
  }


  std::vector<FLASHDeconvAlgorithm::PeakGroup> FLASHDeconvAlgorithm::getPeakGroupsFromSpectrum(std::vector<LogMzPeak> &logMzPeaks, double *filter, double **harmonicFilter, std::vector<std::vector<Size>> &prevMassBinVector, std::vector<double> &prevMinBinLogMassVector, PrecalcularedAveragine &averagines, const Parameter &param, int &specCntr)
  {
    int sn = 1;
    double massDelta = (param.maxMass - param.minMass) / sn;

    double minMass = param.minMass + massDelta * (specCntr % sn);
    double maxMass = minMass + massDelta;

    double massBinMaxValue = std::min(
        logMzPeaks[logMzPeaks.size() - 1].logMz -
        filter[param.chargeRange - param.minContinuousChargePeakCount - 1],
        log(param.maxMass));

    double massBinMinValue = logMzPeaks[0].logMz - filter[param.minContinuousChargePeakCount];
    double mzBinMinValue = logMzPeaks[0].logMz;
    double mzBinMaxValue = logMzPeaks[logMzPeaks.size() - 1].logMz;
    Size massBinNumber = getBinNumber(massBinMaxValue, massBinMinValue, param.binWidth) + 1;

    long *binOffsets = new long[param.chargeRange];

    for (int i = 0; i < param.chargeRange; i++)
    {
      binOffsets[i] = (long) round((mzBinMinValue - filter[i] - massBinMinValue) *
                                   param.binWidth);//rand() %10000 + 100;
    }

    auto hBinOffsets = new long *[param.hCharges.size()];
    for (Size k = 0; k < param.hCharges.size(); k++)
    {
      hBinOffsets[k] = new long[param.chargeRange];
      for (int i = 0; i < param.chargeRange; i++)
      {
        hBinOffsets[k][i] = (long) round((mzBinMinValue - harmonicFilter[k][i] - massBinMinValue) *
                                         param.binWidth);
      }
    }


    if (param.jitter > 0)
    {

      for (int i = 0; i < param.chargeRange - 1; i++)
      {
        int diff = binOffsets[i + 1] - binOffsets[i];

        binOffsets[i] += rand() % diff;
      }

      binOffsets[param.chargeRange - 1] += rand() % 50;
      //random[i] = rand() % 400 - 200;
      //
      //param.jitter >0 ? random[i] : 0)
    }

    Size mzBinNumber = getBinNumber(mzBinMaxValue, mzBinMinValue, param.binWidth) + 1;
    float *intensities = new float[mzBinNumber];

    auto mzBins = getMzBins(logMzPeaks, mzBinMinValue, mzBinNumber, param.binWidth, intensities);

    boost::dynamic_bitset<> massBins(massBinNumber);
    float *sumLogIntensities = new float[massBinNumber];

    auto unionMassBins = getUnionMassBin(massBins, massBinMinValue, prevMassBinVector, prevMinBinLogMassVector,
                                         param);
    auto perMassChargeRanges = getMassBins(massBins, mzBins, massBinMinValue,
                                           sumLogIntensities,
                                           binOffsets,
                                           hBinOffsets,
                                           unionMassBins,
                                           intensities, param,
                                           minMass, maxMass);

    //cout<<1<<endl;
    auto peakGroups = getPeakGroupsWithMassBins(unionMassBins,
                                                logMzPeaks,
                                                mzBinMinValue, massBinMinValue, sumLogIntensities,
                                                binOffsets,
                                                perMassChargeRanges,
                                                param);

    //cout<<2<<endl;
    // cout<<peakGroups.size() << " "; //
    auto filteredPeakGroups = scoreAndFilterPeakGroups(peakGroups, averagines, param);
    if (prevMassBinVector.size() > 0 && prevMassBinVector.size() >= (Size) param.numOverlappedScans)
    {
      prevMassBinVector.erase(prevMassBinVector.begin());
      prevMinBinLogMassVector.erase(prevMinBinLogMassVector.begin());
    }
    //cout<<3<<endl;
    //cout<<filteredPeakGroups.size() << endl; //

    std::vector<Size> mb;
    mb.reserve(filteredPeakGroups.size());
    for (auto &pg : filteredPeakGroups)//filteredPeakGroups
    {
      pg.peaks.shrink_to_fit();
      if (massBins[pg.massBinIndex])
      {
        mb.push_back(pg.massBinIndex);
      }
    }

    prevMassBinVector.push_back(mb); //
    prevMinBinLogMassVector.push_back(massBinMinValue);

    prevMassBinVector.shrink_to_fit();
    prevMinBinLogMassVector.shrink_to_fit();

    //cout<<4<<endl;

    for (Size k = 0; k < param.hCharges.size(); k++)
    {
      delete[] hBinOffsets[k];
    }
    delete[] hBinOffsets;

    delete[] binOffsets;
    for (int i = 0; i < 3; i++)
    {
      delete[] perMassChargeRanges[i]; // delete array within matrix
    }// delete actual matrix
    delete[] perMassChargeRanges;
    delete[] intensities;
    delete[] sumLogIntensities;
    return filteredPeakGroups;
  }

  boost::dynamic_bitset<> FLASHDeconvAlgorithm::getUnionMassBin(boost::dynamic_bitset<> &massBins, double &massBinMinValue,
                                                       std::vector<std::vector<Size>> &prevMassBinVector,
                                                       std::vector<double> &prevMassBinMinValue, const Parameter &param)
  {
    boost::dynamic_bitset<> u(massBins.size());
    if (u.size() == 0)
    {
      return u;
    }
    for (Size i = 0; i < prevMassBinVector.size(); i++)
    {
      auto &pmb = prevMassBinVector[i];
      if (pmb.empty())
      {
        continue;
      }
      long shift = (long) (round((massBinMinValue - prevMassBinMinValue[i]) * param.binWidth));

      for (auto &index : pmb)
      {
        long j = index - shift;
        if (j < 0)
        {
          continue;
        }
        if ((Size) j >= u.size())
        {
          break;
        }
        u[j] = true;
      }
    }
    return u;
  }

  std::vector<FLASHDeconvAlgorithm::PeakGroup> FLASHDeconvAlgorithm::getPeakGroupsWithMassBins(boost::dynamic_bitset<> &unionedMassBins,
      //boost::dynamic_bitset<> &massBins,
                                                     std::vector<LogMzPeak> &logMzPeaks,
                                                     double &mzBinMinValue,
                                                     double &massBinMinValue,
                                                     float *sumLogIntensities,
                                                     long *binOffsets,
                                                     Byte **chargeRanges,
                                                     const Parameter &param)
  {
    double binWidth = param.binWidth;
    double tol = param.tolerance * 2;
    int minCharge = param.minCharge;
    int chargeRange = param.chargeRange;
    int maxIsotopeCount = param.maxIsotopeCount;

    int logMzPeakSize = (int) logMzPeaks.size();
    Size massBinSize = unionedMassBins.size();
    int *currentPeakIndex = new int[param.chargeRange];
    std::fill_n(currentPeakIndex, param.chargeRange, 0); //

    std::vector<PeakGroup> peakGroups;
    peakGroups.reserve(unionedMassBins.count());
    auto &minChargeRanges = chargeRanges[0];
    auto &maxChargeRanges = chargeRanges[1];
    auto &mzChargeRanges = chargeRanges[2];

    //vector<Size> toRemove;
    //toRemove.reserve(1000);
    //boost::dynamic_bitset<> tmp(logMzPeaks.size());

    auto massBinIndex = unionedMassBins.find_first();
    // Size lastSetMassBinIndex = unionedMassBins.size();
    Size *peakBinNumbers = new Size[logMzPeakSize];
    for (int i = 0; i < logMzPeakSize; i++)
    {
      peakBinNumbers[i] = getBinNumber(logMzPeaks[i].logMz, mzBinMinValue, binWidth);
    }

    while (massBinIndex != unionedMassBins.npos)
    {
      double logM = getBinValue(massBinIndex, massBinMinValue, binWidth);
      double diff = Constants::C13C12_MASSDIFF_U / exp(logM);
      double isoLogM1 = logM - diff;
      double isoLogM2 = logM + diff;

      auto b1 = getBinNumber(isoLogM1, massBinMinValue, binWidth);
      if (b1 > 0)
      {
        if (sumLogIntensities[massBinIndex] < sumLogIntensities[b1])
        {
          massBinIndex = unionedMassBins.find_next(massBinIndex);
          continue;
        }
      }

      auto b2 = getBinNumber(isoLogM2, massBinMinValue, binWidth);
      if (b2 < unionedMassBins.size())
      {
        if (sumLogIntensities[massBinIndex] < sumLogIntensities[b2])
        {
          massBinIndex = unionedMassBins.find_next(massBinIndex);
          continue;
        }
      }

      if (sumLogIntensities[b1] == 0 && sumLogIntensities[b2] == 0)
      {
        massBinIndex = unionedMassBins.find_next(massBinIndex);
        continue;
      }

      int isoOff = 0;
      PeakGroup pg;
      pg.reserve(chargeRange * 30);

      for (auto j = minChargeRanges[massBinIndex]; j <= maxChargeRanges[massBinIndex]; j++)
      {
        long &binOffset = binOffsets[j];
        auto bi = massBinIndex - binOffset;
        if (mzChargeRanges[bi] < chargeRange && mzChargeRanges[bi] != j)
        {
          continue;
        }

        int charge = j + minCharge;
        auto &cpi = currentPeakIndex[j];
        double maxIntensity = 0.0;
        int maxIntensityPeakIndex = -1;

        while (cpi < logMzPeakSize - 1)
        {
          //auto bi = peakBinNumbers[cpi] + binOffset;
          if (peakBinNumbers[cpi] == bi)
          {
            auto intensity = logMzPeaks[cpi].orgPeak->getIntensity();
            if (intensity > maxIntensity)
            {
              maxIntensity = intensity;
              maxIntensityPeakIndex = cpi;
            }
          }
          else if (peakBinNumbers[cpi] > bi)
          {
            break;
          }
          cpi++;
        }

        if (maxIntensityPeakIndex >= 0)
        {
          const double mz = logMzPeaks[maxIntensityPeakIndex].orgPeak->getMZ();
          //double &logMz = logMzPeaks[maxIntensityPeakIndex].logMz;
          const double isof = Constants::C13C12_MASSDIFF_U / charge;
          double mzDelta = tol * mz;
          //cout<<mzDelta<<endl;
          int maxI = 0;
          for (int d = -1; d <= 1; d += 2)
          { // negative then positive direction.
            int peakIndex = maxIntensityPeakIndex + (d < 0 ? d : 0);
            int lastPeakIndex = -100;
            for (int i = 0; i < maxIsotopeCount && (peakIndex >= 0 && peakIndex < logMzPeakSize); i++)
            {
              maxI = std::max(maxI, i);
              const double centerMz = mz + isof * i * d;
              const double centerMzMin = centerMz - mzDelta;
              const double centerMzMax = centerMz + mzDelta;
              bool isotopePeakPresent = false;
              if (lastPeakIndex >= 0)
              {
                peakIndex = lastPeakIndex;
              }//maxIntensityPeakIndex + (d < 0 ? d : 0);
              for (; peakIndex >= 0 && peakIndex < logMzPeakSize; peakIndex += d)
              {
                const double observedMz = logMzPeaks[peakIndex].orgPeak->getMZ();
                if (observedMz < centerMzMin)
                {
                  if (d < 0)
                  {
                    break;
                  }
                  else
                  {
                    continue;
                  }
                }
                if (observedMz > centerMzMax)
                {
                  if (d < 0)
                  {
                    continue;
                  }
                  else
                  {
                    break;
                  }
                }

                isotopePeakPresent = true;
                if (peakIndex != lastPeakIndex)
                {
                  // auto &mzBin = ;
                  const auto bin = peakBinNumbers[peakIndex] + binOffset;

                  // if (mzChargeRanges[peakBinNumbers[peakIndex]] < chargeRange && mzChargeRanges[peakBinNumbers[peakIndex]] != j)
                  // {
                  //  continue;
                  //}

                  if (bin < massBinSize) //
                  {
                    LogMzPeak p(*logMzPeaks[peakIndex].orgPeak, charge, i * d);
                    pg.peaks.push_back(p);
                    //isoOff = min(isoOff, p.isotopeIndex);
                    lastPeakIndex = peakIndex;
                    //}
                    // if (massBinIndex != bin && unionedMassBins[bin])
                    //{
                    //unionedMassBins[bin] = false;
                    //if (bin>0) unionedMassBins[bin-1] = false;
                    //if (bin < massBinSize -1)unionedMassBins[bin+1] = false;
                    // toRemove.push_back(bin);
                    //   }
                  }
                }
              }
              if (!isotopePeakPresent)
              {
                break;
              }
            }
          }


          //int minIsoIndex = maxI;
          for (auto &p : pg.peaks)
          {// assign the nearest isotope index..
            if (p.charge != charge)
            {
              continue;
            }
            for (int d = -1; d <= 1; d += 2)
            { // negative then positive direction.
              int maxId = 0;
              double minMzDelta = maxI;

              for (int i = 0; i <= maxI; i++)
              {
                double centerMz = mz + isof * (p.isotopeIndex + i * d);
                double delta = abs(centerMz - p.orgPeak->getMZ());

                if (delta > minMzDelta)
                {
                  break;
                }
                maxId = i * d;
                minMzDelta = delta;
              }
              //if (maxId != 0 )cout<< maxId <<endl;
              p.isotopeIndex += maxId;
            }
            isoOff = std::min(isoOff, p.isotopeIndex);
          }

        }
      }
      //cout<<pg.peaks.size()<<endl;
      //pg.peaks.shrink_to_fit();
      if (!pg.peaks.empty())
      {
        int minIi = 10000;
        int maxIi = -10000;
        for (auto &p : pg.peaks)
        {
          minIi = std::min(minIi, p.isotopeIndex);
          maxIi = std::max(maxIi, p.isotopeIndex);
          p.isotopeIndex -= isoOff;
        }
        if (minIi != maxIi)
        {

          //pg.updateMassesAndIntensity();
          pg.massBinIndex = massBinIndex;
          /*if(lastSetMassBinIndex == massBinIndex-1){// remove duplicate
              auto prevIntensity = peakGroups[peakGroups.size()-1].intensity;
              auto currentIntensity = pg.intensity;
              if(prevIntensity < currentIntensity){
                  peakGroups[peakGroups.size()-1] = pg;
              }
          }else*/
          peakGroups.push_back(pg);
        }
        //lastSetMassBinIndex = massBinIndex;
      }

      //if (massBinIndex < unionedMassBins.size() - 1 && !unionedMassBins[massBinIndex + 1])
      //{
      // for (auto &b : toRemove)
      // {
      //unionedMassBins[b] = false;
      //}
      //toRemove.clear();
      //}

      massBinIndex = unionedMassBins.find_next(massBinIndex);
    }
    delete[] currentPeakIndex;
    delete[] peakBinNumbers;
    return peakGroups;
  }

  boost::dynamic_bitset<> FLASHDeconvAlgorithm::getMzBins(std::vector<LogMzPeak> &logMzPeaks, double &mzBinMinValue, Size &binNumber, double binWidth,
            float *intensities
  )
  {
    boost::dynamic_bitset<> mzBins(binNumber);
    // double *intensities = new double[binNumber];
    std::fill_n(intensities, binNumber, .0);

    for (auto &p : logMzPeaks)
    {
      Size bi = getBinNumber(p.logMz, mzBinMinValue, binWidth);
      if (bi >= binNumber)
      {
        continue;
      }
      mzBins.set(bi);
      intensities[bi] += p.orgPeak->getIntensity();

      auto delta = (p.logMz - getBinValue(bi, mzBinMinValue, binWidth));

      if (delta > 0)
      {
        //add bin + 1
        if (bi < binNumber - 1)
        {
          mzBins.set(bi + 1);
          intensities[bi + 1] += p.orgPeak->getIntensity();
        }
      }
      else if (delta < 0)
      {
        if (bi > 0)
        {
          mzBins.set(bi - 1);
          intensities[bi - 1] += p.orgPeak->getIntensity();
        }
      }
    }


    //  delete[] intensities;
    return mzBins;
  }

  Byte** FLASHDeconvAlgorithm::getMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                            double &massBinMinValue,
                            float *sumLogIntensities,
                            long *binOffsets,
                            long **hBinOffsets,
                            boost::dynamic_bitset<> &unionMassBins,
                            float *intensities,
                            const Parameter &param, double &minMass, double &maxMass)
  {
    long binThresholdMinMass = (long) getBinNumber(log(minMass), massBinMinValue, param.binWidth);
    long binThresholdMaxMass = (long) std::min(massBins.size(),
                                          1 + getBinNumber(log(maxMass), massBinMinValue, param.binWidth));
    boost::dynamic_bitset<> isQualified(massBins.size());

    //cout<<minMass << " " << maxMass<<endl;
    //fill_n(sumLogIntensities, massBins.size(), 0);
    //cout<<0<<endl;
    getInitialMassBins(massBins, mzBins, isQualified,
                       sumLogIntensities,
        //continuousChargePeakPairCount,
                       hBinOffsets, binOffsets,
                       intensities,
                       param);
    //cout<<1<<endl;
    //printMasses(isQualified, massBinMinValue, continuousChargePeakPairCount, param);
    //auto toSkip = (isQualified | isQualified).flip();
    auto perMassChargeRanges = getFinalMassBins(massBins, mzBins, isQualified, unionMassBins,
                                                sumLogIntensities,
        //     massBinMinValue,
        // continuousChargePeakPairCount,
                                                binOffsets,
                                                param,
                                                binThresholdMinMass,
                                                binThresholdMaxMass);
    //cout<<2<<endl;
    //printMasses(massBins, massBinMinValue, continuousChargePeakPairCount, param);

    //cout<<isQualified.count()<<" " << massBins.count()<<" ";//
    //delete[] continuousChargePeakPairCount;
    //delete[] sumLogIntensities;
    return perMassChargeRanges;
  }

//  void FLASHDeconvAlgorithm::printMasses(boost::dynamic_bitset<> &massBins, double &massBinMinValue, Byte *continuousChargePeakPairCount,
//              const Parameter &param)
//  {
//    auto index = massBins.find_first();
//    while (index != massBins.npos)
//    {
//      auto m = exp(getBinValue(index, massBinMinValue, param.binWidth));
//      if (m < 50500 && m > 50400)
//      {
//
//        cout << m <<
//             " " << (int) continuousChargePeakPairCount[index] <<
//             //" " << (int) noneContinuousChargePeakPairCount[index] <<
//             endl;
//      }
//      index = massBins.find_next(index);
//    }
//    //cout << endl;
//  }


  void FLASHDeconvAlgorithm::getInitialMassBins(boost::dynamic_bitset<> &massBins,
                                 boost::dynamic_bitset<> &mzBins,
                                 boost::dynamic_bitset<> &isQualified,
                                 float *signal,
                                 long **hBinOffsets,
                                 long *binOffsets,
                                 float *intensities,
                                 const Parameter &param
  )
  {
    int chargeRange = param.chargeRange;
    int hChargeSize = (int) param.hCharges.size();
    int minContinuousChargePeakCount = param.minContinuousChargePeakCount;
    //    long mzBinSize = (long) mzBins.size();
    long binEnd = (long) massBins.size();

    Byte *continuousChargePeakPairCount = new Byte[massBins.size()];
    std::fill_n(continuousChargePeakPairCount, massBins.size(), 0);

    Byte *prevCharges = new Byte[massBins.size()];
    std::fill_n(prevCharges, massBins.size(), (Byte) (chargeRange + 2));

    float *prevIntensities = new float[massBins.size()];
    std::fill_n(prevIntensities, massBins.size(), 1);

    //auto prevHcharges = new Byte*[hChargeSize];
    //for(auto k=0;k<hChargeSize;k++){
    //  prevHcharges[k] = new Byte[massBins.size()];
    //  fill_n(prevHcharges[k], massBins.size(), (Byte) (chargeRange + 2));
    // }

    //float* signal = new float[massBins.size()];
    std::fill_n(signal, massBins.size(), 0);
    auto noise = new float *[hChargeSize + 1];
    for (auto k = 0; k < hChargeSize + 1; k++)
    {
      noise[k] = new float[massBins.size()];
      std::fill_n(noise[k], massBins.size(), 0);
    }

    //double th = 1;//
    float factor = 4.0;
    auto mzBinIndex = mzBins.find_first();
    while (mzBinIndex != mzBins.npos)
    {
      auto &intensity = intensities[mzBinIndex];
      for (Byte j = 0; j < chargeRange; j++)
      {
        long massBinIndex = mzBinIndex + binOffsets[j];
        if (massBinIndex < 0)
        {
          continue;
        }
        if (massBinIndex >= binEnd)
        {
          break;
        }
        auto cd = prevCharges[massBinIndex] - j;

        auto &prevIntensity = prevIntensities[massBinIndex];
        float minInt = std::min(intensity, prevIntensity);
        float maxInt = std::max(intensity, prevIntensity);

        float id = maxInt / minInt;

        if (prevCharges[massBinIndex] < chargeRange && cd != 1 && id<factor)
        {
          noise[hChargeSize][massBinIndex] += minInt;
        }

        if (cd != 1 || id > factor)
        {
          continuousChargePeakPairCount[massBinIndex] = 0;
        }
        else
        {
          int maxHcharge = -1;
          float maxHint = .0;
          for (auto k = 0; k < hChargeSize; k++)
          {
            auto hmzBinIndex = massBinIndex - hBinOffsets[k][j];
            if (hmzBinIndex > 0 && hmzBinIndex < (long) mzBins.size() && mzBins[hmzBinIndex])
            {
              auto &hintensity = intensities[hmzBinIndex];
              if (hintensity > minInt
                  &&  hintensity < factor * maxInt
                //&&  hintensity > maxInt/factor
                  )
              {
                noise[k][massBinIndex] += hintensity;

                if (hintensity < maxHint)
                {
                  continue;
                }
                maxHint = hintensity;
                maxHcharge = k;
              }
            }
          }
          if (maxHcharge >= 0)
          {
            //  cout<<1<<endl;
            continuousChargePeakPairCount[massBinIndex] = 0;
          }
          else
          {
            signal[massBinIndex] += intensity;
            if (!isQualified[massBinIndex])
            {
              isQualified[massBinIndex] =
                  ++continuousChargePeakPairCount[massBinIndex] >= minContinuousChargePeakCount;
            }
          }
        }
        prevIntensity = intensity;
        prevCharges[massBinIndex] = j;
      }
      mzBinIndex = mzBins.find_next(mzBinIndex);
    }

    auto mindex = isQualified.find_first();
    while (mindex != isQualified.npos)
    {
      auto &s = signal[mindex];
      // auto msnr = s / (noise[0][mindex]);
      float maxNoise = .0;
      for (auto k = 0; k < hChargeSize + 1; k++)
      {
        maxNoise = std::max(maxNoise, noise[k][mindex]);
        // msnr = min(snr, msnr);
      }
      //s -= maxNoise;
      s -= maxNoise;
      mindex = isQualified.find_next(mindex);
    }

    delete[] prevIntensities;
    delete[] continuousChargePeakPairCount;
    delete[] prevCharges;
    for (auto k = 0; k < hChargeSize + 1; k++)
    {
      delete[] noise[k];
    }
    //for(auto k=0;k<hChargeSize; k++)
    //{
    //  delete[] prevHcharges[k];
    //}
    //delete[] prevHcharges;
    delete[] noise;
    //delete[] signal;

  }

  Byte** FLASHDeconvAlgorithm::getFinalMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                                 boost::dynamic_bitset<> &isQualified,
                                 boost::dynamic_bitset<> &unionMassBins,
                                 float *sumLogIntensities,
      // double &massBinMinValue,
                                 long *binOffsets,
                                 const Parameter &param,
                                 long &binStart, long &binEnd
  )
  {
    int chargeRange = param.chargeRange;
    // auto binWidth = param.binWidth;
    Byte *maxChargeRanges = new Byte[massBins.size()];
    std::fill_n(maxChargeRanges, massBins.size(), 0);

    Byte *minChargeRanges = new Byte[massBins.size()];
    std::fill_n(minChargeRanges, massBins.size(), chargeRange + 1);

    Byte *mzChargeRanges = new Byte[mzBins.size()];
    std::fill_n(mzChargeRanges, mzBins.size(), chargeRange + 1);
    //Byte *selected = new Byte[massBins.size()];
    //fill_n(selected, massBins.size(), 0);

    auto mzBinIndex = mzBins.find_first();
    long binSize = (long) massBins.size();
    //int minContinuousChargePeakCount = param.minContinuousChargePeakCount;

    auto toSkip = (isQualified | unionMassBins).flip();
    unionMassBins.reset();


    while (mzBinIndex != mzBins.npos)
    {
      long maxIndex = -1;
      float maxCount = -1e11;
      Byte charge = 0;

      //vector<Byte> setCharges;
      //setCharges.reserve(chargeRange);
      for (Byte j = 0; j < chargeRange; j++)
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

        //if (continuousChargePeakPairCount[massBinIndex] < round(log2(j+param.minCharge)))
        //  continue;

        auto &t = sumLogIntensities[massBinIndex];// + noneContinuousChargePeakPairCount[massBinIndex];//
        if (t == 0) { // no signal
          continue;
        }
        if (maxCount < t)
        {
          maxCount = t;
          maxIndex = massBinIndex;
          charge = j;
        }
      }

      // if (maxIndex >= 0)
      //{
      if (maxIndex > binStart && maxIndex < binEnd)
      {
        //if (sumLogIntensities[maxIndex] >  sumLogIntensities[maxIndex-1] &&
        //    sumLogIntensities[maxIndex] >  sumLogIntensities[maxIndex+1])
        {
          maxChargeRanges[maxIndex] = std::max(maxChargeRanges[maxIndex], charge);
          minChargeRanges[maxIndex] = std::min(minChargeRanges[maxIndex], charge);
          //  ++selected[maxIndex];
          //int t = max(continuousChargePeakPairCount[maxIndex],minContinuousChargePeakCount);
          // if (++selected[maxIndex] >= 3)
          //&& selected[maxIndex] >= continuousChargePeakPairCount[maxIndex])
          //{
          //massBins[maxIndex] = isQualified[maxIndex];

          massBins[maxIndex] = isQualified[maxIndex];

          mzChargeRanges[mzBinIndex] = charge;//minChargeRanges[maxIndex];//...

          unionMassBins[maxIndex] = true;
        }
        // }
      }
      //}
      mzBinIndex = mzBins.find_next(mzBinIndex);
    }


    // bool isChargeWellDistributed = checkChargeDistribution(perChargeIntensity,
    //                                                       param.chargeRange,
    //                                                      3);

    //cout<<" " <<massBins.count()<<endl; //
    Byte **chargeRanges = new Byte *[3];
    chargeRanges[0] = minChargeRanges;
    chargeRanges[1] = maxChargeRanges;
    chargeRanges[2] = mzChargeRanges;
    // delete[] selected;
    return chargeRanges;
  }

  std::vector<FLASHDeconvAlgorithm::PeakGroup> FLASHDeconvAlgorithm::scoreAndFilterPeakGroups(std::vector<PeakGroup> &peakGroups,
                                                    PrecalcularedAveragine &averagines,
                                                    const Parameter &param)
  {
    std::vector<FLASHDeconvAlgorithm::PeakGroup> filteredPeakGroups;
    filteredPeakGroups.reserve(peakGroups.size());
    double threshold = .0;

    Size mc = (Size) param.maxMassCount;
    if (mc > 0)
    {
      std::vector<double> intensities;
      intensities.reserve(peakGroups.size());

      for (auto &pg : peakGroups)
      {
        pg.updateMassesAndIntensity(averagines);
        intensities.push_back(pg.intensity);
      }

      if (intensities.size() > mc)
      {
        sort(intensities.begin(), intensities.end());
        threshold = intensities[intensities.size() - mc];
      }
      std::vector<double>().swap(intensities);
    }
    //filterPeakGroupsByIntensity(peakGroups, intensities, param);


    //auto perIsotopeMaxCharge = new int[param.maxIsotopeCount];
    //auto perIsotopeMinCharge = new int[param.maxIsotopeCount];
    //auto perChargeMaxIsotope = new int[param.chargeRange];
    //auto perChargeMinIsotope = new int[param.chargeRange];// polish later

    auto perIsotopeIntensity = new double[param.maxIsotopeCount];
    auto perChargeIntensity = new double[param.chargeRange];

    for (auto &pg : peakGroups)
    {
      if (pg.intensity < threshold)
      {
        continue;
      }
      updatePerChargeIsotopeIntensity(//perIsotopeMinCharge, perIsotopeMaxCharge,
          //perChargeMinIsotope, perChargeMaxIsotope,
          perIsotopeIntensity, perChargeIntensity,
          pg, param);

      pg.chargeCosineScore = getChargeFitScore(perChargeIntensity, param.chargeRange);

      if (pg.peaks.empty() || (pg.chargeCosineScore < param.minChargeCosineSpec))
      {
        continue;
      }

      /*bool isIsotopeSpanGood = checkSpanDistribution(perChargeMinIsotope, perChargeMaxIsotope,
                                                     param.chargeRange, 3);

      if (!isIsotopeSpanGood)
      {
        continue;
      }


      bool isChargeSpanGood = checkSpanDistribution(perIsotopeMinCharge, perIsotopeMaxCharge,
                                                    param.maxIsotopeCount, 3);

      if (!isChargeSpanGood)
      {
        continue;
      }
*/

      bool isChargeWellDistributed = checkChargeDistribution(perChargeIntensity,
                                                             param.chargeRange,
                                                             param.minContinuousChargePeakCount);

      if (!isChargeWellDistributed)
      {
        continue;
      }

      int offset = 0;
      pg.isotopeCosineScore = getIsotopeCosineAndDetermineIsotopeIndex(pg.peaks[0].getMass(),
                                                                       perIsotopeIntensity,
                                                                       param.maxIsotopeCount,
                                                                       averagines, offset);

      //double isotopeCosineThreshold = param.minIsotopeCosineSpec;// getIsotopeCosineThreshold(pg.peaks[0].getMass(),
      //   param.minIsotopeCosine, param.maxIsotopeCosine,
      //  0, 100000);

      if (pg.peaks.empty() || (pg.isotopeCosineScore < param.minIsotopeCosineSpec))
      {
        continue;
      }


      pg.updateMassesAndIntensity(averagines, offset, param.maxIsotopeCount);
      //cout<<"mass "<< pg.monoisotopicMass <<endl;


      /*updatePerChargeIsotopeIntensity(perIsotopeMinCharge, perIsotopeMaxCharge,
                                      perChargeMinIsotope, perChargeMaxIsotope,
                                      perIsotopeIntensity, perChargeIntensity,
                                      pg, param);
      if (pg.peaks[0].mass>60000){
        cout<<"iso=[";
        for(int i=0;i<param.maxIsotopeCount;i++){
          //if (perIsotopeIntensity[i]<=0){
          //  break;
          //}
          cout<<perIsotopeIntensity[i]<<" ";
        }
        cout<<"];";

        averagines.print(pg.peaks[0].mass);
      }
 */
      filteredPeakGroups.push_back(pg);

    }

    std::vector<PeakGroup>().swap(peakGroups);

    removeOverlappingPeakGroups(filteredPeakGroups, param.tolerance);
    //removeHarmonicPeakGroups(filteredPeakGroups, param);

    delete[] perIsotopeIntensity;
    delete[] perChargeIntensity;
    //delete[] perIsotopeMinCharge;
    //delete[] perIsotopeMaxCharge;
    //delete[] perChargeMinIsotope;
    //delete[] perChargeMaxIsotope;


    //filteredPeakGroups.shrink_to_fit();
    return filteredPeakGroups;
  }


  void FLASHDeconvAlgorithm::removeOverlappingPeakGroups(std::vector<PeakGroup> &pgs, double tol)
  { // pgs are sorted
    std::vector<PeakGroup> merged;
    merged.reserve(pgs.size());

    for (Size i = 0; i < pgs.size(); i++)
    {
      bool select = true;
      auto &pg = pgs[i];
      double massTol = pg.monoisotopicMass * tol * 2;

      for (Size j = i + 1; j < pgs.size(); j++)
      {
        auto &pgo = pgs[j];
        if (!select || pgo.monoisotopicMass - pg.monoisotopicMass > massTol)
        {
          break;
        }
        select &= pg.intensity > pgo.intensity;
      }

      if (!select)
      {
        continue;
      }
      for (int j = i - 1; j >= 0; j--)
      {
        auto &pgo = pgs[j];
        if (!select || pg.monoisotopicMass - pgo.monoisotopicMass > massTol)
        {
          break;
        }
        select &= pg.intensity > pgo.intensity;
      }
      if (!select)
      {
        continue;
      }
      merged.push_back(pg);
    }
    //pgs = merged;
    //merged.shrink_to_fit();
    std::vector<PeakGroup>().swap(pgs);
    merged.swap(pgs);
    //vector<PeakGroup>().swap(merged);
  }

  void FLASHDeconvAlgorithm::updatePerChargeIsotopeIntensity(//int *perIsotopeMinCharge, int *perIsotopeMaxCharge,
      //int *perChargeMinIsotope, int *perChargeMaxIsotope,
      double *perIsotopeIntensity,
      double *perChargeIntensity,
      PeakGroup &pg,
      const Parameter &param)
  {


    //fill_n(perIsotopeMinCharge, param.maxIsotopeCount, param.chargeRange);
    //fill_n(perIsotopeMaxCharge, param.maxIsotopeCount, -1);

    //fill_n(perChargeMinIsotope, param.chargeRange, param.maxIsotopeCount);
    //fill_n(perChargeMaxIsotope, param.chargeRange, -1);

    std::fill_n(perIsotopeIntensity, param.maxIsotopeCount, 0);
    std::fill_n(perChargeIntensity, param.chargeRange, 0);

    // bool *tmp = new bool[param.chargeRange * param.maxIsotopeCount];
    int minCharge = param.chargeRange + param.minCharge + 1;
    int maxCharge = 0;

    for (auto &p : pg.peaks)
    {
      if (p.isotopeIndex < 0 || p.isotopeIndex >= param.maxIsotopeCount)
      {
        continue;
      }
      minCharge = std::min(minCharge, p.charge);
      maxCharge = std::max(maxCharge, p.charge);

      int index = p.charge - param.minCharge;
      perIsotopeIntensity[p.isotopeIndex] += p.orgPeak->getIntensity();
      perChargeIntensity[index] += p.orgPeak->getIntensity();

      //perChargeMinIsotope[index] = min(perChargeMinIsotope[index], p.isotopeIndex);
      //perChargeMaxIsotope[index] = max(perChargeMaxIsotope[index], p.isotopeIndex);

      //perIsotopeMinCharge[p.isotopeIndex] = min(perIsotopeMinCharge[p.isotopeIndex], index);
      //perIsotopeMaxCharge[p.isotopeIndex] = max(perIsotopeMaxCharge[p.isotopeIndex], index);
      //  tmp[p.isotopeIndex + param.maxIsotopeCount * index]=true;
    }
    pg.maxCharge = maxCharge;
    pg.minCharge = minCharge;

  }


  double FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                         double *perIsotopeIntensities,
                                                         int perIsotopeIntensitiesSize,
                                                         PrecalcularedAveragine &averagines,
                                                         int &offset)
  {
    auto iso = averagines.get(mass);
    auto isoNorm = averagines.getNorm(mass);

    int isoSize = (int) iso.size();

    //double isoDiff = Constants::C13C12_MASSDIFF_U;// OpenMS::Math::mean(diffs.begin(), diffs.end());

    offset = 0;
    double maxCosine = -1;
    int maxIsotopeIndex = 0, minIsotopeIndex = -1;

    for (int i = 0; i < perIsotopeIntensitiesSize; i++)
    {
      if (perIsotopeIntensities[i] <= 0)
      {
        continue;
      }
      maxIsotopeIndex = i;
      if (minIsotopeIndex < 0)
      {
        minIsotopeIndex = i;
      }
    }

    for (int f = -isoSize + minIsotopeIndex; f <= maxIsotopeIndex; f++)
    {
      auto cos = FLASHDeconvAlgorithm::getCosine(perIsotopeIntensities, minIsotopeIndex, maxIsotopeIndex, iso, isoSize, isoNorm, f);

      if (maxCosine <= cos)
      {
        maxCosine = cos;
        offset = f;
      }
    }

    return maxCosine;
  }


  bool FLASHDeconvAlgorithm::checkSpanDistribution(int *mins, int *maxs, int range, int threshold)
  {
    int nonZeroStart = -1, nonZeroEnd = 0;
    int maxSpan = 0;

    for (int i = 0; i < range; i++)
    {
      if (maxs[i] >= 0)
      {
        if (nonZeroStart < 0)
        {
          nonZeroStart = i;
        }
        nonZeroEnd = i;
        maxSpan = std::max(maxSpan, maxs[i] - mins[i]);
      }
    }
    if (maxSpan <= 0)
    {
      return false;
    }

    int prevCharge = nonZeroStart;
    int n_r = 0;

    double spanThreshold = maxSpan / 1.5;//

    for (int k = nonZeroStart + 1; k <= nonZeroEnd; k++)
    {
      if (maxs[k] < 0)
      {
        continue;
      }


      if (k - prevCharge == 1)
      {
        int intersectSpan = std::min(maxs[prevCharge], maxs[k])
                            - std::max(mins[prevCharge], mins[k]);

        if (spanThreshold <= intersectSpan)
        { //
          n_r++;
        }
      }
      prevCharge = k;
    }

    if (n_r < threshold)
    {
      return -100.0;
    }

    for (int i = 2; i < std::min(12, range); i++)
    {
      for (int l = 0; l < i; l++)
      {
        int t = 0;
        prevCharge = nonZeroStart + l;
        for (int k = prevCharge + i; k <= nonZeroEnd; k += i)
        {
          if (maxs[k] < 0)
          {
            continue;
          }


          if (k - prevCharge == i)
          {
            int intersectSpan = std::min(maxs[prevCharge], maxs[k])
                                - std::max(mins[prevCharge], mins[k]);

            if (spanThreshold <= intersectSpan)
            {
              t++;
            }
          }
          prevCharge = k;
        }
        if (n_r <= t)
        {
          return false;
        }
      }
    }

    return true;
  }

  double FLASHDeconvAlgorithm::getChargeFitScore(double *perChargeIntensity,
                                  int range)
  {
    double maxPerChargeIntensity = .0;
    std::vector<double> xs;
    std::vector<double> ys;

    xs.reserve(range + 2);
    ys.reserve(range + 2);

    for (int i = 0; i < range; i++)
    {
      maxPerChargeIntensity = std::max(maxPerChargeIntensity, perChargeIntensity[i]);
    }

    double th = maxPerChargeIntensity * .02;//2%
    int first = -1, last = 0;
    for (int i = 0; i < range; i++)
    {
      if (perChargeIntensity[i] <= th)
      {
        continue;
      }
      if (first < 0)
      {
        first = i;
        // xs.push_back(first-1000);
        //ys.push_back(minInt/10);
      }

      last = i;
      //cout<<log(perChargeIntensity[i])<<endl;
      //  xs.push_back(i);
      // ys.push_back((1+perChargeIntensity[i]));
    }
    if (last - first < 2)
    {
      return 0;
    }
    //xs.push_back(first-2);
    //ys.push_back(1);
    //xs.push_back(first-1);
    //ys.push_back(1);

    for (int i = first; i <= last; i++)
    {
      //if (perChargeIntensity[i]<=0)
      //{
      //  continue;
      //}
      xs.push_back(i);
      ys.push_back((1 + perChargeIntensity[i]));
    }
    //xs.push_back(last+1);
    //ys.push_back(1);
    //xs.push_back(last+2);
    //ys.push_back(1);


    // xs.push_back(last + 1000);
    // ys.push_back((minInt/10));

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

    //cout<<m<<endl;
    //cout<<v<<endl;
    auto im = m.inverse();
    //cout<<im<<endl;
    // cout<<m<<endl;
    v(0) = t0;
    v(1) = t1;
    v(2) = t2;
    //cout<<v<<endl;
    auto abc = im * v;
    //cout<<abc<<endl;
    double mu = -abc(1) / abc(2) / 2;
    double omega = -1 / abc(2) / 2;
    //double a = (abc(0)-abc(1)*abc(1))/4/abc(2);

    if (omega <= 0)
    {
      //omega = -omega;
      //mu = -mu;
      return 0;
    }
    //cout<<"**"<<endl;
    //cout<<mu <<" * " << omega<<endl;
    std::vector<double> tys;

    //cout<<endl;
    for (Size i = 0; i < ys.size(); i++)
    {
      double ty = exp(-(xs[i] - mu) * (xs[i] - mu) / 2 / omega);
      tys.push_back(ty);
      //  cout<< xs[i]<<" "<<ys[i] << " " << ty <<endl;
    }
    // cout<<endl;
    // cout<<"cos " << getCosine(ys, tys) <<endl;
    return FLASHDeconvAlgorithm::getCosine(ys, tys);

  }


  bool FLASHDeconvAlgorithm::checkChargeDistribution(double *perChargeIntensity,
                                      int range,
                                      int threshold)
  {
    double maxPerChargeIntensity = .0;
    int nonZeroStart = -1, nonZeroEnd = 0;
    for (int i = 0; i < range; i++)
    {
      if (perChargeIntensity[i] > 0)
      {
        // intensities.push_back(intensity);
        maxPerChargeIntensity = std::max(maxPerChargeIntensity, perChargeIntensity[i]);
        if (nonZeroStart < 0)
        {
          nonZeroStart = i;
        }
        nonZeroEnd = i;
      }
    }

    int prevCharge = nonZeroStart;

    int n_r = .0;

    double intensityThreshold = maxPerChargeIntensity / 4.0;//intensities[intensities.size()*95/100] / 5.0;
    for (int k = prevCharge + 1; k <= nonZeroEnd; k++)
    {
      if (perChargeIntensity[k] <= intensityThreshold)
      {
        continue;
      }

      if (k - prevCharge == 1)
      {
        n_r++;
      }
      prevCharge = k;
    }

    if (n_r < threshold)
    {
      return false;
    }


    for (int i = 2; i < std::min(7, range); i++)
    {
      for (int l = 0; l < i; l++)
      {
        int t = 0;
        prevCharge = nonZeroStart + l;
        for (int k = prevCharge + i; k <= nonZeroEnd; k += i)
        {
          if (perChargeIntensity[k] <= intensityThreshold)
          { //
            continue;
          }
          if (k - prevCharge == i)
          {
            t++;
          }
          prevCharge = k;
        }
        if (n_r <= t)
        {
          return false;
        }
      }
    }

    return true;
  }

  double
  FLASHDeconvAlgorithm::getCosine(double *a, int &aStart, int &aEnd, IsotopeDistribution &b, int &bSize, double &bNorm, int offset)
  {
    double n = .0, d1 = .0;

    for (int j = aStart; j < aEnd; j++)
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

  void
  FLASHDeconvAlgorithm::filterPeakGroupsByIntensity(std::vector<PeakGroup> &peakGroups, std::vector<double> &intensities, const Parameter &param)
  {
    if (param.maxMassCount < 0 || intensities.size() <= (Size) param.maxMassCount)
    {
      return;
    }
    Size mc = (Size) param.maxMassCount;
    sort(intensities.begin(), intensities.end());
    auto threshold = intensities[intensities.size() - mc];
    for (auto pg = peakGroups.begin(); pg != peakGroups.end();)
    {
      if (peakGroups.size() <= mc)
      {
        break;
      }
      if (pg->intensity < threshold)
      {
        pg = peakGroups.erase(pg);
        continue;
      }
      ++pg;
    }
  }

}

