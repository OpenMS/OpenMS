//
// Created by Kyowon Jeong on 9/25/19.
//

#include "OpenMS/ANALYSIS/TOPDOWN/SpectrumDeconvolution.h"

namespace OpenMS
{

  // constructor
  SpectrumDeconvolution::SpectrumDeconvolution(MSSpectrum &s, Parameter &p)
  {
    spec = s;
    param = p;
    //averagines = a;
    setFilters();
    updateLogMzPeaks();
  }

  /// default destructor
  SpectrumDeconvolution::~SpectrumDeconvolution()
  {
    delete[] binOffsets;
    delete[] filter;

    for (Size k = 0; k < param.hCharges.size(); k++)
    {
      delete[] hBinOffsets[k];
      delete[] harmonicFilter[k];
    }
    delete[] hBinOffsets;
    delete[] harmonicFilter;
  }

  void SpectrumDeconvolution::setFilters()
  {

    filter = new double[param.chargeRange];
    harmonicFilter = new double *[param.hCharges.size()];
//    int *random = new int[param.chargeRange];
    //std::srand(std::time(nullptr));

    for (int i = 0; i < param.chargeRange; i++)
    {
      filter[i] = log(
          1.0 / (i + param.minCharge));
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
  }

  void SpectrumDeconvolution::updateLogMzPeaks()
  {
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
  }


  double SpectrumDeconvolution::getBinValue(Size bin, double minV, double binWidth)
  {
    return minV + bin / binWidth;
  }

  Size SpectrumDeconvolution::getBinNumber(double v, double minV, double binWidth)
  {
    if (v < minV)
    {
      return 0;
    }
    return (Size) (((v - minV) * binWidth) + .5);

  }

  void SpectrumDeconvolution::updateMzBins(double &mzBinMinValue, Size &binNumber, double binWidth,
                                           float *mzBinIntensities
  )
  {
    mzBins = boost::dynamic_bitset<>(binNumber);
    // double *intensities = new double[binNumber];
    std::fill_n(mzBinIntensities, binNumber, .0);

    for (auto &p : logMzPeaks)
    {
      Size bi = getBinNumber(p.logMz, mzBinMinValue, binWidth);
      if (bi >= binNumber)
      {
        continue;
      }
      mzBins.set(bi);
      mzBinIntensities[bi] += p.intensity;

      auto delta = (p.logMz - getBinValue(bi, mzBinMinValue, binWidth));

      if (delta > 0)
      {
        if (bi < binNumber - 1)
        {
          mzBins.set(bi + 1);
          mzBinIntensities[bi + 1] += p.intensity;
        }
      }
      else if (delta < 0)
      {
        if (bi > 0)
        {
          mzBins.set(bi - 1);
          mzBinIntensities[bi - 1] += p.intensity;
        }
      }
    }
  }

  void SpectrumDeconvolution::unionPrevMassBins(double &massBinMinValue,
                                                std::vector<std::vector<Size>> &prevMassBinVector,
                                                std::vector<double> &prevMassBinMinValue)
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
      long shift = (long) (round((massBinMinValue - prevMassBinMinValue[i]) * param.binWidth));

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


  boost::dynamic_bitset<> SpectrumDeconvolution::getCandidateMassBinsForThisSpectrum(float *massIntensitites,
                                                                                     float *mzIntensities)
  {
    int chargeRange = param.chargeRange;
    int hChargeSize = (int) param.hCharges.size();
    int minContinuousChargePeakCount = param.minContinuousChargePeakCount;

    long binEnd = (long) massBins.size();
    auto candidateMassBinsForThisSpectrum = boost::dynamic_bitset<>(massBins.size());

    Byte *continuousChargePeakPairCount = new Byte[massBins.size()];
    std::fill_n(continuousChargePeakPairCount, massBins.size(), 0);

    Byte *prevCharges = new Byte[massBins.size()];
    std::fill_n(prevCharges, massBins.size(), (Byte) (chargeRange + 2));

    float *prevIntensities = new float[massBins.size()];
    std::fill_n(prevIntensities, massBins.size(), 1);

    std::fill_n(massIntensitites, massBins.size(), 0);
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
      auto &intensity = mzIntensities[mzBinIndex];
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

        if (prevCharges[massBinIndex] < chargeRange && cd != 1 && id < factor)
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
              auto &hintensity = mzIntensities[hmzBinIndex];
              if (hintensity > minInt
                  && hintensity < factor * maxInt
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
            continuousChargePeakPairCount[massBinIndex] = 0;
          }
          else
          {
            massIntensitites[massBinIndex] += intensity;
            if (!candidateMassBinsForThisSpectrum[massBinIndex])
            {
              candidateMassBinsForThisSpectrum[massBinIndex] =
                  ++continuousChargePeakPairCount[massBinIndex] >= minContinuousChargePeakCount;
            }
          }
        }
        prevIntensity = intensity;
        prevCharges[massBinIndex] = j;
      }
      mzBinIndex = mzBins.find_next(mzBinIndex);
    }
    auto mindex = candidateMassBinsForThisSpectrum.find_first();
    while (mindex != candidateMassBinsForThisSpectrum.npos)
    {
      auto &s = massIntensitites[mindex];
      // auto msnr = s / (noise[0][mindex]);
      float maxNoise = .0;
      for (auto k = 0; k < hChargeSize + 1; k++)
      {
        maxNoise = std::max(maxNoise, noise[k][mindex]);
        // msnr = min(snr, msnr);
      }
      //s -= maxNoise;
      s -= maxNoise;
      mindex = candidateMassBinsForThisSpectrum.find_next(mindex);
    }

    delete[] prevIntensities;
    delete[] continuousChargePeakPairCount;
    delete[] prevCharges;
    for (auto k = 0; k < hChargeSize + 1; k++)
    {
      delete[] noise[k];
    }
    delete[] noise;
    return candidateMassBinsForThisSpectrum;
  }


  Byte **SpectrumDeconvolution::updateMassBins_(boost::dynamic_bitset<> candidateMassBinsForThisSpectrum,
                                                float *massIntensities,
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

    auto mzBinIndex = mzBins.find_first();
    long binSize = (long) massBins.size();

    massBinsForThisSpectrum = boost::dynamic_bitset<>(massBins.size());
    auto toSkip = (candidateMassBinsForThisSpectrum | massBins).flip();
    massBins.reset();
    //massBinsForThisSpectrum.reset();

    while (mzBinIndex != mzBins.npos)
    {
      long maxIndex = -1;
      float maxCount = -1e11;
      Byte charge = 0;

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
          maxChargeRanges[maxIndex] = std::max(maxChargeRanges[maxIndex], charge);
          minChargeRanges[maxIndex] = std::min(minChargeRanges[maxIndex], charge);

          massBinsForThisSpectrum[maxIndex] = candidateMassBinsForThisSpectrum[maxIndex];

          mzChargeRanges[mzBinIndex] = charge;//minChargeRanges[maxIndex];//...
          massBins[maxIndex] = true;
        }
        // }
      }
      //}
      mzBinIndex = mzBins.find_next(mzBinIndex);
    }

    Byte **chargeRanges = new Byte *[3];
    chargeRanges[0] = minChargeRanges;
    chargeRanges[1] = maxChargeRanges;
    chargeRanges[2] = mzChargeRanges;
    // delete[] selected;
    return chargeRanges;
  }


  Byte **SpectrumDeconvolution::updateMassBins(double &massBinMinValue,
                                               float *massIntensities,
                                               float *mzIntensities
  )
  {
    long binThresholdMinMass = (long) getBinNumber(log(param.minMass), massBinMinValue, param.binWidth);
    long binThresholdMaxMass = (long) std::min(massBins.size(),
                                               1 + getBinNumber(log(param.maxMass), massBinMinValue, param.binWidth));
    auto candidateMassBins = getCandidateMassBinsForThisSpectrum(massIntensities, mzIntensities);

    auto perMassChargeRanges = updateMassBins_(candidateMassBins,
                                               massIntensities,
                                               binThresholdMinMass,
                                               binThresholdMaxMass);
    return perMassChargeRanges;
  }


  void SpectrumDeconvolution::getCandidatePeakGroups(double &mzBinMinValue, double &massBinMinValue,
                                                     float *massIntensities,
                                                     Byte **chargeRanges
  )
  {
    double binWidth = param.binWidth;
    double tol = param.tolerance * 2;
    int minCharge = param.minCharge;
    int chargeRange = param.chargeRange;
    int maxIsotopeCount = param.maxIsotopeCount;

    int logMzPeakSize = (int) logMzPeaks.size();
    Size massBinSize = massBins.size();
    int *currentPeakIndex = new int[param.chargeRange];
    std::fill_n(currentPeakIndex, param.chargeRange, 0);

    peakGroups.reserve(massBins.count());
    auto &minChargeRanges = chargeRanges[0];
    auto &maxChargeRanges = chargeRanges[1];
    auto &mzChargeRanges = chargeRanges[2];

    auto massBinIndex = massBins.find_first();
    // Size lastSetMassBinIndex = unionedMassBins.size();
    Size *peakBinNumbers = new Size[logMzPeakSize];
    for (int i = 0; i < logMzPeakSize; i++)
    {
      peakBinNumbers[i] = getBinNumber(logMzPeaks[i].logMz, mzBinMinValue, binWidth);
    }

    while (massBinIndex != massBins.npos)
    {
      double logM = getBinValue(massBinIndex, massBinMinValue, binWidth);
      double diff = Constants::C13C12_MASSDIFF_U / exp(logM);
      double isoLogM1 = logM - diff;
      double isoLogM2 = logM + diff;

      auto b1 = getBinNumber(isoLogM1, massBinMinValue, binWidth);
      if (b1 > 0)
      {
        if (massIntensities[massBinIndex] < massIntensities[b1])
        {
          massBinIndex = massBins.find_next(massBinIndex);
          continue;
        }
      }

      auto b2 = getBinNumber(isoLogM2, massBinMinValue, binWidth);
      if (b2 < massBins.size())
      {
        if (massIntensities[massBinIndex] < massIntensities[b2])
        {
          massBinIndex = massBins.find_next(massBinIndex);
          continue;
        }
      }

      if (massIntensities[b1] == 0 && massIntensities[b2] == 0)
      {
        massBinIndex = massBins.find_next(massBinIndex);
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
            auto intensity = logMzPeaks[cpi].intensity;
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
          const double mz = logMzPeaks[maxIntensityPeakIndex].mz;
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
                const double observedMz = logMzPeaks[peakIndex].mz;
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

                  if (bin < massBinSize) //
                  {
                    LogMzPeak p(logMzPeaks[peakIndex], charge, i * d);
                    pg.peaks.push_back(p);
                    lastPeakIndex = peakIndex;
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
                double delta = abs(centerMz - p.mz);

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
          peakGroups.push_back(pg);
        }
      }

      massBinIndex = massBins.find_next(massBinIndex);
    }
    delete[] currentPeakIndex;
    delete[] peakBinNumbers;
  }


  void SpectrumDeconvolution::removeOverlappingPeakGroups()
  { // pgs are sorted
    double tol = param.tolerance;

    std::vector<PeakGroup> merged;
    merged.reserve(peakGroups.size());

    for (Size i = 0; i < peakGroups.size(); i++)
    {
      bool select = true;
      auto &pg = peakGroups[i];
      double massTol = pg.monoisotopicMass * tol * 2;

      for (Size j = i + 1; j < peakGroups.size(); j++)
      {
        auto &pgo = peakGroups[j];
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
        auto &pgo = peakGroups[j];
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
    std::vector<PeakGroup>().swap(peakGroups);
    merged.swap(peakGroups);
    //vector<PeakGroup>().swap(merged);
  }


  void SpectrumDeconvolution::updatePerChargeIsotopeIntensity(
      double *perIsotopeIntensity,
      double *perChargeIntensity,
      PeakGroup &pg)
  {
    std::fill_n(perIsotopeIntensity, param.maxIsotopeCount, 0);
    std::fill_n(perChargeIntensity, param.chargeRange, 0);

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
      perIsotopeIntensity[p.isotopeIndex] += p.intensity;
      perChargeIntensity[index] += p.intensity;
    }
    pg.maxCharge = maxCharge;
    pg.minCharge = minCharge;

  }


  double SpectrumDeconvolution::getChargeFitScore(double *perChargeIntensity, int range)
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

    double th = maxPerChargeIntensity * .02;// as recommended in the original paper...
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
      }

      last = i;
    }
    if (last - first < 2)
    {
      return 0;
    }

    for (int i = first; i <= last; i++)
    {
      xs.push_back(i);
      ys.push_back((1 + perChargeIntensity[i]));
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
  }


  double
  SpectrumDeconvolution::getCosine(double *a,
                                   int &aStart,
                                   int &aEnd,
                                   IsotopeDistribution &b,
                                   int &bSize,
                                   double &bNorm,
                                   int offset)
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

  double SpectrumDeconvolution::getCosine(std::vector<double> &a, std::vector<double> &b, int off)
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

  double SpectrumDeconvolution::getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                                         double *perIsotopeIntensities,
                                                                         int perIsotopeIntensitiesSize,
                                                                         int &offset,
                                                                         FLASHDeconvHelperStructs::PrecalcularedAveragine &avg)
  {
    auto iso = avg.get(mass);
    auto isoNorm = avg.getNorm(mass);

    int isoSize = (int) iso.size();

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
      auto cos = SpectrumDeconvolution::getCosine(perIsotopeIntensities,
                                                  minIsotopeIndex,
                                                  maxIsotopeIndex,
                                                  iso,
                                                  isoSize,
                                                  isoNorm,
                                                  f);

      if (maxCosine <= cos)
      {
        maxCosine = cos;
        offset = f;
      }
    }

    return maxCosine;
  }


  bool SpectrumDeconvolution::checkChargeDistribution(double *perChargeIntensity,
                                                      int range,
                                                      int threshold)
  {
    double maxPerChargeIntensity = .0;
    int nonZeroStart = -1, nonZeroEnd = 0;
    for (int i = 0; i < range; i++)
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


  void SpectrumDeconvolution::scoreAndFilterPeakGroups(FLASHDeconvHelperStructs::PrecalcularedAveragine &avg)
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
        pg.updateMassesAndIntensity(avg);
        intensities.push_back(pg.intensity);
      }

      if (intensities.size() > mc)
      {
        sort(intensities.begin(), intensities.end());
        threshold = intensities[intensities.size() - mc];
      }
      std::vector<double>().swap(intensities);
    }

    auto perIsotopeIntensity = new double[param.maxIsotopeCount];
    auto perChargeIntensity = new double[param.chargeRange];

    for (auto &pg : peakGroups)
    {
      if (pg.intensity < threshold)
      {
        continue;
      }
      updatePerChargeIsotopeIntensity(
          perIsotopeIntensity, perChargeIntensity,
          pg);

      pg.chargeCosineScore = getChargeFitScore(perChargeIntensity, param.chargeRange);

      if (pg.peaks.empty() || (pg.chargeCosineScore < param.minChargeCosineSpec))
      {
        continue;
      }

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
                                                                       offset,
                                                                       avg);


      if (pg.peaks.empty() || (pg.isotopeCosineScore < param.minIsotopeCosineSpec))
      {
        continue;
      }


      pg.updateMassesAndIntensity(avg, offset, param.maxIsotopeCount);
      filteredPeakGroups.push_back(pg);

    }

    peakGroups.swap(filteredPeakGroups);
    std::vector<PeakGroup>().swap(filteredPeakGroups);


    removeOverlappingPeakGroups();
    //removeHarmonicPeakGroups(filteredPeakGroups, param);

    delete[] perIsotopeIntensity;
    delete[] perChargeIntensity;

  //  return filteredPeakGroups;
  }


  std::vector<SpectrumDeconvolution::PeakGroup> SpectrumDeconvolution::getPeakGroupsFromSpectrum(std::vector<std::vector<Size>> &prevMassBinVector,
                                                                                                std::vector<double> &prevMinBinLogMassVector,
                                                                                                FLASHDeconvHelperStructs::PrecalcularedAveragine &avg)
  {
    double massBinMaxValue = std::min(
        logMzPeaks[logMzPeaks.size() - 1].logMz -
        filter[param.chargeRange - param.minContinuousChargePeakCount - 1],
        log(param.maxMass));

    double massBinMinValue = logMzPeaks[0].logMz - filter[param.minContinuousChargePeakCount];
    double mzBinMinValue = logMzPeaks[0].logMz;
    double mzBinMaxValue = logMzPeaks[logMzPeaks.size() - 1].logMz;
    Size massBinNumber = getBinNumber(massBinMaxValue, massBinMinValue, param.binWidth) + 1;

    binOffsets = new long[param.chargeRange];

    for (int i = 0; i < param.chargeRange; i++)
    {
      binOffsets[i] = (long) round((mzBinMinValue - filter[i] - massBinMinValue) *
                                   param.binWidth);
    }

    hBinOffsets = new long *[param.hCharges.size()];
    for (Size k = 0; k < param.hCharges.size(); k++)
    {
      hBinOffsets[k] = new long[param.chargeRange];
      for (int i = 0; i < param.chargeRange; i++)
      {
        hBinOffsets[k][i] = (long) round((mzBinMinValue - harmonicFilter[k][i] - massBinMinValue) *
                                         param.binWidth);
      }
    }

    Size mzBinNumber = getBinNumber(mzBinMaxValue, mzBinMinValue, param.binWidth) + 1;
    float *mzBinIntensities = new float[mzBinNumber];

    updateMzBins(mzBinMinValue, mzBinNumber, param.binWidth, mzBinIntensities);

    float *massIntensities = new float[massBinNumber];

    massBins = boost::dynamic_bitset<>(massBinNumber);
    massBinsForThisSpectrum = boost::dynamic_bitset<>(massBinNumber);

    unionPrevMassBins(massBinMinValue, prevMassBinVector, prevMinBinLogMassVector);

    auto perMassChargeRanges = updateMassBins(massBinMinValue, massIntensities,
                                              mzBinIntensities);
    getCandidatePeakGroups(mzBinMinValue, massBinMinValue,
                           massIntensities,
                           perMassChargeRanges
    );
    scoreAndFilterPeakGroups(avg);
    if (prevMassBinVector.size() > 0 && prevMassBinVector.size() >= (Size) param.numOverlappedScans)
    {
      prevMassBinVector.erase(prevMassBinVector.begin());
      prevMinBinLogMassVector.erase(prevMinBinLogMassVector.begin());
    }

    std::vector<Size> mb;
    mb.reserve(peakGroups.size());
    for (auto &pg : peakGroups)//filteredPeakGroups
    {
      pg.peaks.shrink_to_fit();
      if (massBinsForThisSpectrum[pg.massBinIndex])
      {
        mb.push_back(pg.massBinIndex);
      }
    }

    prevMassBinVector.push_back(mb); //
    prevMinBinLogMassVector.push_back(massBinMinValue);

    prevMassBinVector.shrink_to_fit();
    prevMinBinLogMassVector.shrink_to_fit();

    for (int i = 0; i < 3; i++)
    {
      delete[] perMassChargeRanges[i]; // delete array within matrix
    }// delete actual matrix
    delete[] perMassChargeRanges;
    delete[] mzBinIntensities;
    delete[] massIntensities;
    return peakGroups;
  }
}
