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

#include "OpenMS/ANALYSIS/TOPDOWN/PeakGroupScoring.h"

namespace OpenMS
{
  //Once a peak group is defined it is scored using this class.
  PeakGroupScoring::PeakGroupScoring(std::vector<PeakGroup> &pg,
                                     int minCharge,
                                     int maxCharge,
                                     double minChargeCosine,
                                     DoubleList &tolerance,
                                     IntList &maxMassCount,
                                     IntList &minContinuousChargePeakCount,
                                     DoubleList &minIsotopeCosine)
      :
      peakGroups(pg), minCharge(minCharge), minChargeCosine(minChargeCosine), chargeRange(maxCharge - minCharge),
      maxMassCount(maxMassCount), tolerance(tolerance), minContinuousChargePeakCount(minContinuousChargePeakCount),
      minIsotopeCosine(minIsotopeCosine)
  {

    //defaultsToParam_();
  }

  PeakGroupScoring::~PeakGroupScoring()
  {
  }


  double PeakGroupScoring::getChargeFitScore(double *perChargeIntensity, int range)
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
  }

  double
  PeakGroupScoring::getCorrelation(const double *a,
                                   int &aStart,
                                   int &aEnd,
                                   IsotopeDistribution &b,
                                   int &bSize,
                                   int offset)
  {
    double n = .0, d1 = .0;
    double ma = .0, mb = .0;
    for (int j = aStart; j < aEnd; j++)
    {
      ma += a[j];
    }
    ma /= (aEnd - aStart);

    for (int j = 0; j < bSize; j++)
    {
      mb += b[j].getIntensity();
    }
    mb /= bSize;

    double bNorm = .0;
    for (int j = 0; j < bSize; j++)
    {
      bNorm += (b[j].getIntensity() - mb) * (b[j].getIntensity() - mb);
    }

    for (int j = aStart; j < aEnd; j++)
    {
      d1 += (a[j] - ma) * (a[j] - ma);
      int i = j - offset;
      if (i < 0 || i >= bSize)
      {
        continue;
      }
      n += (a[j] - ma) * (b[i].getIntensity() - mb); //
    }

    double d = (d1 * bNorm);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }


  double
  PeakGroupScoring::getCosine(const double *a,
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
      //c++;
      n += a[j] * b[i].getIntensity(); //
    }
    //if (c < 2)
    //{
    //   return 0;
    //}
    double d = (d1 * bNorm);
    if (d <= 0)
    {
      return 0;
    }
    return n / sqrt(d);
  }

  double PeakGroupScoring::getCosine(const double *a, double *b, Size size)
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


  double PeakGroupScoring::getCosine(std::vector<double> &a, std::vector<double> &b, int off)
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

  double PeakGroupScoring::getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                                    double *perIsotopeIntensities,
                                                                    int &offset,
                                                                    FLASHDeconvHelperStructs::PrecalculatedAveragine &avg)
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
      auto cos = PeakGroupScoring::getCosine(perIsotopeIntensities,
                                             minIsotopeIndex,
                                             maxIsotopeIndex,
                                             iso,
                                             isoSize,
                                             isoNorm,
                                             f);

      /* auto cos = PeakGroupScoring::getCorrelation(perIsotopeIntensities,
                                              minIsotopeIndex,
                                              maxIsotopeIndex,
                                              iso,
                                              isoSize,
                                              f);
 */
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


  bool PeakGroupScoring::checkChargeDistribution(double *perChargeIntensity)
  {
    double maxPerChargeIntensity = .0;
    int nonZeroStart = -1, nonZeroEnd = 0;
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
    double intensityThreshold = maxPerChargeIntensity / factor;//intensities[intensities.size()*95/100] / 5.0;
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
      //maxn_r = maxn_r < n_r? n_r : maxn_r;
      if (n_r >= minContinuousChargePeakCount[0])
      {
        return true;
      }
      prevCharge = k;
    }

    return false;
  }

  std::vector<PeakGroup> &PeakGroupScoring::scoreAndFilterPeakGroups(int &msLevel,
                                                                     FLASHDeconvHelperStructs::PrecalculatedAveragine &avg)
  {
    std::vector<PeakGroup> filteredPeakGroups;
    filteredPeakGroups.reserve(peakGroups.size());
    double threshold = .0;

    auto mc = maxMassCount.size() > msLevel - 1 ? maxMassCount[msLevel - 1] : -1;
    if (mc > 0)
    {
      std::vector<double> intensities;
      intensities.reserve(peakGroups.size());

      for (auto &pg : peakGroups)
      {
        pg.updateMassesAndIntensity(avg);
        intensities.push_back(pg.intensity);
      }

      if (intensities.size() > (Size) mc)
      {
        sort(intensities.begin(), intensities.end());
        threshold = intensities[intensities.size() - mc];
      }
      std::vector<double>().swap(intensities);
    }

    auto perIsotopeIntensity = new double[avg.maxIsotopeIndex];
    auto perChargeIntensity = new double[chargeRange];
    auto perChargeMaxIntensity = new double[chargeRange];

    for (auto &pg : peakGroups)
    {
      if (pg.intensity < threshold)
      {
        //delete[] pg.perChargeSNR;
        continue; //
      }

      auto indices = updatePerChargeIsotopeIntensity(
          perIsotopeIntensity, perChargeIntensity,
          avg.maxIsotopeIndex, pg);

      pg.chargeCosineScore = getChargeFitScore(perChargeIntensity, chargeRange);

      if (msLevel == 1)
      {
        if (pg.empty() ||
            pg.chargeCosineScore <= minChargeCosine)
        {
          // delete[] pg.perChargeSNR;
          continue;
        }

        bool isChargeWellDistributed = checkChargeDistribution(perChargeIntensity);

        if (!isChargeWellDistributed)
        {
          continue;
        }

      }
      else
      {
        if (pg.empty() || pg.chargeCosineScore < 0.1)//
        {
          continue;
        }
      }

      int offset = 0;
      pg.isotopeCosineScore = getIsotopeCosineAndDetermineIsotopeIndex(pg[0].getUnchargedMass(),
                                                                       perIsotopeIntensity,
                                                                       offset,
                                                                       avg);

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
      std::fill_n(perChargeMaxIntensity, chargeRange, 0);

      for (auto charge = pg.minCharge; charge <= pg.maxCharge; charge++)
      {
        int j = charge - minCharge;
        auto perIsotopeIntensities = new double[avg.maxIsotopeIndex];
        std::fill_n(perIsotopeIntensities, avg.maxIsotopeIndex, .0);

        int minIsotopeIndex = avg.maxIsotopeIndex;
        int maxIsotopeIndex = 0;

        double minMz = pg.monoisotopicMass * 2;
        double maxMz = 0;
        double maxIntensity = .0;
        double sumIntensity = .0;
        double sp = .0;

        for (auto &p:pg)
        {
          if (p.charge != charge)
          {
            continue;
          }
          //

          if (p.isotopeIndex > isoSize)
          {
            continue;
          }

          perIsotopeIntensities[p.isotopeIndex] += p.intensity;
          sumIntensity += p.intensity;
          minIsotopeIndex = minIsotopeIndex < p.isotopeIndex ? minIsotopeIndex : p.isotopeIndex;
          maxIsotopeIndex = maxIsotopeIndex < p.isotopeIndex ? p.isotopeIndex : maxIsotopeIndex;

          minMz = minMz < p.mz ? minMz : p.mz;
          maxMz = maxMz > p.mz ? maxMz : p.mz;
          if (maxIntensity < p.intensity)
          {
            maxIntensity = p.intensity;
            perChargeMaxIntensity[j] = maxIntensity;
          }
          // sp += p.intensity * p.intensity;
        }
        if (maxIntensity <= 0)
        {
          delete[] perIsotopeIntensities;
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

        auto cos = PeakGroupScoring::getCosine(perIsotopeIntensities,
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

        delete[] perIsotopeIntensities;
      }

      pg.totalSNR = totalSignal / totalNoise;
      pg.qScore = -10000;

      for (auto& item : pg.perChargeInfo)
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
      if(pg.maxQScoreMzStart > pg.maxQScoreMzEnd) {
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

    removeHarmonicPeakGroups(tolerance[msLevel - 1]); //
    removeOverlappingPeakGroups(tolerance[msLevel - 1]);

    //std::cout<< "* " << mc<<std::endl;
    //(param.currentMaxMassCount); //

    if (msLevel > 1)
    {
      filterPeakGroupsByIsotopeCosine(mc);
    }
    else
    {
      filterPeakGroupsByQScore(mc);
    }
    delete[] perIsotopeIntensity;
    delete[] perChargeIntensity;
    delete[] perChargeMaxIntensity;
    return peakGroups;
  }


  void PeakGroupScoring::filterPeakGroupsByIsotopeCosine(int currentMaxMassCount)
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

  void PeakGroupScoring::filterPeakGroupsByQScore(int currentMaxMassCount)
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


  void PeakGroupScoring::removeHarmonicPeakGroups(double tol)
  {
    sort(peakGroups.begin(), peakGroups.end());
    std::vector<PeakGroup> merged;
    merged.reserve(peakGroups.size());

    std::vector<double> masses;
    masses.reserve(peakGroups.size());
    for (auto &peakGroup : peakGroups)
    {
      masses.push_back(peakGroup.monoisotopicMass);
    }
    for (auto &pg : peakGroups)
    {
      bool select = true;
      for (int h = 2; h <= 3; h++)
      {
        for (int k = 0; k < 2; k++)
        {
          auto hmass = k == 0 ? pg.monoisotopicMass * h : pg.monoisotopicMass / h;
          double massTol = hmass * tol;

          auto iter = std::lower_bound(masses.begin(), masses.end(), hmass - massTol);
          Size j = iter - masses.begin();

          if (j >= 0 && j < peakGroups.size())
          {
            for (; j < peakGroups.size(); j++)
            {
              auto &pgo = peakGroups[j];
              if (hmass - pgo.monoisotopicMass > massTol)
              {
                continue;
              }

              if (!select || pgo.monoisotopicMass - hmass > massTol)
              {
                break;
              }
              select &= pg.totalSNR >= pgo.totalSNR;
            }
          }
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

  void PeakGroupScoring::removeOverlappingPeakGroups(double tol)
  { // pgs are sorted
    int isoLength = 1; // inclusive
    sort(peakGroups.begin(), peakGroups.end());
    std::vector<PeakGroup> filtered;
    filtered.reserve(peakGroups.size());

    for (Size i = 0; i < peakGroups.size(); i++)
    {
      bool select = true;
      auto &pg = peakGroups[i];

      if (pg.monoisotopicMass <= 0)
      {
        continue;
      }
      if (i > 0 && abs(pg.monoisotopicMass - peakGroups[i - 1].monoisotopicMass) < 1e-3)
      {
        continue;
      }

      double massTol = pg.monoisotopicMass * tol;

      int j = i + 1;
      for (int l = 0; l <= isoLength; l++)
      {
        auto off = Constants::ISOTOPE_MASSDIFF_55K_U * l;
        for (; j < peakGroups.size(); j++)
        {
          auto &pgo = peakGroups[j];

          if (l != 0 && pgo.monoisotopicMass - pg.monoisotopicMass < off - massTol)
          {
            continue;
          }

          if (!select || pgo.monoisotopicMass - pg.monoisotopicMass > off + massTol)
          {
            break;
          }
          select &= pg.isotopeCosineScore > pgo.isotopeCosineScore;
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

          if (l != 0 && pg.monoisotopicMass - pgo.monoisotopicMass < off - massTol)
          {
            continue;
          }

          if (!select || pg.monoisotopicMass - pgo.monoisotopicMass > off + massTol)
          {
            break;
          }
          select &= pg.isotopeCosineScore > pgo.isotopeCosineScore;
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

  std::vector<int> PeakGroupScoring::updatePerChargeIsotopeIntensity(
      double *perIsotopeIntensity,
      double *perChargeIntensity,
      int maxIsotopeCount,
      PeakGroup &pg)
  {
    std::fill_n(perIsotopeIntensity, maxIsotopeCount, 0);
    std::fill_n(perChargeIntensity, chargeRange, 0);
    /*for (int j = 0; j < param.currentChargeRange; ++j)
    {
      std::fill_n(intensityGrid[j], param.maxIsotopeIndex, 0);
    }

    for (int j = 0; j < param.maxIsotopeIndex; ++j)
    {
      std::fill_n(intensityGrid2[j], param.currentChargeRange, 0);
    }*/

    int minPgCharge = chargeRange + minCharge + 1;
    int maxPgCharge = minCharge - 1;
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
}