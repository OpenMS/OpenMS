//
// Created by Kyowon Jeong on 10/8/19.
//

#include "OpenMS/ANALYSIS/TOPDOWN/PeakGroupScoring.h"
namespace OpenMS
{

  PeakGroupScoring::PeakGroupScoring(std::vector<PeakGroup> &pg, Parameter &p): peakGroups(pg), param(p)
  {
  }

  PeakGroupScoring::~PeakGroupScoring(){

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
  PeakGroupScoring::getCosine(double *a,
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
      auto cos = PeakGroupScoring::getCosine(perIsotopeIntensities,
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


  bool PeakGroupScoring::checkChargeDistribution(double *perChargeIntensity,
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

  std::vector<FLASHDeconvHelperStructs::PeakGroup> & PeakGroupScoring::scoreAndFilterPeakGroups(FLASHDeconvHelperStructs::PrecalcularedAveragine &avg)
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

    return peakGroups;
  }


  void PeakGroupScoring::removeOverlappingPeakGroups()
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


  void PeakGroupScoring::updatePerChargeIsotopeIntensity(
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

}
