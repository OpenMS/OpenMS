// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/ZhangSimilarityScore.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <boost/math/special_functions/erf.hpp>

using namespace std;

namespace OpenMS
{
  ZhangSimilarityScore::ZhangSimilarityScore() :
    PeakSpectrumCompareFunctor()
  {
    setName("ZhangSimilarityScore");
    defaults_.setValue("tolerance", 0.2, "defines the absolute (in Da) or relative (in ppm) tolerance");
    defaults_.setValue("is_relative_tolerance", "false", "If set to true, the tolerance is interpreted as relative");
    defaults_.setValidStrings("is_relative_tolerance", {"true","false"});
    defaults_.setValue("use_linear_factor", "false", "if true, the intensities are weighted with the relative m/z difference");
    defaults_.setValidStrings("use_linear_factor", {"true","false"});
    defaults_.setValue("use_gaussian_factor", "false", "if true, the intensities are weighted with the relative m/z difference using a gaussian");
    defaults_.setValidStrings("use_gaussian_factor", {"true","false"});
    defaultsToParam_();
  }

  ZhangSimilarityScore::ZhangSimilarityScore(const ZhangSimilarityScore & source) = default;

  ZhangSimilarityScore::~ZhangSimilarityScore() = default;

  ZhangSimilarityScore & ZhangSimilarityScore::operator=(const ZhangSimilarityScore & source)
  {
    if (this != &source)
    {
      PeakSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double ZhangSimilarityScore::operator()(const PeakSpectrum & spec) const
  {
    return operator()(spec, spec);
  }

  double ZhangSimilarityScore::operator()(const PeakSpectrum & s1, const PeakSpectrum & s2) const
  {
    const double tolerance = (double)param_.getValue("tolerance");
    bool use_linear_factor = param_.getValue("use_linear_factor").toBool();
    bool use_gaussian_factor = param_.getValue("use_gaussian_factor").toBool();
    double score(0), sum(0), sum1(0), sum2(0) /*, squared_sum1(0), squared_sum2(0)*/;

    // TODO remove parameter 
    if (param_.getValue("is_relative_tolerance").toBool() )
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    for (const Peak1D& it1 : s1)
    {
      sum1 += it1.getIntensity();
      /*
for (PeakSpectrum::ConstIterator it2 = s1.begin(); it2 != s1.end(); ++it2)
{
  if (abs(it1->getPosition()[0] - it2->getPosition()[0]) <= 2 * tolerance)
  {
    squared_sum1 += it1->getIntensity() * it2->getIntensity();
  }
}*/
    }

/*
        UInt i_left(0);
        for (Size i = 0; i != s1.size(); ++i)
        {
            sum1 += s1[i].getIntensity();
            for (Size j = i_left; j != s1.size(); ++j)
            {
                double pos1(s1[i].getPosition()[0]), pos2(s1[j].getPosition()[0]);
                if (abs(pos1 - pos2) <= 2 * tolerance)
                {
                    squared_sum1 += s1[i].getIntensity() * s1[j].getIntensity();
                }
                else
                {
                    if (pos2 > pos1)
                    {
                        break;
                    }
                    else
                    {
                        i_left = i;
                    }
                }
            }
        }*/

/*
    i_left = 0;
    for (Size i = 0; i != s2.size(); ++i)
    {
      sum2 += s2[i].getIntensity();
      for (Size j = i_left; j != s2.size(); ++j)
      {
        double pos1(s2[i].getPosition()[0]), pos2(s2[j].getPosition()[0]);
        if (abs(pos1 - pos2) <= 2 * tolerance)
        {
          squared_sum1 += s2[i].getIntensity() * s2[j].getIntensity();
        }
        else
        {
          if (pos2 > pos1)
          {
            break;
          }
          else
          {
            i_left = i;
          }
        }
      }
    }*/

    for (const Peak1D& it1 : s2)
    {
      sum2 += it1.getIntensity();
      /*
for (PeakSpectrum::ConstIterator it2 = s2.begin(); it2 != s2.end(); ++it2)
{
  if (abs(it1->getPosition()[0] - it2->getPosition()[0]) <= 2 * tolerance)
  {
    squared_sum2 += it1->getIntensity() * it2->getIntensity();
  }
}
      */
    }

    Size j_left(0);
    for (Size i = 0; i != s1.size(); ++i)
    {
      for (Size j = j_left; j != s2.size(); ++j)
      {
        double pos1(s1[i].getMZ()), pos2(s2[j].getMZ());
        if (fabs(pos1 - pos2) < tolerance)
        {
          //double factor((tolerance - fabs(pos1 - pos2)) / tolerance);
          double factor = 1.0;

          if (use_linear_factor || use_gaussian_factor)
          {
            factor = getFactor_(tolerance, fabs(pos1 - pos2), use_gaussian_factor);
          }
          sum += sqrt(s1[i].getIntensity() * s2[j].getIntensity() * factor);
        }
        else
        {
          if (pos2 > pos1)
          {
            break;
          }
          else
          {
            j_left = j;
          }
        }
      }
    }


    /*
for (PeakSpectrum::ConstIterator it1 = s1.begin(); it1 != s1.end(); ++it1)
{
  for (PeakSpectrum::ConstIterator it2 = s2.begin(); it2 != s2.end(); ++it2)
  {
    if (abs(it1->getPosition()[0] - it2->getPosition()[0]) <= 2 * tolerance)
    {
      sum += sqrt(it1->getIntensity() * it2->getIntensity());
    }
  }
}*/

    score = sum / (sqrt(sum1 * sum2));

    return score;

  }

  double ZhangSimilarityScore::getFactor_(double mz_tolerance, double mz_difference, bool is_gaussian) const
  {
    double factor(0.0);

    if (is_gaussian)
    {
      static const double denominator = mz_tolerance * 3.0 * sqrt(2.0);
      factor = std::erfc(mz_difference / denominator);
      //cerr << "Factor: " << factor << " " << mz_tolerance << " " << mz_difference << endl;
    }
    else
    {
      factor = (mz_tolerance - mz_difference) / mz_tolerance;
    }
    return factor;
  }

}
