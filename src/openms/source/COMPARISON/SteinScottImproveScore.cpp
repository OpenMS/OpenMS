// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SteinScottImproveScore.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
  /// default constructor
  SteinScottImproveScore::SteinScottImproveScore() :
    PeakSpectrumCompareFunctor()
  {
    setName("SteinScottImproveScore");
    defaults_.setValue("tolerance", 0.2, "defines the absolute error of the mass spectrometer");
    defaults_.setValue("threshold", 0.2, "if the calculated score is smaller than the threshold, a zero is given back");
    defaultsToParam_();
  }

  /// copy constructor
  SteinScottImproveScore::SteinScottImproveScore(const SteinScottImproveScore & source) = default;

  /// destructor
  SteinScottImproveScore::~SteinScottImproveScore() = default;

  /// assignment operator
  SteinScottImproveScore & SteinScottImproveScore::operator=(const SteinScottImproveScore & source)
  {
    if (this != &source)
    {
      PeakSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  /**
  @brief Similarity pairwise score itself

  This function return the similarity score of itself based on SteinScott.

  @param spec  const PeakSpectrum Spectrum 1
  @see SteinScottImproveScore()
  */
  double SteinScottImproveScore::operator()(const PeakSpectrum & spec) const
  {
    return operator()(spec, spec);
  }

  /**
  @brief Similarity pairwise score

  This function return the similarity score of two spectra based on SteinScott.

  @param s1  const PeakSpectrum Spectrum 1
  @param s2  const PeakSpectrum Spectrum 2
  @see SteinScottImproveScore()
  */
  double SteinScottImproveScore::operator()(const PeakSpectrum & s1, const PeakSpectrum & s2) const
  {
    const double epsilon = (double)param_.getValue("tolerance");
    const double constant = epsilon / 10000;

    //const double c(0.0004);
    double score(0), sum(0), sum1(0), sum2(0), sum3(0), sum4(0);
    /* std::cout << s1 << std::endl;
    std::cout << std::endl;
    std::cout << s2 << std::endl;*/

    for (const Peak1D& it1 : s1)
    {
      double temp = it1.getIntensity();
      sum1 += temp * temp;
      sum3 += temp;
    }

    for (const Peak1D& it1 : s2)
    {
      double temp = it1.getIntensity();
      sum2 += temp * temp;
      sum4 += temp;
    }
    double z = constant * (sum3 * sum4);
    Size j_left(0);
    for (Size i = 0; i != s1.size(); ++i)
    {
      for (Size j = j_left; j != s2.size(); ++j)
      {
        double pos1(s1[i].getMZ()), pos2(s2[j].getMZ());
        if (std::abs(pos1 - pos2) <= 2 * epsilon)
        {
          sum += s1[i].getIntensity() * s2[j].getIntensity();
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
    //std::cout<< sum << " Sum " << z << " z " << std::endl;
    score = (sum - z) / (std::sqrt((sum1 * sum2)));
    // std::cout<<score<< " score" << std::endl;
    if (score < (float)param_.getValue("threshold"))
    {
      score = 0;
    }

    return score;
  }

}
