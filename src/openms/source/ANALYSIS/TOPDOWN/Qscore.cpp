// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qscore.h>
#include <iomanip>

namespace OpenMS
{
  //==================================== CV 0
  // Att0         -2.9763
  // Att1         -1.5465
  // Att2         -0.3062
  // Att3          0.2173
  // Intercept     4.0676

  std::vector<double> Qscore::weight_CV_0_ {-2.9763, -1.5465, -0.3062, 0.2173, 4.0676};

  //====================================== CV 40
  // Att0              15.6788
  // Att1              -3.4887
  // Att2              -0.2782
  // Att3               0.3794
  // Intercept        -11.3143

  std::vector<double> Qscore::weight_CV_40_ {-15.6788, 3.4887, 0.2782, -0.3794, 11.3143};

  // ====================================== CV 50
  // Att0                 21.7578
  // Att1                  -4.259
  // Att2                 -0.1171
  // Att3                  0.2773
  // Intercept           -16.2634

  std::vector<double> Qscore::weight_CV_50_ {-21.7578, 4.259, 0.1171, -0.2773, 16.2634};

  //====================================== CV 60
  // Att0              20.7225
  // Att1              -2.3573
  // Att2                -0.29
  // Att3               0.6051
  // Intercept        -17.4618

  std::vector<double> Qscore::weight_CV_60_ {-20.7225, 2.3573, 0.29, -0.6051, 17.4618};

  //====================================== Normal
  // Att0               17.7589
  // Att1               -3.1289
  // Att2                -0.113
  // Att3                0.3179
  // Intercept         -13.7447

  std::vector<double> Qscore::weight_centroid_ {-17.7589, 3.1289, 0.113, -0.3179, 13.7447};
  std::vector<double> Qscore::weight_profile_ {-17.7589, 3.1289, 0.113, -0.3179, 13.7447};


  double Qscore::getQscore(const PeakGroup* pg, const MSSpectrum& spectrum)
  {
    if (pg->empty())
    { // all zero
      return .0;
    }

    bool is_profile = spectrum.getType(false) != SpectrumSettings::CENTROID;
    auto filter_str = spectrum.getMetaValue("filter string").toString();
    Size pos = filter_str.find("cv=");
    double cv = 1;

    if (pos != String::npos)
    {
      Size end = filter_str.find(" ", pos);
      if (end == String::npos)
        end = filter_str.length() - 1;
      cv = std::stod(filter_str.substr(pos + 3, end - pos));
    }
    auto weights = (is_profile ? weight_profile_ : weight_centroid_);
    if (cv <= 0)
    {
      const std::vector<double> cvs {.0, -40.0, -50.0, -60.0};
      double min_val = cvs.back();
      for (double i : cvs)
      {
        double diff = std::abs(cv - i);
        if (diff > std::abs(cv - min_val))
          continue;
        min_val = i;
      }

      if (min_val == .0)
      {
        weights = weight_CV_0_;
      }
      else if (min_val == -40.0)
      {
        weights = weight_CV_40_;
      }
      else if (min_val == -50.0)
      {
        weights = weight_CV_50_;
      }
      else
      {
        weights = weight_CV_60_;
      }
    }

    double score = weights.back();
    auto fv = toFeatureVector_(pg);

    for (Size i = 0; i < weights.size() - 1; i++)
    {
      score += fv[i] * weights[i];
    }
    double qscore = 1.0 / (1.0 + exp(score));

    return qscore;
  }

  std::vector<double> Qscore::toFeatureVector_(const PeakGroup* pg)
  {
    std::vector<double> fvector(4); // length of weights vector - 1, excluding the intercept weight.

    int index = 0;
    fvector[index++] = pg->getIsotopeCosine(); // (log2(a + d));

    // a = pg->getSNR();
    fvector[index++] = pg->getChargeIsotopeCosine(pg->getRepAbsCharge()); // (log2(d + a / (d + a)));

    // a = pg->getAvgPPMError();
    fvector[index++] = log2(1 + pg->getSNR()); //(log2(d + a / (d + a)));

    // a = pg->getChargeScore();
    fvector[index++] = log2(1 + pg->getChargeSNR(pg->getRepAbsCharge())); //(log2(a + d));

    return fvector;
  }

  void Qscore::writeAttCsvForQscoreTrainingHeader(std::fstream& f)
  {
    Size att_count = weight_centroid_.size() - 1;
    for (Size i = 0; i < att_count; i++)
      f << "Att" << i << ",";
    f << "Class\n";
  }

  void Qscore::writeAttCsvForQscoreTraining(const DeconvolvedSpectrum& deconvolved_spectrum, std::fstream& f)
  {
    for (auto& pg : deconvolved_spectrum)
    {
      auto fv = toFeatureVector_(&pg);
      bool target = pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::target;
      for (auto& item : fv)
      {
        f << item << ",";
      }
      f << (target? "T" : "F") << "\n";
    }
  }
} // namespace OpenMS
