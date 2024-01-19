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
  std::vector<double> Qscore::weight_CV_0_ {-9.2654, 0.558, 0.014, 0.0249, 8.2536};

  //====================================== CV 40
  std::vector<double> Qscore::weight_CV_40_ {-22.48, 0.9741, -0.0082, -0.4234, 20.5086};

  // ====================================== CV 50
  std::vector<double> Qscore::weight_CV_50_ {-18.1896, -0.6636, -0.0075, -0.5833, 16.5708};

  //====================================== CV 60
  std::vector<double> Qscore::weight_CV_60_ {-22.1279, 0.342, -0.0484, -0.6139, 20.3207};

  //====================================== Normal
  std::vector<double> Qscore::weight_centroid_ {-22.7763, 0.9107, -0.1328, -0.4266, 22.2197}; // apr23 all
  std::vector<double> Qscore::weight_profile_(weight_centroid_);                              //{-6.0783, -1.3585, -0.1004, -0.3308, 5.9575}; // in silico profile


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

    double score = weights.back() + .5;
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
    std::vector<double> fvector(4, .0); // length of weights vector - 1, excluding the intercept weight.
    if (pg->empty())
      return fvector;
    int index = 0;
    fvector[index++] = pg->getIsotopeCosine(); // (log2(a + d));

    // a = pg->getSNR();
    fvector[index++] = pg->getIsotopeCosine() - pg->getChargeIsotopeCosine(pg->getRepAbsCharge()); // (log2(d + a / (d + a)));

    // a = pg->getChargeSNR();
    fvector[index++] = log2(1 + pg->getChargeSNR(pg->getRepAbsCharge())); //(log2(d + a / (d + a)));

    // a = pg->getChargeScore();
    fvector[index++] = log2(1 + pg->getChargeSNR(pg->getRepAbsCharge())) - log2(1 + pg->getSNR()); //(log2(a + d));

    return fvector;
  }

  void Qscore::writeAttCsvForQscoreTrainingHeader(std::fstream& f)
  {
    PeakGroup pg;
    Size att_count = toFeatureVector_(&pg).size();
    for (Size i = 0; i < att_count; i++)
      f << "Att" << i << ",";
    f << "Class\n";
  }

  void Qscore::writeAttCsvForQscoreTraining(const DeconvolvedSpectrum& deconvolved_spectrum, std::fstream& f)
  {
    DeconvolvedSpectrum dspec;
    dspec.reserve(deconvolved_spectrum.size());
    for (auto& pg : deconvolved_spectrum)
    {
      dspec.push_back(pg);
    }

    if (dspec.empty())
      return;

    for (auto& pg : dspec)
    {
      bool target = pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::target;
      auto fv = toFeatureVector_(&pg);

      for (auto& item : fv)
      {
        f << item << ",";
      }
      f << (target ? "T" : "F") << "\n";
    }
  }
} // namespace OpenMS
