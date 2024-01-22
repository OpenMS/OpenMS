// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
  std::vector<double> Qscore::weight_CV_0_ {-9.351, 0.6602, 0.0386, 0.052, 8.278};

  //====================================== CV 40
  std::vector<double> Qscore::weight_CV_40_ {-22.427, 0.8907, 0.0064, -0.4009, 20.4231};

  // ====================================== CV 50
  std::vector<double> Qscore::weight_CV_50_ {-17.9098, -0.6368, 0.0153, -0.5496, 16.2575};

  //====================================== CV 60
  std::vector<double> Qscore::weight_CV_60_ {-20.3958, 0.1383, -0.0095, -0.4995, 18.6557};

  //====================================== Normal
  std::vector<double> Qscore::weight_centroid_ {-22.3725, 0.9744, -0.1092, -0.3888, 20.7859}; // apr23 all
  std::vector<double> Qscore::weight_profile_ {-5.5533, 0.2767, 0.052, -0.027, 4.8713}; // in silico profile


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
