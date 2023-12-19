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
  // IsotopeCosine          16.7614
  // ChargeCosine            -0.8677
  // MassSNR1                 0.2665
  // ChargeSNR1              -0.3478
  // Intercept              -14.1557

  std::vector<double> Qscore::weight_CV_0_ {-16.7614, 0.8677, -0.2665, 0.3478, 14.1557};

  //====================================== CV 40
  // IsotopeCosine         13.9501
  // ChargeCosine            -0.705
  // MassSNR1                0.2562
  // ChargeSNR1             -0.3645
  // Intercept             -11.6675

  std::vector<double> Qscore::weight_CV_40_ {-13.9501, 0.705, -0.2562, 0.3645, 11.6675};

  // ====================================== CV 50
  // IsotopeCosine           19.0932
  // ChargeCosine             -0.4785
  // MassSNR1                  0.2784
  // ChargeSNR1               -0.3956
  // Intercept               -16.5475

  std::vector<double> Qscore::weight_CV_50_ {-19.0932, 0.4785, -0.2784, 0.3956, 16.5475};

  //====================================== CV 60
  // IsotopeCosine           20.2079
  // ChargeCosine             -1.4063
  // MassSNR1                  0.1948
  // ChargeSNR1               -0.2952
  // Intercept               -16.1873

  std::vector<double> Qscore::weight_CV_60_ {-20.2079, 1.4063, -0.1948, 0.2952, 16.1873};

  //====================================== Normal
  // IsotopeCosine        12.778
  // ChargeCosine         -1.6103
  // MassSNR1             -0.0726
  // ChargeSNR1           -0.0815
  // Intercept            -9.7626

  std::vector<double> Qscore::weight_centroid_ {-12.778, 1.6103, 0.0726, 0.0815, 9.7626};
  std::vector<double> Qscore::weight_profile_ {-12.778, 1.6103, 0.0726, 0.0815, 9.7626};


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
      for (auto& item : fv)
      {
        f << item << ",";
      }
      f << (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::target? "T" : "F") << "\n";
    }
  }
} // namespace OpenMS
