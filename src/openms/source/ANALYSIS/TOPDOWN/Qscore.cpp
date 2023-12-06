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
  // IsotopeCosine                   42.2496 centroid + profile
  // ChargeCosine                      2.8767
  // MassSNR1                         -0.1523
  // ChargeSNR1                        0.4797
  // Intercept                       -45.5284

  //==================================== CV 0
  // IsotopeCosine                27.6125
  // ChargeCosine                  3.1915
  // MassSNR1                      0.1087
  // ChargeSNR1                   -0.5081
  // Intercept                   -26.9646

  std::vector<double> Qscore::weight_CV_0_ {-27.6125, -3.1915, -0.1087, 0.5081, 26.9646};

  //====================================== CV 40
  // IsotopeCosine                  24.3741
  // ChargeCosine                    4.2927
  // MassSNR1                        0.1893
  // ChargeSNR1                     -0.5764
  // Intercept                     -25.1341

  std::vector<double> Qscore::weight_CV_40_ {-24.3741, -4.2927, -0.1893, 0.5764, 25.1341};

  // ====================================== CV 50
  // IsotopeCosine                  28.0786
  // ChargeCosine                      3.89
  // MassSNR1                        0.1007
  // ChargeSNR1                     -0.5481
  // Intercept                     -27.8257

  std::vector<double> Qscore::weight_CV_50_ {-28.0786, -3.89, -0.1007, 0.5481, 27.8257};

  //====================================== CV 60
  // IsotopeCosine                  32.4126
  // ChargeCosine                     1.463
  // MassSNR1                        0.0073
  // ChargeSNR1                     -0.4006
  // Intercept                     -29.7409

  std::vector<double> Qscore::weight_CV_60_ {-32.4126, -41.463, -0.0073, 0.4006, 29.7409};

  std::vector<double> Qscore::weight_centroid_ {-42.2496, -2.8767, 0.1523, -0.4797, 45.5284};
  std::vector<double> Qscore::weight_profile_ {-42.2496, -2.8767, 0.1523, -0.4797, 45.5284};


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

  void Qscore::writeAttCsvFromDummyHeader(std::fstream& f)
  {
    f << "MSLevel,Cos,SNR,AvgPPMError,ChargeScore,Class\n";
  }

  void Qscore::writeAttCsvFromDummy(const DeconvolvedSpectrum& deconvolved_spectrum, std::fstream& f)
  {
    uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
    String cns[] = {"T", "D", "D", "D"};
    for (auto& pg : deconvolved_spectrum)
    {
      if (pg.getChargeSNR(pg.getRepAbsCharge()) < .5) // remove masses with too low SNRs - they act as outliers.
      {
        continue;
      }
      if (pg.getChargeIsotopeCosine(pg.getRepAbsCharge()) < .85) // remove masses with too low SNRs - they act as outliers.
      {
        continue;
      }

      auto fv = toFeatureVector_(&pg);
      f << ms_level << ",";
      for (auto& item : fv)
      {
        f << item << ",";
      }
      f << cns[pg.getTargetDecoyType()] << "\n";
    }
  }
} // namespace OpenMS
