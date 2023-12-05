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
  // ﻿IsotopeCosine                   42.2496
  // ChargeCosine                      2.8767
  // MassSNR1                         -0.1523
  // ChargeSNR1                        0.4797
  // Intercept                       -45.5284

 // std::vector<double> Qscore::weight_centroid_{-40.7425, -.205, 0.1984, -0.213, 40.3701};
  std::vector<double> Qscore::weight_profile_{-42.2496, -2.8767, 0.1523, -0.4797,  45.5284};
  //std::vector<double> Qscore::weight_CV_{-55.8387, 0.0253, 0.2473, -0.6765,  55.8594};

  double Qscore::getQscore(const PeakGroup* pg, bool is_profile, double cv)
  {
    if (pg->empty())
    { // all zero
      return .0;
    }

    auto weights = weight_profile_;//cv > 0 ? weight_CV_ : (is_profile? weight_profile_ : weight_centroid_);
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
