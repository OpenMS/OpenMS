// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qscore.h>
#include <iomanip>

namespace OpenMS
{
//-28.0219, -0.1854, -0.1084, 0.0312, 0, 25.2087
  std::vector<double> Qscore::weight_centroid_ { -21.0476, 1.5045, -0.1303, 0.183, 0.1834, 17.804};
  // Att0                21.0476
  // Att1                -1.5045
  // Att2                 0.1303
  // Att3                 -0.183
  // Att4                -0.1834
  // Intercept           -17.804

  /// calculate Qscore using PeakGroup attributes
  double Qscore::getQscore(const PeakGroup* pg)
  {
    if (pg->empty())
    { // all zero
      return .0;
    }

    double score = weight_centroid_.back() + .5;
    auto fv = toFeatureVector_(pg);

    for (Size i = 0; i < weight_centroid_.size() - 1; i++)
    {
      score += fv[i] * weight_centroid_[i];
    }
    double qscore = 1.0 / (1.0 + exp(score));

    return qscore;
  }

  /// convert PeakGroup into feature (attribute) vector
  std::vector<double> Qscore::toFeatureVector_(const PeakGroup* pg)
  {
    std::vector<double> fvector(5, .0); // length of weights vector - 1, excluding the intercept weight.
    if (pg->empty())
      return fvector;
    int index = 0;
    fvector[index++] = pg->getIsotopeCosine(); // (log2(a + d));

    fvector[index++] = pg->getIsotopeCosine() - pg->getChargeIsotopeCosine(pg->getRepAbsCharge()); // (log2(d + a / (d + a)));

    fvector[index++] = log2(1 + pg->getChargeSNR(pg->getRepAbsCharge())); //(log2(d + a / (d + a)));

    fvector[index++] = log2(1 + pg->getChargeSNR(pg->getRepAbsCharge())) - log2(1 + pg->getSNR()); //(log2(a + d));

    fvector[index++] = pg->getAvgPPMError();

    return fvector;
  }

  /// to write down training csv file header.
  void Qscore::writeAttCsvForQscoreTrainingHeader(std::fstream& f)
  {
    PeakGroup pg;
    Size att_count = toFeatureVector_(&pg).size();
    for (Size i = 0; i < att_count; i++)
      f << "Att" << i << ",";
    f << "Class\n";
  }

  /// to write down training csv file rows
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