// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroupScoring.h>
#include <iomanip>

namespace OpenMS
  {
    std::vector<double> PeakGroupScoring::weight_ { -21.0476, 1.5045, -0.1303, 0.183, 0.1834, 17.804};
    // Att0                21.0476
    // Att1                -1.5045
    // Att2                 0.1303
    // Att3                 -0.183
    // Att4                -0.1834
    // Intercept           -17.804

    int PeakGroupScoring::getQscoreFeatureCount()
    {
      return weight_.size() - 1;
    }

    /// calculate PeakGroupScoring using PeakGroup attributes
    double PeakGroupScoring::getQscore(const PeakGroup* pg)
    {
      if (pg->empty())
      { // all zero
        return .0;
      }

      double score = weight_.back() + .5;
      auto fv = toFeatureVector(pg);

      for (Size i = 0; i < weight_.size() - 1; i++)
      {
        score += fv[i] * weight_[i];
      }
      double qscore = 1.0 / (1.0 + exp(score));

      return qscore;
    }

    /// convert PeakGroup into feature (attribute) vector
    std::vector<double> PeakGroupScoring::toFeatureVector(const PeakGroup* pg)
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
    void PeakGroupScoring::writeAttCsvForQscoreTrainingHeader(std::fstream& f)
    {
      PeakGroup pg;
      Size att_count = toFeatureVector(&pg).size();
      for (Size i = 0; i < att_count; i++)
        f << "Att" << i << ",";
      f << "Class\n";
    }

    /// to write down training csv file rows
    void PeakGroupScoring::writeAttCsvForQscoreTraining(const DeconvolvedSpectrum& deconvolved_spectrum, std::fstream& f)
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
        auto fv = toFeatureVector(&pg);

        for (auto& item : fv)
        {
          f << item << ",";
        }
        f << (target ? "T" : "F") << "\n";
      }
    }

    double PeakGroupScoring::getDLscore(PeakGroup* pg, const MSSpectrum& spec, const FLASHHelperClasses::PrecalculatedAveragine& avg, double tol)
    {
      const auto& [sig, noise]
                  = pg->getDLVector(spec, charge_count_for_DL_scoring_, iso_count_for_DL_scoring_, avg, tol);

      /// calculate score with sig and  noise


      return 0;
    }

  } // namespace OpenMS