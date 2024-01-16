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
  // Att0                       58.4269
  // Att1                        7.3944
  // Att2                            -0
  // Att3                       -0.8218
  // Att4                        6.9182
  // Att5                      -19.6314
  // Intercept                 -53.5639

  std::vector<double> Qscore::weight_CV_0_ {-58.4269, -7.3944, 0, 0.8218, -6.9182, 19.6314, 53.5639};

  //====================================== CV 40
  // Att0              15.6911
  // Att1               9.2787
  // Att2               0.1691
  // Att3              -0.4813
  // Att4               2.0662
  // Intercept        -16.6082

  std::vector<double> Qscore::weight_CV_40_ {-15.6911, -9.2787, -0.1691, 0.4813, -2.0662, 16.6082};

  // ====================================== CV 50
  // Att0               16.1422
  // Att1                 9.224
  // Att2                0.2562
  // Att3               -0.5149
  // Att4                 1.777
  // Intercept         -16.9095

  std::vector<double> Qscore::weight_CV_50_ {-16.1422, -9.224, -0.2562, 0.5149, -1.777, 16.9095};

  //====================================== CV 60
  // Att0               18.1732
  // Att1                5.3769
  // Att2                0.3479
  // Att3               -0.6145
  // Att4                0.9822
  // Intercept         -18.4009

  std::vector<double> Qscore::weight_CV_60_ {-18.1732, -5.3769, -0.3479, 0.6145, -0.9822,  18.4009};

  //====================================== Normal
  // Att0                      28.9649
  // Att1                      -5.8543
  // Att2                      -0.2382
  // Att3                        0.532
  // Intercept                -23.1794

  // Att0            10.9947
  // Att1            -2.7864
  // Att2            -0.2119
  // Att3             0.5438
  // Intercept       -8.2189

  std::vector<double> Qscore::weight_centroid_ {- 28.9649, 5.8543, 0.2382, -0.532, 23.1794};
  std::vector<double> Qscore::weight_profile_ {- 28.9649, 5.8543, 0.2382, -0.532, 23.1794};


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
    std::vector<double> fvector(4, .0); // length of weights vector - 1, excluding the intercept weight.
    if (pg->empty())
      return fvector;
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
    PeakGroup pg;
    Size att_count = toFeatureVector_(&pg).size();
    for (Size i = 0; i < att_count; i++)
      f << "Att" << i << ",";
    f << "Class\n";
  }

  void Qscore::writeAttCsvForQscoreTraining(const DeconvolvedSpectrum& deconvolved_spectrum, std::fstream& f, const std::vector<double>& weights)
  {
    std::map<int, int> per_charge_count;
    for (auto& pg : deconvolved_spectrum)
    {
      int z = pg.getRepAbsCharge();
      if (per_charge_count.find(z) == per_charge_count.end())
        per_charge_count[z] = 0;
      per_charge_count[z]++;
    }

    std::vector<int> counts;
    counts.reserve(per_charge_count.size());
    for (auto& [z, c] : per_charge_count)
    {
      counts.push_back(c);
      c = 0; // initialize
    }

    std::sort(counts.begin(), counts.end());
    int median = counts[counts.size() / 2];

    auto dspec = deconvolved_spectrum;
 /*   std::vector<PeakGroup> targets, filtered;
    std::vector<double> decoy_masses;
    dspec.sort();
    for (const auto& pg : dspec)
    {
      if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::target) targets.push_back(pg);
      else
      {
        if (pg.getTargetDecoyType() != PeakGroup::TargetDecoyType::noise_decoy) decoy_masses.push_back(pg.getMonoMass());
        filtered.push_back(pg);
      }
    }
    std::sort(decoy_masses.begin(), decoy_masses.end());

    for (const auto& pg : targets)
    {
      double delta = pg.getMonoMass() * 1e-5 * 2;
      auto upper = std::upper_bound(decoy_masses.begin(), decoy_masses.end(), pg.getMonoMass() + delta);
      bool exclude = false;
      while (!exclude)
      {
        if (upper != decoy_masses.end())
        {
          if (std::abs(*upper - pg.getMonoMass()) < delta)
          {
            exclude = true;
          }
          if (pg.getMonoMass() - *upper > delta)
          {
            break;
          }
        }
        if (upper == decoy_masses.begin())
        {
          break;
        }
        --upper;
      }
      if (exclude)
      {
        continue;
      }
      filtered.push_back(pg);
    }
    dspec.setPeakGroups(filtered);*/

    std::sort(dspec.begin(), dspec.end(), [](const PeakGroup& p1, const PeakGroup& p2) { return p1.getIsotopeCosine() > p2.getIsotopeCosine(); });
    std::vector<String> decoy_class  {"F", "F", "F"};
    double max_weight = *std::max_element(weights.begin(), weights.end() - 1);

    for (auto& pg : dspec)
    {
      bool target = pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::target;
      auto fv = toFeatureVector_(&pg);
      int w_index = 0;

      if (!target)
      {
        if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::noise_decoy) w_index = 1;
        else if(pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::isotope_decoy) w_index = 2;
        double w = weights.size() >= w_index? (weights[w_index] / max_weight) : 1.0;
        w = std::max(w, .1);
        if (w < 1.0)
        {
          double random_value = (double)std::rand()/RAND_MAX;
          if (random_value > w) continue;
        }
      }

      int z = pg.getRepAbsCharge();
      per_charge_count[z]++;
      if (per_charge_count[z] > median * 3)
        continue;
      for (auto& item : fv)
      {
        f << item << ",";
      }
      f << (target ? "T" : decoy_class[w_index]) << "\n";
    }
  }
} // namespace OpenMS
