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
  // Att0           1.5594
  // Att1           9.0293
  // Att2          -0.0391
  // Att3          -0.6219
  // Att4           2.7049
  // Intercept     -3.8496

  std::vector<double> Qscore::weight_CV_0_ {-1.5594, -9.0293, 0.0391, 0.6219, -2.7049, 3.8496};

  //====================================== CV 40
  // Att0              14.4268
  // Att1              10.6621
  // Att2               0.1747
  // Att3              -0.5034
  // Att4               2.1574
  // Intercept        -15.5248

  std::vector<double> Qscore::weight_CV_40_ {-14.4268, -10.6621, -0.1747, 0.5034, -2.1574, 15.5248};

  // ====================================== CV 50
  // Att0              14.8552
  // Att1              11.0904
  // Att2               0.2581
  // Att3              -0.5179
  // Att4               1.8335
  // Intercept        -15.7714

  std::vector<double> Qscore::weight_CV_50_ {-14.8552, -11.0904, -0.2581, 0.5179, -1.8335, 15.7714};

  //====================================== CV 60
  // Att0               16.5106
  // Att1                6.1827
  // Att2                0.3472
  // Att3                -0.814
  // Att4                1.0387
  // Intercept         -16.8964

  std::vector<double> Qscore::weight_CV_60_ {-16.5106, -6.1827, -0.3472, 0.814, -1.0387,  16.8964};

  //====================================== Normal
  // Att0             13.5141
  // Att1              9.2568
  // Att2              0.2671
  // Att3             -0.5303
  // Att4              1.8682
  // Intercept       -14.6228

  std::vector<double> Qscore::weight_centroid_ {-13.5141, -9.2568, -0.2671, 0.5303, -1.8682, 14.6228};
  std::vector<double> Qscore::weight_profile_ {-13.5141, -9.2568, -0.2671, 0.5303, -1.8682, 14.6228};


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
    std::vector<double> fvector(5, .0); // length of weights vector - 1, excluding the intercept weight.
    if (pg->empty())
      return fvector;
    int index = 0;
    fvector[index++] = pg->getIsotopeCosine(); // (log2(a + d));

    // a = pg->getChargeIsotopeCosine();
    fvector[index++] = pg->getIsotopeCosine() - pg->getChargeIsotopeCosine(pg->getRepAbsCharge()); // (log2(d + a / (d + a)));

    // a = pg->getSNR();
    fvector[index++] = log2(1 + pg->getSNR()); //(log2(d + a / (d + a)));

    // a = pg->getChargeScore();
    fvector[index++] = log2(1 + pg->getSNR()) - log2(1 + pg->getChargeSNR(pg->getRepAbsCharge())); //(log2(a + d));

    auto [z, Z] = pg->getAbsChargeRange();
    int min_i = -1, max_i = 0;
    for (const auto& p : *pg)
    {
      int i = p.isotopeIndex;
      max_i = std::max(max_i, i);
      if (min_i < 0)
        min_i = i;

      min_i = std::min(min_i, i);
    }

    auto used = std::vector<bool>((Z - z + 1) * (max_i - min_i + 1), false);
    for (const auto& p : *pg)
    {
      used[(p.abs_charge - z + 1) * (p.isotopeIndex - min_i + 1) - 1] = true;
    }
    int count = 0;
    for (const auto& b : used)
      if (b)
        count++;

    fvector[index++] = (double)count / used.size();
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
    std::sort(dspec.begin(), dspec.end(), [](const PeakGroup& p1, const PeakGroup& p2) { return p1.getIsotopeCosine() > p2.getIsotopeCosine(); });

    for (auto& pg : dspec)
    {
      int z = pg.getRepAbsCharge();
      per_charge_count[z]++;
      if (per_charge_count[z] > median * 2)
        continue;

      auto fv = toFeatureVector_(&pg);
      bool target = pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::target;
      // if (target && pg.getQscore2D() < .5) continue;
      for (auto& item : fv)
      {
        f << item << ",";
      }
      f << (target ? "T" : "F") << "\n";
    }
  }
} // namespace OpenMS
