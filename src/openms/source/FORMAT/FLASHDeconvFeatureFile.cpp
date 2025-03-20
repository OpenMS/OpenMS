// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FLASHDeconvFeatureFile.h>

namespace OpenMS
{
  /**
    @brief FLASHDeconv Spectrum level output *.tsv, *.msalign (for TopPIC) file formats
     @ingroup FileIO
**/

  void FLASHDeconvFeatureFile::writeHeader(std::fstream& fs)
  {
    fs << "FeatureIndex\tFileName\tMonoisotopicMass\tAverageMass\tMassCount\tStartRetentionTime"
          "\tEndRetentionTime\tRetentionTimeDuration\tApexRetentionTime"
          "\tSumIntensity\tMaxIntensity\tFeatureQuantity\tMinCharge\tMaxCharge\tChargeCount\tIsotopeCosineScore\tMaxQscore\tPerChargeIntensity\tPerIsotopeIntensity"
          "\n";
  }

  void FLASHDeconvFeatureFile::writePromexHeader(std::fstream& fs)
  {
    fs << "FeatureID\tMinScan\tMaxScan\tMinCharge\tMaxCharge\t"
          "MonoMass\tRepScan\tRepCharge\tRepMz\tAbundance\tApexScanNum\tApexIntensity\tMinElutionTime\tMaxElutionTime\t"
          "ElutionLength\tEnvelope\tLikelihoodRatio"
          "\n";
  }

  void FLASHDeconvFeatureFile::writeTopFDFeatureHeader(std::vector<std::fstream>& fs)
  {
    for (Size i = 0; i < fs.size(); ++i)
    {
      if (i == 0)
      {
        fs[i] << "Sample_ID\tID\tMass\tIntensity\tTime_begin\tTime_end\tTime_apex\tMinimum_charge_state\tMaximum_charge_state\tMinimum_fraction_id\tMaximum_fraction_id\n";
      }
      else
      {
        fs[i] << "Spec_ID\tFraction_ID\tFile_name\tScans\tMS_one_ID\tMS_one_scans\tPrecursor_mass\tPrecursor_intensity\tFraction_feature_ID\tFraction_feature_intensity\tFraction_feature_"
                 "score\tFraction_feature_time_apex\tSample_feature_ID\tSample_feature_intensity\n";
      }
    }
  }

  void FLASHDeconvFeatureFile::writeFeatures(const std::vector<FLASHDeconvHelperStructs::MassFeature>& mass_features, const String& file_name, std::fstream& fs)
  {
    int feature_index = 0;
    for (auto& mass_feature : mass_features)
    {
      auto mt = mass_feature.mt;
      double mass = mt.getCentroidMZ() + mass_feature.iso_offset * Constants::ISOTOPE_MASSDIFF_55K_U;
      double avg_mass = mass_feature.avg_mass;
      double sum_intensity = .0;

      for (auto& p : mt)
      {
        sum_intensity += p.getIntensity();
      }

      fs << feature_index++ << "\t" << file_name << "\t" << std::to_string(mass) << "\t" << std::to_string(avg_mass) << "\t" // massdiff
         << mt.getSize() << "\t" << mt.begin()->getRT() << "\t" << mt.rbegin()->getRT() << "\t" << mt.getTraceLength() << "\t" << mt[mt.findMaxByIntPeak()].getRT() << "\t" << sum_intensity << "\t"
         << mt.getMaxIntensity(false) << "\t" << mt.computePeakArea() << "\t" << mass_feature.min_charge << "\t" << mass_feature.max_charge << "\t" << mass_feature.charge_count << "\t"
         << mass_feature.isotope_score << "\t" << mass_feature.qscore << "\t";

      for (int i = mass_feature.min_charge; i <= mass_feature.max_charge; i++)
      {
        fs << mass_feature.per_charge_intensity[abs(i)];
        if (i < mass_feature.max_charge)
        {
          fs << ";";
        }
      }

      fs << "\t";
      int iso_end_index = 0;

      for (Size i = 0; i < mass_feature.per_isotope_intensity.size(); i++)
      {
        if (mass_feature.per_isotope_intensity[i] == 0)
        {
          continue;
        }
        iso_end_index = (int)i;
      }
      for (int i = 0; i <= iso_end_index; i++)
      {
        fs << mass_feature.per_isotope_intensity[i];

        if (i < iso_end_index)
        {
          fs << ";";
        }
      }
      fs << "\n";
    }
  }

  void FLASHDeconvFeatureFile::writeTopFDFeatures(const std::vector<FLASHDeconvHelperStructs::MassFeature>& mass_features, const std::map<int, PeakGroup>& precursor_peak_groups,
                                                  const std::map<int, double>& scan_rt_map, const String& file_name, std::vector<std::fstream>& fs)
  {
    int topid = 1;
    std::unordered_map<int, int> mtid_topid;

    for (Size l = 0; l < mass_features.size(); l++)
    {
      auto mass_feature = mass_features[l];
      double sum_intensity = .0;
      for (auto& m : mass_feature.mt)
      {
        sum_intensity += m.getIntensity();
      }
      for (Size i = 0; i < fs.size(); i++)
      {
        if (i == 0)
        {
          fs[i] << "0\t" << topid << "\t" << mass_feature.mt.getCentroidMZ() << "\t" << sum_intensity << "\t" << mass_feature.mt.begin()->getRT() << "\t" << mass_feature.mt.rbegin()->getRT() << "\t"
                << mass_feature.mt[mass_feature.mt.findMaxByIntPeak()].getRT() << "\t" << mass_feature.min_charge << "\t" << mass_feature.max_charge << "\t0\t0\n";
          mtid_topid[(int)l] = topid;
        }
      }
      topid++;
    }

    for (auto& precursor : precursor_peak_groups)
    {
      int ms1_scan_number = precursor.second.getScanNumber();
      int ms2_scan_number = precursor.first;
      double rt = scan_rt_map.at(ms2_scan_number);
      bool selected = false;
      int selected_index = -1;
      for (Size l = 0; l < mass_features.size(); l++)
      {
        auto mass_feature = mass_features[l];
        auto mt = mass_feature.mt;
        if (abs(precursor.second.getMonoMass() - mt.getCentroidMZ()) > 1.5)
        {
          continue;
        }
        if (rt < mt.begin()->getRT() || rt > mt.rbegin()->getRT())
        {
          continue;
        }
        selected = true;
        selected_index = (int)l;
        break;
      }

      if (selected)
      {
        for (Size i = 1; i < fs.size(); i++)
        {
          double sum_intensity = .0;
          for (auto& m : mass_features[selected_index].mt)
          {
            sum_intensity += m.getIntensity();
          }
          fs[i] << ms2_scan_number << "\t0\t" << file_name << "\t" << ms2_scan_number << "\t" << ms1_scan_number << "\t" << ms1_scan_number << "\t" << precursor.second.getMonoMass() << "\t"
                << precursor.second.getIntensity() << "\t" << mtid_topid[selected_index] << "\t" << sum_intensity << "\t-1000\t"
                << mass_features[selected_index].mt[mass_features[selected_index].mt.findMaxByIntPeak()].getRT() << "\t" << topid << "\t" << sum_intensity << "\n";
        }
        continue;
      }

      auto crange = precursor.second.getAbsChargeRange();

      for (Size i = 0; i < fs.size(); i++)
      {
        if (i == 0)
        {
          fs[i] << "0\t" << topid << "\t" << precursor.second.getMonoMass() << "\t" << precursor.second.getIntensity() << "\t" << rt - 1 << "\t" << rt + 1 << "\t" << rt << "\t"
                << (precursor.second.isPositive() ? std::get<0>(crange) : -std::get<1>(crange)) << "\t" << (precursor.second.isPositive() ? std::get<1>(crange) : -std::get<0>(crange)) << "\t0\t0\n";
        }
        else
        {
          fs[i] << ms2_scan_number << "\t0\t" << file_name << "\t" << ms2_scan_number << "\t" << ms1_scan_number << "\t" << ms1_scan_number << "\t" << precursor.second.getMonoMass() << "\t"
                << precursor.second.getIntensity() << "\t" << topid << "\t" << precursor.second.getIntensity() << "\t-1000\t" << rt << "\t" << topid << "\t" << precursor.second.getIntensity() << "\n";
        }
      }
      topid++;
    }
  }

  void FLASHDeconvFeatureFile::writePromexFeatures(const std::vector<FLASHDeconvHelperStructs::MassFeature>& mass_features, const std::map<int, PeakGroup>& precursor_peak_groups,
                                                   const std::map<int, double>& scan_rt_map, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, std::fstream& fs)
  {
    int promexid = 1;

    auto per_isotope_intensity = std::vector<double>(avg.getMaxIsotopeIndex(), .0);
    std::map<double, int> rt_scan_map;
    for (auto const& item : scan_rt_map)
    {
      rt_scan_map[item.second] = item.first;
    }

    for (auto& mass_feature : mass_features)
    {
      auto mt = mass_feature.mt;
      double sum_intensity = .0;

      int min_scan_num = -1;
      int max_scan_num = 0;

      for (auto& m : mt)
      {
        auto iter = rt_scan_map.lower_bound(m.getRT());
        if (iter != rt_scan_map.end())
        {
          int scan = iter->second;
          if (min_scan_num < 0)
          {
            min_scan_num = scan;
          }
          min_scan_num = std::min(min_scan_num, scan);
          max_scan_num = std::max(max_scan_num, scan);
        }
        sum_intensity += m.getIntensity();
      }

      fs << promexid << "\t" << min_scan_num << "\t" << max_scan_num << "\t" << mass_feature.min_charge << "\t" << mass_feature.max_charge << "\t" << std::to_string(mt.getCentroidMZ()) << "\t"
         << std::fixed << std::setprecision(2) << mass_feature.scan_number << "\t" << mass_feature.rep_charge << "\t" << mass_feature.rep_mz << "\t" << sum_intensity << "\t"
         << mass_feature.scan_number << "\t" << sum_intensity << "\t" << mt.begin()->getRT() / 60.0 << "\t" << mt.rbegin()->getRT() / 60.0 << "\t" << mt.getTraceLength() / 60.0 << "\t";

      int iso_end_index = 0;

      for (Size i = 0; i < mass_feature.per_isotope_intensity.size(); i++)
      {
        if (mass_feature.per_isotope_intensity[i] == 0)
        {
          continue;
        }
        iso_end_index = (int)i;
      }
      for (int i = 0; i <= iso_end_index; i++)
      {
        fs << i << "," << mass_feature.per_isotope_intensity[i];

        if (i < iso_end_index)
        {
          fs << ";";
        }
      }

      fs << "\t" << mass_feature.isotope_score << "\n";
      fs << std::setprecision(0);

      promexid++;
    }

    for (auto& precursor : precursor_peak_groups)
    {
      int ms2_scan_number = precursor.first;
      auto rt = scan_rt_map.at(ms2_scan_number);
      bool selected = false;
      double max_intensity = .0;

      for (auto& mass_feature : mass_features)
      {
        auto mt = mass_feature.mt;
        if (abs(precursor.second.getMonoMass() - mt.getCentroidMZ()) > 1.5)
        {
          continue;
        }
        if (rt < mt.begin()->getRT() || rt > mt.rbegin()->getRT())
        {
          continue;
        }

        double sum_intensity = .0;

        for (auto& p : mt)
        {
          sum_intensity += p.getIntensity();
        }

        if (max_intensity < sum_intensity)
        {
          max_intensity = sum_intensity;
        }
        selected = true;
      }

      if (selected)
      {
        continue;
      }

      auto crange = precursor.second.getAbsChargeRange();
      bool is_positive = precursor.second.isPositive();
      auto mz_range = precursor.second.getRepMzRange();
      double rmz = (std::get<0>(mz_range) + std::get<1>(mz_range)) / 2.0;

      for (auto& pp : precursor.second)
      {
        if (pp.isotopeIndex < 0 || pp.isotopeIndex >= (int)avg.getMaxIsotopeIndex())
        {
          continue;
        }
        per_isotope_intensity[pp.isotopeIndex] += pp.intensity;
      }

      fs << promexid << "\t" << precursor.second.getScanNumber() << "\t" << precursor.second.getScanNumber() << "\t" << (is_positive ? std::get<0>(crange) : -std::get<1>(crange)) << "\t"
         << (is_positive ? std::get<1>(crange) : -std::get<0>(crange)) << "\t" << std::to_string(precursor.second.getMonoMass()) << "\t" << std::fixed << std::setprecision(2)
         << precursor.second.getScanNumber() << "\t" << (is_positive ? precursor.second.getRepAbsCharge() : -precursor.second.getRepAbsCharge()) << "\t" << rmz << "\t"
         << precursor.second.getIntensity() << "\t" << precursor.second.getScanNumber() << "\t" << precursor.second.getIntensity() << "\t" << (rt - 1) / 60.0 << "\t" << (rt + 1) / 60.0 << "\t"
         << 2.0 / 60.0 << "\t";

      double isotope_score = precursor.second.getIsotopeCosine();

      for (size_t j = 0; j < avg.getMaxIsotopeIndex(); ++j)
      {
        if (per_isotope_intensity[j] <= 0)
        {
          continue;
        }
        fs << j << "," << per_isotope_intensity[j] << ";";
      }
      fs << "\t" << isotope_score << "\n";
      fs << std::setprecision(0);
      promexid++;
    }
  }
} // namespace OpenMS