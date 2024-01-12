// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

  void FLASHDeconvFeatureFile::writeHeader(std::fstream& fs, bool report_decoy)
  {
    fs << "FeatureIndex\tFileName\tMSLevel";
    if (report_decoy) fs << "\tIsDecoy";

    fs << "\tMonoisotopicMass\tAverageMass\tMassCount\tStartRetentionTime"
          "\tEndRetentionTime\tRetentionTimeDuration\tApexRetentionTime"
          "\tSumIntensity\tMaxIntensity\tFeatureQuantity\tMinCharge\tMaxCharge\tChargeCount\tIsotopeCosineScore\tQscore2D\tPerChargeIntensity\tPerIsotopeIntensity"
          "\n";
  }

  void FLASHDeconvFeatureFile::writeTopFDFeatureHeader(std::fstream& fs, uint ms_level)
  {
    if (ms_level == 1)
    {
      fs << "Sample_ID\tID\tMass\tIntensity\tTime_begin\tTime_end\tTime_apex\tMinimum_charge_state\tMaximum_charge_state\tMinimum_fraction_id\tMaximum_fraction_id\n";
    }
    else
    {
      fs << "Spec_ID\tFraction_ID\tFile_name\tScans\tMS_one_ID\tMS_one_scans\tPrecursor_mass\tPrecursor_intensity\tFraction_feature_ID\tFraction_feature_intensity\tFraction_feature_"
               "score\tFraction_feature_time_apex\tSample_feature_ID\tSample_feature_intensity\n";
    }

  }

  void FLASHDeconvFeatureFile::writeFeatures(const std::vector<FLASHDeconvHelperStructs::MassFeature>& mass_features, const String& file_name, std::fstream& fs, bool report_decoy)
  {
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

      fs << mass_feature.index << "\t" << file_name << "\t" << mass_feature.ms_level;

      if (report_decoy)
      {
        fs << "\t" << (mass_feature.is_decoy? 1 : 0);
      }

      fs << "\t" << std::to_string(mass) << "\t" << std::to_string(avg_mass) << "\t" // massdiff
         << mt.getSize() << "\t" << mt.begin()->getRT() << "\t" << mt.rbegin()->getRT() << "\t" << mt.getTraceLength() << "\t" << mt[mt.findMaxByIntPeak()].getRT() << "\t" << sum_intensity << "\t"
         << mt.getMaxIntensity(false) << "\t" << mt.computePeakArea() << "\t" << mass_feature.min_charge << "\t" << mass_feature.max_charge << "\t" << mass_feature.charge_count << "\t"
         << mass_feature.isotope_score << "\t" << std::setprecision (15) << (mass_feature.qscore) << std::setprecision (-1) << "\t";

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
                                                  const std::map<int, double>& scan_rt_map, const String& file_name, std::fstream& fs, uint ms_level)
  {
    int topid = 1;
    std::unordered_map<int, int> mtid_topid;

    for (Size l = 0; l < mass_features.size(); l++)
    {
      auto mass_feature = mass_features[l];
      if (mass_feature.is_decoy) continue;
      double sum_intensity = .0;
      for (auto& m : mass_feature.mt)
      {
        sum_intensity += m.getIntensity();
      }

      if (ms_level == 1)
      {
        fs << "0\t" << topid << "\t" << mass_feature.mt.getCentroidMZ() << "\t" << sum_intensity << "\t" << mass_feature.mt.begin()->getRT() << "\t" << mass_feature.mt.rbegin()->getRT() << "\t"
              << mass_feature.mt[mass_feature.mt.findMaxByIntPeak()].getRT() << "\t" << mass_feature.min_charge << "\t" << mass_feature.max_charge << "\t0\t0\n";
        mtid_topid[(int)l] = topid;
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
        if (mass_feature.is_decoy) continue;
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
        if (ms_level > 1)
        {
          double sum_intensity = .0;
          for (auto& m : mass_features[selected_index].mt)
          {
            sum_intensity += m.getIntensity();
          }
          fs << ms2_scan_number << "\t0\t" << file_name << "\t" << ms2_scan_number << "\t" << ms1_scan_number << "\t" << ms1_scan_number << "\t" << precursor.second.getMonoMass() << "\t"
                << precursor.second.getIntensity() << "\t" << mtid_topid[selected_index] << "\t" << sum_intensity << "\t-1000\t"
                << mass_features[selected_index].mt[mass_features[selected_index].mt.findMaxByIntPeak()].getRT() << "\t" << topid << "\t" << sum_intensity << "\n";
        }
        continue;
      }

      auto crange = precursor.second.getAbsChargeRange();

      if (ms_level == 1)
      {
        fs << "0\t" << topid << "\t" << precursor.second.getMonoMass() << "\t" << precursor.second.getIntensity() << "\t" << rt - 1 << "\t" << rt + 1 << "\t" << rt << "\t"
              << (precursor.second.isPositive() ? std::get<0>(crange) : -std::get<1>(crange)) << "\t" << (precursor.second.isPositive() ? std::get<1>(crange) : -std::get<0>(crange)) << "\t0\t0\n";
      }
      else
      {
        fs << ms2_scan_number << "\t0\t" << file_name << "\t" << ms2_scan_number << "\t" << ms1_scan_number << "\t" << ms1_scan_number << "\t" << precursor.second.getMonoMass() << "\t"
              << precursor.second.getIntensity() << "\t" << topid << "\t" << precursor.second.getIntensity() << "\t-1000\t" << rt << "\t" << topid << "\t" << precursor.second.getIntensity() << "\n";
      }
      topid++;
    }
  }
} // namespace OpenMS