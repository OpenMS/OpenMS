// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>


namespace OpenMS
{
  MassFeatureTrace::MassFeatureTrace() : DefaultParamHandler("MassFeatureTrace")
  {
    Param mtd_defaults = MassTraceDetection().getDefaults();
    mtd_defaults.setValue("min_sample_rate", .1, "Minimum fraction of scans along the feature trace that must contain a peak. To raise feature detection sensitivity, lower this value close to 0.");
    mtd_defaults.setValue("min_trace_length", 10.0);

    mtd_defaults.setValue("chrom_peak_snr", .0);
    mtd_defaults.addTag("chrom_peak_snr", "advanced");
    mtd_defaults.setValue("reestimate_mt_sd", "false");
    mtd_defaults.addTag("reestimate_mt_sd", "advanced");
    mtd_defaults.setValue("noise_threshold_int", .0);
    mtd_defaults.addTag("noise_threshold_int", "advanced");

    mtd_defaults.setValue("quant_method", "area");
    mtd_defaults.addTag("quant_method", "advanced"); // hide entry

    defaults_.insert("", mtd_defaults);
    defaults_.setValue("min_isotope_cosine", .75, "cosine threshold between avg. and observed isotope pattern for MS1");
    defaultsToParam_();
  }

  std::vector<FLASHDeconvHelperStructs::MassFeature> MassFeatureTrace::findFeatures(const PrecalculatedAveragine& averagine)
  {
    MSExperiment map;
    std::map<int, MSSpectrum> index_spec_map;
    int min_abs_charge = INT_MAX;
    int max_abs_charge = INT_MIN;
    bool is_positive = true;
    std::vector<FLASHDeconvHelperStructs::MassFeature> mass_features;
    for (auto& item : peak_group_map_)
    {
      double rt = item.first;
      MSSpectrum deconv_spec;
      deconv_spec.setRT(rt);
      for (auto& pg : item.second)
      {
        is_positive = pg.second.isPositive();
        auto [z1, z2] = pg.second.getAbsChargeRange();
        max_abs_charge = max_abs_charge > z2 ? max_abs_charge : z2;
        min_abs_charge = min_abs_charge < z1 ? min_abs_charge : z1;

        Peak1D tp(pg.first, (float)pg.second.getIntensity());
        deconv_spec.push_back(tp);
      }
      map.addSpectrum(deconv_spec);
    }

    // when map size is less than 3, MassTraceDetection aborts - too few spectra for mass tracing.
    if (map.size() < 3)
    {
      return mass_features;
    }

    MassTraceDetection mtdet;
    Param mtd_param = getParameters().copy("");
    mtd_param.remove("min_isotope_cosine");

    mtdet.setParameters(mtd_param);
    std::vector<MassTrace> m_traces;

    mtdet.run(map, m_traces); // m_traces : output of this function
    int charge_range = max_abs_charge - min_abs_charge + 1;

    for (auto& mt : m_traces)
    {
      double max_qscore = .0;
      int min_feature_abs_charge = INT_MAX; // min feature charge
      int max_feature_abs_charge = INT_MIN; // max feature charge

      auto per_isotope_intensity = std::vector<float>(averagine.getMaxIsotopeIndex(), .0f);
      auto per_charge_intensity = std::vector<float>(charge_range + min_abs_charge + 1, .0f);

      double mass = mt.getCentroidMZ();
      double max_iso = 0;
      int max_iso_off = 0;
      boost::dynamic_bitset<> charges(charge_range + 1);
      std::vector<PeakGroup> pgs;
      pgs.reserve(mt.getSize());

      for (auto& p2 : mt)
      {
        auto& pg_map = peak_group_map_[p2.getRT()];
        auto& pg = pg_map[p2.getMZ()];
        auto [z1, z2] = pg.getAbsChargeRange();
        min_feature_abs_charge = min_feature_abs_charge < z1 ? min_feature_abs_charge : z1;
        max_feature_abs_charge = max_feature_abs_charge > z2 ? max_feature_abs_charge : z2;

        if (pg.getIsotopeCosine() > max_iso)
        {
          max_iso = pg.getIsotopeCosine();
        }
        pgs.push_back(pg);
      }

      for (auto& pg : pgs)
      {
        for (size_t z = min_abs_charge; z < per_charge_intensity.size(); z++)
        {
          float zint = pg.getChargeIntensity((int)z);
          if (zint <= 0)
          {
            continue;
          }
          charges[z - min_abs_charge] = true;
          per_charge_intensity[z] += zint;
        }

        int iso_off = int(.5 + (pg.getMonoMass() - mass) / pg.getIsotopeDaDistance());
        max_iso_off = std::max(max_iso_off, abs(iso_off));
        auto iso_int = pg.getIsotopeIntensities();
        for (size_t i = 0; i < per_isotope_intensity.size() - iso_off; i++)
        {
          if ((int)i + iso_off < 0 || i >= iso_int.size())
          {
            continue;
          }
          per_isotope_intensity[i + iso_off] += iso_int[i];
        }

        max_qscore = max_qscore < pg.getQscore() ? pg.getQscore() : max_qscore;
      }

      int offset = 0;
      double isotope_score = FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(mass, per_isotope_intensity, offset, averagine, 0, 0);

      if (isotope_score < min_isotope_cosine_)
      {
        continue;
      }

      FLASHDeconvHelperStructs::MassFeature mass_feature;
      mass_feature.iso_offset = offset;
      mass += offset * Constants::ISOTOPE_MASSDIFF_55K_U;

      mass_feature.avg_mass = averagine.getAverageMassDelta(mass) + mass;
      mass_feature.mt = mt;
      mass_feature.charge_count = (int)charges.count();
      mass_feature.isotope_score = isotope_score;
      mass_feature.min_charge = (is_positive ? min_feature_abs_charge : -max_feature_abs_charge);
      mass_feature.max_charge = (is_positive ? max_feature_abs_charge : -min_feature_abs_charge);
      mass_feature.qscore = max_qscore;

      mass_feature.per_charge_intensity = per_charge_intensity;
      mass_feature.per_isotope_intensity = per_isotope_intensity;

      auto apex = mt[mt.findMaxByIntPeak()];
      auto& sub_pg_map = peak_group_map_[apex.getRT()];
      auto& rep_pg = sub_pg_map[apex.getMZ()];
      mass_feature.rep_mz = apex.getMZ();
      mass_feature.scan_number = rep_pg.getScanNumber();
      mass_feature.rep_charge = rep_pg.getRepAbsCharge();
      mass_features.push_back(mass_feature);
    }
    return mass_features;
  }

  void MassFeatureTrace::storeInformationFromDeconvolvedSpectrum(DeconvolvedSpectrum& deconvolved_spectrum)
  {
    double rt = deconvolved_spectrum.getOriginalSpectrum().getRT();
    if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != 1)
    {
      return;
    }
    else
    {
      peak_group_map_[rt] = std::map<double, PeakGroup>();
      auto& sub_pg_map = peak_group_map_[rt];
      for (auto& pg : deconvolved_spectrum)
      {
        sub_pg_map[pg.getMonoMass()] = pg;
      }
    }
  }

  void MassFeatureTrace::updateMembers_()
  {
    min_isotope_cosine_ = param_.getValue("min_isotope_cosine");
  }
} // namespace OpenMS
