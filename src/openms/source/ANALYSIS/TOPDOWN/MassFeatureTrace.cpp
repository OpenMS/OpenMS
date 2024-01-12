// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>


namespace OpenMS
{
  MassFeatureTrace::MassFeatureTrace() : DefaultParamHandler("MassFeatureTrace")
  {
    Param mtd_defaults = MassTraceDetection().getDefaults();
    mtd_defaults.setValue("min_sample_rate", .1, "Minimum fraction of scans along the feature trace that must contain a peak. To raise feature detection sensitivity, lower this value close to 0.");
    mtd_defaults.setValue(
      "min_trace_length", 10.0,
      "Minimum expected length of a mass trace (in seconds). Only for MS1 (or minimum MS level in the dataset) feature tracing. For MSn, all traces are kept regardless of this value.");

    mtd_defaults.setValue("chrom_peak_snr", .0);
    mtd_defaults.addTag("chrom_peak_snr", "advanced");
    mtd_defaults.setValue("reestimate_mt_sd", "false");
    mtd_defaults.addTag("reestimate_mt_sd", "advanced");
    mtd_defaults.setValue("noise_threshold_int", .0);
    mtd_defaults.addTag("noise_threshold_int", "advanced");

    mtd_defaults.setValue("quant_method", "area");
    mtd_defaults.addTag("quant_method", "advanced"); // hide entry

    defaults_.insert("", mtd_defaults);
    defaults_.setValue("min_cos", .75, "Cosine similarity threshold between avg. and observed isotope pattern.");

    defaultsToParam_();
  }

  std::vector<FLASHDeconvHelperStructs::MassFeature> MassFeatureTrace::findFeaturesAndUpdateQscore2D(const PrecalculatedAveragine& averagine, std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                                                                                                     int ms_level, bool is_decoy)
  {
    static uint findex = 1;
    MSExperiment map;
    std::map<int, MSSpectrum> index_spec_map;
    int min_abs_charge = INT_MAX;
    int max_abs_charge = INT_MIN;
    bool is_positive = true;
    std::vector<FLASHDeconvHelperStructs::MassFeature> mass_features;
    std::map<double, Size> rt_index_map;

    std::map<int, int> prev_scans;
    int prev_scan = 0;
    for (Size i = 0; i < deconvolved_spectra.size(); i++)
    {
      auto deconvolved_spectrum = deconvolved_spectra[i];
      if (deconvolved_spectrum.empty() || is_decoy != deconvolved_spectrum.isDecoy())
        continue;
      if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != ms_level)
        continue;
      int scan = deconvolved_spectrum.getScanNumber();

      if (scan > prev_scan)
        prev_scans[scan] = prev_scan;

      prev_scan = scan;
      double rt = deconvolved_spectrum.getOriginalSpectrum().getRT();
      rt_index_map[rt] = i;
      MSSpectrum deconv_spec;
      deconv_spec.setRT(rt);
      for (auto& pg : deconvolved_spectrum)
      {
        is_positive = pg.isPositive();
        auto [z1, z2] = pg.getAbsChargeRange();
        max_abs_charge = max_abs_charge > z2 ? max_abs_charge : z2;
        min_abs_charge = min_abs_charge < z1 ? min_abs_charge : z1;

        Peak1D tp(pg.getMonoMass(), (float)pg.getIntensity());
        deconv_spec.push_back(tp);
      }
      map.addSpectrum(deconv_spec);
    }
    map.sortSpectra();
    // when map size is less than 3, MassTraceDetection aborts - too few spectra for mass tracing.
    if (map.size() < 3)
    {
      return mass_features;
    }

    MassTraceDetection mtdet;
    Param mtd_param = getParameters().copy("");
    mtd_param.remove("min_cos");
    mtdet.setParameters(mtd_param);
    std::vector<MassTrace> m_traces;

    mtdet.setLogType(ProgressLogger::NONE);
    mtdet.run(map, m_traces); // m_traces : output of this function
    int charge_range = max_abs_charge - min_abs_charge + 1;

    for (auto& mt : m_traces)
    {
      double qscore_2D = 1.0;
      double tmp_qscore_2D = 1.0;
      int min_feature_abs_charge = INT_MAX; // min feature charge
      int max_feature_abs_charge = INT_MIN; // max feature charge

      auto per_isotope_intensity = std::vector<float>(averagine.getMaxIsotopeIndex(), .0f);
      auto per_charge_intensity = std::vector<float>(charge_range + min_abs_charge + 1, .0f);

      double mass = mt.getCentroidMZ();

      boost::dynamic_bitset<> charges(charge_range + 1);
      std::vector<std::vector<PeakGroup>::iterator> pgs;
      pgs.reserve(mt.getSize());
      std::vector<double> qscores;

      prev_scan = 0;
      for (auto& p2 : mt)
      {
        auto& dspec = deconvolved_spectra[rt_index_map[p2.getRT()]];
        if (dspec.empty() || is_decoy != dspec.isDecoy())
          continue;
        PeakGroup comp;
        comp.setMonoisotopicMass(p2.getMZ() - 1e-7);
        auto pg = std::lower_bound(dspec.begin(), dspec.end(), comp);
        if (pg == dspec.end() || std::abs(pg->getMonoMass() - p2.getMZ()) > 1e-7)
          continue;

        auto [z1, z2] = pg->getAbsChargeRange();
        min_feature_abs_charge = min_feature_abs_charge < z1 ? min_feature_abs_charge : z1;
        max_feature_abs_charge = max_feature_abs_charge > z2 ? max_feature_abs_charge : z2;
        int scan = dspec.getScanNumber();
        if (prev_scan != 0 && (prev_scans[scan] <= prev_scan)) // only when consecutive scans are connected.
        {
          tmp_qscore_2D *= (1.0 - pg->getQscore());
        }
        else
        {
          tmp_qscore_2D = 1.0 - pg->getQscore();
        }
        qscore_2D = std::min(qscore_2D, tmp_qscore_2D);
        prev_scan = scan;
        pgs.push_back(pg);
      }
      qscore_2D = 1.0 - qscore_2D;
      for (auto& pg : pgs)
      {
        for (size_t z = min_abs_charge; z < per_charge_intensity.size(); z++)
        {
          float zint = pg->getChargeIntensity((int)z);
          if (zint <= 0)
          {
            continue;
          }
          charges[z - min_abs_charge] = true;
          per_charge_intensity[z] += zint;
        }
        int iso_off = int(.5 + (pg->getMonoMass() - mass) / pg->getIsotopeDaDistance());
        auto iso_int = pg->getIsotopeIntensities();
        for (int i = 0; i + iso_off < per_isotope_intensity.size(); i++)
        {
          if ((int)i + iso_off < 0 || i >= iso_int.size())
          {
            continue;
          }
          per_isotope_intensity[i + iso_off] += iso_int[i];
        }
      }

      int offset = 0;
      double isotope_score = SpectralDeconvolution::getIsotopeCosineAndDetermineIsotopeIndex(mass, per_isotope_intensity, offset, averagine, 0, 0);

      if (isotope_score < .5)
      {
        continue;
      }
      double max_int = 0;
      PeakGroup rep_pg = *pgs[0];
      for (auto& pg : pgs)
      {
        if (max_int <= pg->getIntensity())
        {
          rep_pg = *pg;
          max_int = pg->getIntensity();
        }

        pg->setFeatureIndex(findex);
        if (findex > 0)
          pg->setQscore2D(qscore_2D);
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
      mass_feature.qscore = qscore_2D;

      mass_feature.per_charge_intensity = per_charge_intensity;
      mass_feature.per_isotope_intensity = per_isotope_intensity;

      mass_feature.rep_mz = rep_pg.getMonoMass();
      mass_feature.scan_number = rep_pg.getScanNumber();
      mass_feature.rep_charge = rep_pg.getRepAbsCharge();
      mass_feature.index = findex;
      mass_feature.is_decoy = is_decoy;
      mass_feature.ms_level = ms_level;
      mass_features.push_back(mass_feature);
      findex++;
    }
    return mass_features;
  }

  void MassFeatureTrace::updateMembers_()
  {
  }
} // namespace OpenMS
