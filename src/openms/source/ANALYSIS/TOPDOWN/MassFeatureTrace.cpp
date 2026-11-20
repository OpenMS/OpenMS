// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <numeric>

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
    defaults_.setValue("iso_merge_tol_ppm", 1.0,
                       "Traces overlapping in RT whose masses differ by n isotopes (within this ppm tolerance + 3 mDa) are merged into one feature. "
                       "Set to 0 or negative to disable.");

    defaultsToParam_();
  }

  std::vector<FLASHHelperClasses::MassFeature> MassFeatureTrace::findFeaturesAndUpdateQscore2D(const PrecalculatedAveragine& averagine, std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                                                                                                     int ms_level, bool is_decoy)
  {
    static uint findex = 1;
    MSExperiment map;
    std::map<int, MSSpectrum> index_spec_map;
    int min_abs_charge = INT_MAX;
    int max_abs_charge = INT_MIN;
    bool is_positive = true;
    std::vector<FLASHHelperClasses::MassFeature> mass_features;
    std::map<double, Size> rt_index_map;

    std::map<int, int> prev_scans;
    int prev_scan = 0;
    for (Size i = 0; i < deconvolved_spectra.size(); i++)
    {
      auto deconvolved_spectrum = deconvolved_spectra[i];
      if (deconvolved_spectrum.empty())
        continue;
      if ((int)deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != ms_level)
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
        if (is_decoy && pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::target) continue;
        if (!is_decoy && pg.getTargetDecoyType() != PeakGroup::TargetDecoyType::target) continue;

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
    double cos_threshold = mtd_param.getValue("min_cos");
    double merge_tol_ppm = mtd_param.getValue("iso_merge_tol_ppm");
    mtd_param.remove("min_cos");
    mtd_param.remove("iso_merge_tol_ppm");
    mtdet.setParameters(mtd_param);
    std::vector<MassTrace> m_traces;

    mtdet.setLogType(ProgressLogger::NONE);
    mtdet.run(map, m_traces); // m_traces : output of this function
    int charge_range = max_abs_charge - min_abs_charge + 1;

    // peak group behind a trace point (same filter as above); dspec.end() if none
    auto peakGroupAt = [&](const Peak2D& p2, std::vector<PeakGroup>::iterator& out) {
      auto& dspec = deconvolved_spectra[rt_index_map[p2.getRT()]];
      out = dspec.end();
      if (dspec.empty()) return false;
      PeakGroup comp;
      comp.setMonoisotopicMass(p2.getMZ() - 1e-7);
      auto pg = std::lower_bound(dspec.begin(), dspec.end(), comp);
      if (pg == dspec.end() || std::abs(pg->getMonoMass() - p2.getMZ()) > 1e-7) return false;
      if (is_decoy && pg->getTargetDecoyType() == PeakGroup::TargetDecoyType::target) return false;
      if (!is_decoy && pg->getTargetDecoyType() != PeakGroup::TargetDecoyType::target) return false;
      out = pg;
      return true;
    };

    // The isotope index is decided per spectrum; when it flips between scans, one elution yields traces one isotope spacing apart
    // that the ppm-based tracer cannot join. Traces overlapping in RT (>= half of the shorter one) whose masses differ by n isotopes
    // are grouped (single linkage) and become one feature. The spacing is the one the peak groups were built with (decoys may differ).
    double d = Constants::ISOTOPE_MASSDIFF_55K_U;
    for (const auto& mt : m_traces)
    {
      std::vector<PeakGroup>::iterator pg;
      if (mt.getSize() > 0 && peakGroupAt(*mt.begin(), pg)) { d = pg->getIsotopeDaDistance(); break; }
    }
    std::vector<Size> parent(m_traces.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&parent](Size x) {
      while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
      return x;
    };
    if (merge_tol_ppm > 0)
    {
      std::vector<Size> order(m_traces.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](Size a, Size b) { return m_traces[a].getCentroidMZ() < m_traces[b].getCentroidMZ(); });
      std::vector<double> sorted_mass(order.size());
      for (Size a = 0; a < order.size(); a++) sorted_mass[a] = m_traces[order[a]].getCentroidMZ();
      for (Size a = 0; a < order.size(); a++)
      {
        const auto& ta = m_traces[order[a]];
        double ma = sorted_mass[a], sa = ta.begin()->getRT(), ea = ta.rbegin()->getRT();
        // residual tolerance: ppm + 3 mDa, capped at 15 mDa so that a deamidation pair (0.984 Da, 18.4 mDa off) is never merged
        double tol = std::min(ma * merge_tol_ppm * 1e-6 + .003, .015);
        for (Size b = std::lower_bound(sorted_mass.begin(), sorted_mass.end(), ma + .5 * d) - sorted_mass.begin(); b < order.size(); b++)
        {
          double delta = sorted_mass[b] - ma;
          if (delta > 3.5 * d) break;
          int n = (int)std::round(delta / d);
          if (std::abs(delta - n * d) > tol) continue;
          const auto& tb = m_traces[order[b]];
          double sb = tb.begin()->getRT(), eb = tb.rbegin()->getRT();
          if (std::min(ea, eb) - std::max(sa, sb) < .5 * std::min(ea - sa, eb - sb)) continue;
          parent[find(order[a])] = find(order[b]);
        }
      }
    }
    std::vector<std::vector<Size>> by_root(m_traces.size()), groups;
    for (Size i = 0; i < m_traces.size(); i++) by_root[find(i)].push_back(i);
    for (auto& members : by_root)
    {
      if (members.empty()) continue;
      double lo = INFINITY, hi = -INFINITY;
      for (Size ti : members) { lo = std::min(lo, m_traces[ti].getCentroidMZ()); hi = std::max(hi, m_traces[ti].getCentroidMZ()); }
      if (hi - lo > 3.5 * d) // a chain spanning more than three isotopes is not one proteoform: keep its traces separate
        for (Size ti : members) groups.push_back({ti});
      else
        groups.push_back(members);
    }

    for (const auto& members : groups)
    {
      const bool merged = members.size() > 1;
      std::vector<Peak2D> points; // all trace points of the group, in RT order
      for (Size ti : members)
        for (const auto& p2 : m_traces[ti]) points.push_back(p2);
      std::sort(points.begin(), points.end(), [](const Peak2D& x, const Peak2D& y) { return x.getRT() < y.getRT(); });

      // Peak groups per scan. Two labels in one scan are the same raw peaks in most cases: then only the more intense one counts;
      // if they share less than half of their peaks they are distinct signal and both count.
      std::vector<std::vector<PeakGroup>::iterator> pgs, all_pgs;
      std::vector<double> pg_rt; // RT of the scan of each entry of pgs
      for (const auto& p2 : points)
      {
        std::vector<PeakGroup>::iterator pg;
        if (!peakGroupAt(p2, pg)) continue;
        all_pgs.push_back(pg);
        if (!pgs.empty() && pgs.back()->getScanNumber() == pg->getScanNumber())
        {
          auto& prev = pgs.back();
          Size shared = 0;
          for (const auto& x : *pg)
            for (const auto& y : *prev)
              if (x.mz == y.mz) { shared++; break; }
          if (shared * 2 >= std::min(pg->size(), prev->size()))
          {
            if (pg->getIntensity() > prev->getIntensity()) prev = pg;
            continue;
          }
        }
        pgs.push_back(pg);
        pg_rt.push_back(p2.getRT());
      }
      if (pgs.empty()) continue;

      double qscore_2D = 1.0;
      double tmp_qscore_2D = 1.0;
      int min_feature_abs_charge = INT_MAX; // min feature charge
      int max_feature_abs_charge = INT_MIN; // max feature charge
      int min_scan_number = INT_MAX;
      int max_scan_number = INT_MIN;
      prev_scan = 0;
      for (auto& pg : pgs)
      {
        auto [z1, z2] = pg->getAbsChargeRange();
        min_feature_abs_charge = std::min(min_feature_abs_charge, z1);
        max_feature_abs_charge = std::max(max_feature_abs_charge, z2);
        int scan = pg->getScanNumber();
        min_scan_number = std::min(min_scan_number, scan);
        max_scan_number = std::max(max_scan_number, scan);
        if (scan == prev_scan) continue; // second (distinct) label of the same scan
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
      }
      qscore_2D = 1.0 - qscore_2D;

      PeakGroup rep_pg = **std::max_element(pgs.begin(), pgs.end(), [](const auto& x, const auto& y) { return x->getIntensity() < y->getIntensity(); });

      // Reference mass of the feature: a single trace keeps its centroid (unchanged behaviour). A merged group is placed at the isotope
      // label that carries the most intensity; member masses are then averaged in that label's frame.
      double mass = m_traces[members[0]].getCentroidMZ();
      if (merged)
      {
        std::map<int, double> label_intensity;
        for (auto& pg : pgs) label_intensity[(int)std::round((pg->getMonoMass() - rep_pg.getMonoMass()) / d)] += pg->getIntensity();
        int best_label = std::max_element(label_intensity.begin(), label_intensity.end(), [](const auto& x, const auto& y) { return x.second < y.second; })->first;
        double frame = rep_pg.getMonoMass() + best_label * d, num = .0, den = .0;
        for (auto& pg : pgs)
        {
          int off = (int)std::round((pg->getMonoMass() - frame) / d);
          num += pg->getIntensity() * (pg->getMonoMass() - off * d);
          den += pg->getIntensity();
        }
        if (den <= 0) continue;
        mass = num / den;
      }

      // isotope and charge intensities of all members, aligned to the feature mass
      auto per_isotope_intensity = std::vector<float>(averagine.getMaxIsotopeIndex(), .0f);
      auto per_charge_intensity = std::vector<float>(charge_range + min_abs_charge + 1, .0f);
      boost::dynamic_bitset<> charges(charge_range + 1);
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
        int iso_off = (int)std::round((pg->getMonoMass() - mass) / pg->getIsotopeDaDistance());
        const auto& iso_int = pg->getIsotopeIntensities();
        int i_end = std::min((int)iso_int.size(), (int)per_isotope_intensity.size() - iso_off);
        for (int i = std::max(0, -iso_off); i < i_end; i++) per_isotope_intensity[i + iso_off] += iso_int[i];
      }

      // per-scan isotope arrays hold isotope -1 in bin 0, hence the shift of 1; the cosine is the feature quality filter only (no re-labelling)
      int offset = 0;
      float isotope_score = SpectralDeconvolution::getIsotopeCosineAndIsoOffset(mass, per_isotope_intensity, offset, averagine, 1, 0, std::vector<double>{});
      if (isotope_score < cos_threshold)
      {
        continue;
      }

      for (auto& pg : all_pgs)
      {
        pg->setFeatureIndex(findex);
        if (findex > 0)
          pg->setQscore2D(qscore_2D);
      }

      MassTrace mt = m_traces[members[0]];
      if (merged) // one trace point per scan at the feature mass; two distinct labels of a scan contribute their sum
      {
        std::vector<Peak2D> pts;
        pts.reserve(pgs.size());
        for (Size k = 0; k < pgs.size(); k++)
        {
          if (!pts.empty() && pts.back().getRT() == pg_rt[k])
          {
            pts.back().setIntensity(pts.back().getIntensity() + pgs[k]->getIntensity());
            continue;
          }
          Peak2D q;
          q.setRT(pg_rt[k]);
          q.setMZ(mass);
          q.setIntensity(pgs[k]->getIntensity());
          pts.push_back(q);
        }
        mt = MassTrace(pts);
        mt.updateWeightedMeanMZ();
        mt.updateWeightedMeanRT();
        mt.setQuantMethod(m_traces[members[0]].getQuantMethod());
      }

      FLASHHelperClasses::MassFeature mass_feature;
      mass_feature.iso_offset = 0;

      mass_feature.avg_mass = averagine.getAverageMassDelta(mass) + mass;
      mass_feature.mt = mt;
      mass_feature.charge_count = (int)charges.count();
      mass_feature.isotope_score = isotope_score;
      mass_feature.min_charge = (is_positive ? min_feature_abs_charge : -max_feature_abs_charge);
      mass_feature.max_charge = (is_positive ? max_feature_abs_charge : -min_feature_abs_charge);
      mass_feature.qscore = qscore_2D;

      mass_feature.per_charge_intensity = per_charge_intensity;
      mass_feature.per_isotope_intensity = per_isotope_intensity;

      mass_feature.rep_mz = mass_feature.avg_mass / rep_pg.getRepAbsCharge();
      mass_feature.scan_number = rep_pg.getScanNumber();
      mass_feature.min_scan_number = min_scan_number;
      mass_feature.max_scan_number = max_scan_number;
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