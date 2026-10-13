// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Sam Wein $
// $Authors: Sam Wein $
// --------------------------------------------------------------------------
#include "RNAPrecursorCorrector.h"
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace OpenMS
{
  RNAPrecursorCorrector::RNAPrecursorCorrector(double tolerance_ppm, int min_abs_charge,
                                              int max_abs_charge, bool positive,
                                              double max_mass, double min_cos, double min_snr) :
    tolerance_ppm_(tolerance_ppm), min_cos_(min_cos), min_snr_(min_snr),
    min_abs_charge_(min_abs_charge), max_abs_charge_(max_abs_charge), positive_(positive)
  {
    if (!std::isfinite(tolerance_ppm) || tolerance_ppm <= 0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "FLASH tolerance_ppm must be finite and positive.");
    }
    if (min_abs_charge < 1 || max_abs_charge < min_abs_charge)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "FLASH charge magnitudes must satisfy 1 <= min_abs_charge <= max_abs_charge.");
    }
    if (!std::isfinite(min_cos))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "FLASH min_cos must be finite.");
    }
    if (!std::isfinite(min_snr))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "FLASH min_snr must be finite.");
    }
    Param p = detector_.getParameters();
    if (!std::isfinite(max_mass) || max_mass < double(p.getValue("min_mass")))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "max_mass must be finite and at least the FLASH minimum mass (50 Da).");
    }
    p.setValue("tol", DoubleList{tolerance_ppm, tolerance_ppm});
    p.setValue("min_charge", positive ? min_abs_charge : -max_abs_charge);
    p.setValue("max_charge", positive ? max_abs_charge : -min_abs_charge);
    p.setValue("max_mass", max_mass);
    p.setValue("min_cos", DoubleList{min_cos, min_cos});
    p.setValue("min_snr", DoubleList{min_snr, min_snr});
    detector_.setParameters(p);
    detector_.calculateAveragine(true);
  }

  RNAPrecursorCorrector::Result RNAPrecursorCorrector::match(
    const DeconvolvedSpectrum& groups, const Precursor& precursor) const
  {
    Result result;
    result.original_mz = result.corrected_mz = precursor.getMZ();
    result.original_charge = result.charge = precursor.getCharge();
    if (!std::isfinite(result.original_mz) || result.original_mz <= 0)
    {
      result.status = "invalid_precursor";
      return result;
    }
    const int reported_charge = std::abs(precursor.getCharge());
    const double tolerance = result.original_mz * tolerance_ppm_ * 1e-6;
    const double carrier = positive_ ? Constants::PROTON_MASS_U : -Constants::PROTON_MASS_U;
    std::vector<Result> matches;
    for (const auto& group : groups)
    {
      // A profile isotope may contain several samples. Take its closest sample;
      // do not count them as competing envelopes or alter their intensities.
      std::map<int, const FLASHHelperClasses::LogMzPeak*> closest;
      std::set<int> observed_mono;
      for (const auto& peak : group)
      {
        const int z = peak.abs_charge;
        if (z < min_abs_charge_ || z > max_abs_charge_ || peak.isotopeIndex < 0 ||
            (reported_charge != 0 && z != reported_charge))
        {
          continue;
        }
        const double cos = group.getChargeIsotopeCosine(z), snr = group.getChargeSNR(z);
        if (!std::isfinite(cos) || !std::isfinite(snr) || cos < min_cos_ || snr <= min_snr_)
        {
          continue;
        }
        if (peak.isotopeIndex == 0) observed_mono.insert(z);
        const double error = std::abs(peak.mz - result.original_mz);
        if (error > tolerance) continue;
        auto old = closest.find(z);
        if (old == closest.end() || error < std::abs(old->second->mz - result.original_mz))
        {
          closest[z] = &peak;
        }
      }
      for (const auto& [z, peak] : closest)
      {
        Result candidate = result;
        // Preserve signed input conventions; unsigned and unknown input keeps
        // the existing positive-magnitude convention.
        candidate.charge = result.original_charge < 0 ? -z : z;
        candidate.isotope = peak->isotopeIndex;
        // FLASH establishes the isotope index. Preserve the instrument's mass
        // precision instead of replacing it with an aggregate envelope mass,
        // which may drift when profile samples or several charge states overlap.
        candidate.corrected_mz = result.original_mz -
          candidate.isotope * Constants::C13C12_MASSDIFF_U / z;
        candidate.mono_mass = (candidate.corrected_mz - carrier) * z;
        candidate.cosine = group.getChargeIsotopeCosine(z);
        candidate.snr = group.getChargeSNR(z);
        candidate.mono_observed = observed_mono.contains(z);
        candidate.selected_error_ppm = (result.original_mz - peak->mz) / peak->mz * 1e6;
        if (!std::isfinite(candidate.corrected_mz) || candidate.corrected_mz <= 0) continue;
        candidate.status = candidate.isotope == 0 ? "monoisotopic" : "corrected";
        // Duplicate detections of the same mass/charge are one candidate.
        auto duplicate = std::find_if(matches.begin(), matches.end(), [&candidate, tolerance](const Result& old)
        {
          return old.charge == candidate.charge && old.isotope == candidate.isotope &&
                 std::abs(old.corrected_mz - candidate.corrected_mz) <= tolerance;
        });
        if (duplicate == matches.end()) matches.push_back(candidate);
        else if (candidate.cosine > duplicate->cosine) *duplicate = candidate;
      }
    }
    if (matches.size() == 1) return matches.front();
    if (matches.size() > 1) result.status = "ambiguous";
    return result;
  }

  std::vector<RNAPrecursorCorrector::Result> RNAPrecursorCorrector::correct(MSExperiment& experiment)
  {
    std::map<std::string, Size> ids;
    std::set<std::string> duplicate_ids;
    for (Size i = 0; i < experiment.size(); ++i)
    {
      if (!experiment[i].getNativeID().empty() && !ids.emplace(experiment[i].getNativeID(), i).second)
      {
        duplicate_ids.insert(experiment[i].getNativeID());
      }
    }
    std::map<Size, DeconvolvedSpectrum> cache;
    std::vector<Result> results;
    Size last_ms1 = experiment.size();
    for (Size i = 0; i < experiment.size(); ++i)
    {
      auto& spec = experiment[i];
      if (spec.getMSLevel() == 1)
      {
        last_ms1 = i;
        continue;
      }
      if (spec.getMSLevel() != 2) continue;
      Result result;
      result.spectrum = spec.getNativeID();
      result.status = "invalid_precursor";
      if (spec.getPrecursors().size() == 1)
      {
        auto& precursor = spec.getPrecursors()[0];
        result.original_mz = result.corrected_mz = precursor.getMZ();
        result.original_charge = result.charge = precursor.getCharge();
        Size parent = last_ms1;
        std::string parent_source = "preceding_ms1";
        if (precursor.metaValueExists("spectrum_ref"))
        {
          const std::string ref = precursor.getMetaValue("spectrum_ref").toString();
          auto found = ids.find(ref);
          if (found != ids.end() && !duplicate_ids.contains(ref) && found->second < i &&
              experiment[found->second].getMSLevel() == 1)
          {
            parent = found->second;
            parent_source = "spectrum_ref";
          }
          else parent_source = "invalid_reference_fallback";
        }
        if (parent == experiment.size()) result.status = "no_ms1";
        else
        {
          if (!cache.contains(parent))
          {
            auto ms1 = experiment[parent];
            ms1.sortByPosition();
            detector_.performSpectrumDeconvolution(ms1, int(parent), PeakGroup());
            auto groups = detector_.getDeconvolvedSpectrum();
            groups.setOriginalSpectrum(MSSpectrum()); // Cache envelopes without a second raw MS1 copy.
            cache.emplace(parent, std::move(groups));
          }
          result = match(cache.at(parent), precursor);
          if (result.accepted())
          {
            // Preserve the acquisition window even if its center equalled the
            // original selected ion. Only the selected-ion m/z is corrected.
            if (!precursor.metaValueExists("isolation window target m/z"))
            {
              precursor.setMetaValue("isolation window target m/z", precursor.getMZ());
            }
            precursor.setMZ(result.corrected_mz);
            precursor.setCharge(result.charge);
            if (precursor.metaValueExists("selected ion m/z"))
            {
              precursor.setMetaValue("selected ion m/z", result.corrected_mz);
            }
          }
        }
        result.spectrum = spec.getNativeID();
        result.parent = parent == experiment.size() ? "" : experiment[parent].getNativeID();
        result.parent_source = parent_source;
      }
      // Replace the owned annotations, including accepted-only fields, so an
      // unresolved rerun cannot retain a previous correction's m/z or quality.
      std::vector<std::string> keys;
      spec.getKeys(keys);
      for (const auto& key : keys)
      {
        if (key.starts_with("flash_precursor_")) spec.removeMetaValue(key);
      }
      spec.setMetaValue("flash_precursor_status", result.status);
      spec.setMetaValue("flash_precursor_original_mz", result.original_mz);
      spec.setMetaValue("flash_precursor_original_charge", result.original_charge);
      spec.setMetaValue("flash_precursor_parent", result.parent);
      spec.setMetaValue("flash_precursor_parent_source", result.parent_source);
      if (result.accepted())
      {
        spec.setMetaValue("flash_precursor_mz", result.corrected_mz);
        spec.setMetaValue("flash_precursor_mono_mass", result.mono_mass);
        spec.setMetaValue("flash_precursor_isotope", result.isotope);
        spec.setMetaValue("flash_precursor_cosine", result.cosine);
        spec.setMetaValue("flash_precursor_snr", result.snr);
        spec.setMetaValue("flash_precursor_mono_observed", result.mono_observed ? "true" : "false");
      }
      results.push_back(result);
    }
    return results;
  }

  void RNAPrecursorCorrector::storeReport(const std::string& filename, const std::vector<Result>& results)
  {
    TextFile report;
    report.addLine("spectrum\tparent\tparent_source\tstatus\toriginal_mz\tcorrected_mz\toriginal_charge\tcharge\tisotope\tmono_mass\tcosine\tsnr\tselected_error_ppm\tmono_observed");
    const auto quote = [](const std::string& value)
    {
      std::string escaped = "\"";
      for (char c : value)
      {
        escaped += c;
        if (c == '"') escaped += c;
      }
      return escaped + "\"";
    };
    for (const auto& r : results)
    {
      report.addLine(quote(r.spectrum) + "\t" + quote(r.parent) + "\t" + r.parent_source + "\t" + r.status + "\t" +
        StringUtils::toStr(r.original_mz) + "\t" + StringUtils::toStr(r.corrected_mz) + "\t" +
        StringUtils::toStr(r.original_charge) + "\t" + StringUtils::toStr(r.charge) + "\t" +
        StringUtils::toStr(r.isotope) + "\t" + StringUtils::toStr(r.mono_mass) + "\t" +
        StringUtils::toStr(r.cosine) + "\t" + StringUtils::toStr(r.snr) + "\t" +
        StringUtils::toStr(r.selected_error_ppm) + "\t" + (r.mono_observed ? "true" : "false"));
    }
    report.store(filename);
  }
}
