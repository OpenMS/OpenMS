// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Sam Wein $
// $Authors: Sam Wein $
// --------------------------------------------------------------------------
#include "RNADeisotoper.h"
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace OpenMS
{
  RNADeisotoper::RNADeisotoper(double tolerance_ppm, int max_abs_charge, bool positive,
                             double max_mass, double min_cos, double min_snr) :
    positive_(positive), max_abs_charge_(max_abs_charge), min_cos_(min_cos), min_snr_(min_snr)
  {
    if (!std::isfinite(tolerance_ppm) || tolerance_ppm <= 0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "FLASH tolerance_ppm must be finite and positive.");
    }
    if (max_abs_charge < 1)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "FLASH max_abs_charge must be positive.");
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
    p.setValue("min_charge", positive ? 1 : -max_abs_charge);
    p.setValue("max_charge", positive ? max_abs_charge : -1);
    p.setValue("max_mass", max_mass);
    p.setValue("min_cos", DoubleList{min_cos, min_cos});
    p.setValue("min_snr", DoubleList{min_snr, min_snr});
    detector_.setParameters(p);
    detector_.calculateAveragine(true);
  }

  void RNADeisotoper::filter(MSSpectrum& spectrum, int scan_number, bool decharge)
  {
    if (spectrum.empty()) return;
    const double carrier = positive_ ? Constants::PROTON_MASS_U : -Constants::PROTON_MASS_U;
    int assigned_max_charge = max_abs_charge_;
    double assigned_max_mass = std::numeric_limits<double>::infinity();
    if (!spectrum.getPrecursors().empty())
    {
      const auto& precursor = spectrum.getPrecursors().front();
      const int z = std::abs(precursor.getCharge());
      if (z > 0)
      {
        assigned_max_charge = std::min(assigned_max_charge, z);
        assigned_max_mass = (precursor.getMZ() - carrier) * z;
      }
    }
    // Keep configured detector bounds unchanged. Within that range, higher-charge
    // envelopes must remain available to reject lower-charge harmonic aliases
    // even when the precursor permits only the lower assignment charge/mass.
    // These competitors can veto assignments, but are never emitted above the
    // user's or precursor's assignment bounds.
    detector_.performSpectrumDeconvolution(spectrum, scan_number, PeakGroup());
    std::vector<bool> consumed(spectrum.size(), false);
    std::vector<std::set<int>> mono_charges(spectrum.size());
    Size missing_mono = 0, accepted_envelopes = 0;
    // Match the original samples recruited by FLASH. For profile data, an isotope
    // consists of several samples; retain every sample assigned to isotope zero.
    const auto index_of = [&spectrum](double mz)
    {
      const auto it = std::lower_bound(spectrum.begin(), spectrum.end(), mz,
        [](const Peak1D& peak, double value) { return peak.getMZ() < value; });
      return (it != spectrum.end() && it->getMZ() == mz) ? Size(it - spectrum.begin()) : spectrum.size();
    };
    for (const auto& group : detector_.getDeconvolvedSpectrum())
    {
      if (group.getMonoMass() > assigned_max_mass) continue;
      std::map<int, std::vector<Size>> by_charge;
      std::map<int, std::vector<Size>> monoisotopes;
      for (const auto& peak : group)
      {
        const int z = peak.abs_charge;
        if (z < 1 || z > assigned_max_charge || peak.isotopeIndex < 0 ||
            group.getChargeIsotopeCosine(z) < min_cos_ || group.getChargeSNR(z) <= min_snr_) continue;
        const Size index = index_of(peak.mz);
        if (index == spectrum.size()) continue;
        by_charge[z].push_back(index);
        if (peak.isotopeIndex == 0) monoisotopes[z].push_back(index);
      }
      for (const auto& [z, indices] : by_charge)
      {
        if (!monoisotopes.contains(z))
        {
          ++missing_mono;
          continue;
        }
        ++accepted_envelopes;
        for (Size index : monoisotopes.at(z)) mono_charges[index].insert(z);
        for (Size index : indices) consumed[index] = true;
      }
    }
    std::vector<Size> retained;
    MSSpectrum::IntegerDataArray charges;
    charges.setName("flash_charge");
    Size ambiguous = 0, retained_mono_samples = 0;
    for (Size i = 0; i < spectrum.size(); ++i)
    {
      if (consumed[i] && mono_charges[i].empty()) continue;
      retained.push_back(i);
      const auto& zs = mono_charges[i];
      const int z = zs.size() == 1 ? *zs.begin() : 0;
      if (!zs.empty()) ++retained_mono_samples;
      if (zs.size() > 1) ++ambiguous;
      charges.push_back(positive_ ? z : -z);
      if (decharge && z > 0)
      {
        spectrum[i].setMZ((spectrum[i].getMZ() - carrier) * z + carrier);
      }
    }
    const Size removed = spectrum.size() - retained.size();
    spectrum.select(retained); // Preserve every existing per-peak data array.
    auto& arrays = spectrum.getIntegerDataArrays();
    arrays.erase(std::remove_if(arrays.begin(), arrays.end(), [](const auto& a)
      { return a.getName() == "flash_charge"; }), arrays.end());
    arrays.push_back(std::move(charges));
    spectrum.setMetaValue("flash_detected_peak_groups", detector_.getDeconvolvedSpectrum().size());
    spectrum.setMetaValue("flash_isotope_peaks_removed", removed);
    spectrum.setMetaValue("flash_accepted_envelopes", accepted_envelopes);
    spectrum.setMetaValue("flash_missing_mono_envelopes", missing_mono);
    spectrum.setMetaValue("flash_ambiguous_mono_peaks", ambiguous);
    spectrum.setMetaValue("flash_retained_mono_samples", retained_mono_samples);
    spectrum.sortByPosition();
  }
}
