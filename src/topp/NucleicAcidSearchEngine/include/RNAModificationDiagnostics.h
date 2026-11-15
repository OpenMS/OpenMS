// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace OpenMS
{
/**
  @brief Spectrum-specific modification evidence from raw, unfiltered MS2 peaks.

  Explicit diagnostic m/z values take precedence over theoretical single-residue
  backbone ions. An explicit rule marks a labile residue; its intact backbone
  mass need not be observed. Diagnostic support permits a chemistry but is not
  evidence for its position or a score bonus. Targets and decoys share a result.
*/
class RNAModificationDiagnostics
{
public:
  struct Evidence
  {
    std::set<std::string> tested;
    std::set<std::string> allowed;
    std::set<std::string> labile;
    std::map<std::string, double> matched_mz;
    std::map<std::string, Int> matched_charge;
    std::map<std::string, double> isotope_mz;
    bool accepts(const NASequence& candidate) const;
  };

  /**
    @brief Prepare theoretical one-nucleotide and explicit diagnostic rules.
    @param[in] modifications Enabled internal variable modifications.
    @param[in] explicit_mz Fixed diagnostic m/z values, keyed by residue code.
    @param[in] negative_mode Whether raw spectra contain negative ions.
    @param[in] tolerance Fragment tolerance.
    @param[in] tolerance_ppm Whether tolerance is in ppm (otherwise Da).
    @param[in] min_relative_intensity Minimum fraction of raw MS2 base-peak intensity.
    @param[in] max_charge Largest absolute charge for generic single-nucleotide ions. Multiply charged evidence requires an isotope partner.
  */
  RNAModificationDiagnostics(const std::set<const Ribonucleotide*>& modifications,
                             const std::map<std::string, std::vector<double>>& explicit_mz,
                             bool negative_mode,
                             double tolerance,
                             bool tolerance_ppm,
                             double min_relative_intensity = 0.0,
                             Int max_charge = 1);

  /// @brief Evaluate a raw spectrum without changing peaks, intensities or annotations.
  /// @param[in] raw Raw MS2, before deisotoping, normalization and noise filters.
  /// @return Supported modification codes and their diagnostic matches.
  Evidence evaluate(const MSSpectrum& raw) const;

private:
  struct Rule
  {
    struct Ion
    {
      double mz;
      Int charge;
    };
    std::string code;
    std::vector<Ion> ions;
    bool labile = false;
  };
  std::vector<Rule> rules_;
  double tolerance_;
  bool tolerance_ppm_;
  double min_relative_intensity_;
};
} // namespace OpenMS
