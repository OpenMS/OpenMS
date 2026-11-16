// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenDIATraces.h>

#include <array>
#include <string>

namespace OpenMS::OpenDIA
{
  /**
    @brief MS1 isotope evidence of one candidate window (PLAN_MS1), label-free and identical for targets and decoys.

    Inputs: the window's M+0..M+3 MS1 traces and its summed fragment trace on the 32-bin grid over apex +- 15 s, and the
    precursor's theoretical isotope envelope. Features, in names() order:
      m1_corr_xic   Pearson(M+0, fragment sum) over bins 6..25
      m1_iso_cos    cosine(central areas M+0..M+3 over bins 12..19, theoretical envelope)
      m1_iso_corr   mean of Pearson(M+0, M+1) and Pearson(M+0, M+2) over bins 6..25
      m1_apex_dt    |argmax M+0 - argmax fragment sum| over bins 6..25 (first maximum)
      m1_log_m0     log1p(central M+0 area, bins 12..19)
    Without MS1 signal in bins 6..25, features 1-3 are 0 and feature 4 is 20. A Pearson of a constant trace is 0.
  */
  struct MS1Evidence
  {
    static constexpr std::size_t FEATURES = 5;
    static constexpr std::size_t BINS = TraceEvidence::RT_BINS;
    static constexpr std::size_t TRACE_BEGIN = 6, TRACE_END = 26, CORE_BEGIN = 12, CORE_END = 20;

    static const std::array<std::string, FEATURES>& names();
    static std::array<double, FEATURES> features(const IsotopeTraces& isotopes, const std::array<double, BINS>& fragment_sum,
                                                 const std::array<double, 4>& envelope);

    /**
      @brief Theoretical M+0..M+3 abundances (summing to 1) of a peptide precursor from its modified sequence.

      Composition: residues plus H2O; UniMod:4 adds C2H3NO, UniMod:35 adds O. Heavy isotopes 13C 0.0107, 2H 0.000115,
      15N 0.00364, 17O 0.00038 and 18O 0.00205, 33S 0.0075 and 34S 0.0425, convolved per atom and truncated at M+3.
      @return All zeros for an unknown residue or modification.
    */
    static std::array<double, 4> envelope(const std::string& modified_sequence);
  };
}
