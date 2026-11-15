// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#pragma once

#include "RNASpectrumGraph.h"

#include <OpenMS/CHEMISTRY/NASequence.h>

namespace OpenMS
{
/**
  @brief Validate RNA graph tags against legal candidates and NASE's annotated theoretical ions.

  The same rules apply to targets and decoys. a-B ions remain available to full
  scoring but cannot seed constant-residue-increment graph paths. Unsupported
  linkage chemistry, short candidates and incomplete graphs use enumeration.
*/
class RNATagCandidateFilter
{
public:
  struct Evidence
  {
    bool supported = false;
    std::string fallback_reason;
    std::string ion_series;
    Size length = 0;
    Size charge = 0;
    Size start_ion = 0;
    bool accepted() const
    {
      return supported || ! fallback_reason.empty();
    }
    bool operator==(const Evidence&) const = default;
  };

  /// Candidate and theoretical-ion information shared by all precursor matches.
  struct PreparedQuery
  {
    struct Ion
    {
      Size position = 0;
      double mz = 0.0;
    };
    struct Group
    {
      char series = 'a';
      Size charge = 1;
      std::vector<Ion> ions;
      bool sorted_finite_mz = false;
    };
    std::shared_ptr<const RNASpectrumGraph::Alphabet> alphabet;
    Size sequence_length = 0;
    std::vector<UInt> codes;
    std::vector<Size> theoretical_charges;
    std::vector<Group> groups;
    bool supported_series = false;
    bool supported_residues = false;
    bool valid_metadata = false;
  };

  /**
    @brief Index the existing theoretical annotations once for a candidate/charge hypothesis.
    @param[in] candidate Fully placed RNA candidate.
    @param[in] alphabet Shared graph alphabet.
    @param[in] theoretical Existing NASE theoretical spectrum after blacklisting.
    @param[in] ions Enabled theoretical ion series.
    @return Immutable information reusable across experimental spectra.
  */
  static PreparedQuery prepare(const NASequence& candidate,
                               std::shared_ptr<const RNASpectrumGraph::Alphabet> alphabet,
                               const MSSpectrum& theoretical,
                               const std::set<std::string>& ions);

  /**
    @brief Match prepared theoretical annotations using sparse observed-ion state.
    @param[in] query Candidate/theory index from prepare().
    @param[in] graph Experimental mass DAG using the same alphabet.
    @param[in] experimental Original processed scoring spectrum.
    @param[in] tolerance Fragment tolerance.
    @param[in] tolerance_ppm Whether tolerance is in ppm (otherwise Da).
    @return The same evidence and bounded fallback as the reference matcher.
  */
  static Evidence
  match(const PreparedQuery& query, const RNASpectrumGraph::Result& graph, const MSSpectrum& experimental, double tolerance, bool tolerance_ppm);

  /**
    @brief Test residue compatibility before theoretical-spectrum generation.
    @param[in] candidate Fully placed RNA candidate; modification labels are preserved.
    @param[in] graph Compact experimental mass DAG.
    @param[in] ions Enabled theoretical ion series.
    @param[in] max_charge Maximum fragment charge magnitude (0 = no additional restriction).
    @return Whether at least one tag is compatible or enumeration fallback is needed.
  */
  static bool
  isCompatible(const NASequence& candidate, const RNASpectrumGraph::Result& graph, const std::set<std::string>& ions, Size max_charge = 0);

  /**
    @brief Conservatively query tags while internal modification placements remain unresolved.
    @param[in] alternatives Allowed residue labels at each sequence position.
    @param[in] graph Compact experimental graph.
    @param[in] ions Enabled ion series.
    @param[in] max_charge Maximum fragment charge magnitude (0 = unrestricted).
    @return True if a completion may have a tag or needs enumeration fallback.
  */
  static bool isCompatible(const std::vector<std::vector<const Ribonucleotide*>>& alternatives,
                           const RNASpectrumGraph::Result& graph,
                           const std::set<std::string>& ions,
                           Size max_charge = 0);

  /**
    @brief Require every observed path peak to match one annotated series, charge and orientation.
    @param[in] candidate Fully placed RNA candidate.
    @param[in] graph Compact experimental mass DAG.
    @param[in] experimental Original processed scoring spectrum (indices referenced by graph).
    @param[in] theoretical Existing NASE theoretical spectrum, with ion names and charges.
    @param[in] ions Enabled theoretical ion series.
    @param[in] tolerance Fragment tolerance.
    @param[in] tolerance_ppm Whether tolerance is in ppm (otherwise Da).
    @return Confirmed ladder evidence, rejection, or an explicit enumeration fallback.
  */
  static Evidence match(const NASequence& candidate,
                        const RNASpectrumGraph::Result& graph,
                        const MSSpectrum& experimental,
                        const MSSpectrum& theoretical,
                        const std::set<std::string>& ions,
                        double tolerance,
                        bool tolerance_ppm);

  /**
    @brief Exhaustive reference implementation for checking prepared matching and its witnesses.
    @copydetails match(const NASequence&, const RNASpectrumGraph::Result&, const MSSpectrum&, const MSSpectrum&, const std::set<std::string>&, double,
    bool)
  */
  static Evidence matchReference(const NASequence& candidate,
                                 const RNASpectrumGraph::Result& graph,
                                 const MSSpectrum& experimental,
                                 const MSSpectrum& theoretical,
                                 const std::set<std::string>& ions,
                                 double tolerance,
                                 bool tolerance_ppm);
};
} // namespace OpenMS
