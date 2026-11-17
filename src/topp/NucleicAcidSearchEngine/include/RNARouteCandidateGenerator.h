// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#pragma once

#include "RNAModificationDiagnostics.h"
#include "RNASpectrumGraph.h"

#include <OpenMS/CHEMISTRY/NASequence.h>
#include <array>
#include <atomic>
#include <memory>

namespace OpenMS
{
/**
  @brief Generate precursor-compatible modified RNAs from anchored mass-graph routes.

  Database prefix/suffix tries prune graph traversal using canonical residue
  identities while preserving modification labels and isobars. One- and two-
  nucleotide anchors and one-residue gaps are supported. Compatible routes fix
  placements before mass-bounded completion of uncovered positions. Diagnostic-
  supported labile modifications may also occupy an otherwise canonical route
  position; the final scorer independently checks the actual modified theory.

  Resource exhaustion is reported explicitly. This search never silently drops
  candidates or switches an unrestricted modification panel to exhaustive
  expansion. The caller must fail the run when Result::complete is false.
*/
class RNARouteCandidateGenerator
{
public:
  struct Options
  {
    double tolerance = 20;
    bool tolerance_ppm = true;
    double precursor_tolerance = 20;
    bool precursor_tolerance_ppm = true;
    Size max_modifications = 25;
    /// Variable modifications without an intact modified graph placement.
    Size max_unanchored_modifications = 25;
    Size max_length = 25;
    Size max_peaks = 200;
    Size min_route_nodes = 2;
    Size max_states = 200000;
    Size max_completion_nodes = 5000000;
    Size max_candidates = 200000;
    Size max_supports = 400000;
    bool complete_unanchored_modifications = false;
    bool compress_equivalent_routes = true;
    bool index_precursor_mass = true;
    bool allow_full_length_routes = true;
    bool compress_route_states = true;
  };
  struct Scaffold
  {
    NASequence sequence;
    std::string label;
  };
  struct Candidate
  {
    Size database_index = 0;
    Size precursor_index = 0;
    NASequence sequence;
    /// Route sequences before diagnostic-supported labile completion, shared across placements.
    std::vector<std::shared_ptr<const Scaffold>> support_scaffolds;
  };
  struct Result
  {
    bool complete = true;
    std::string reason;
    Size graph_nodes = 0, graph_edges = 0, route_states = 0, routes = 0, completion_nodes = 0, support_variants = 0;
    std::vector<Candidate> candidates;
  };

  /**
    @brief Index canonical database prefixes/suffixes and modification choices.
    @param[in] database Digested database; order defines Candidate::database_index.
    @param[in] modifications Internal variable modifications.
    @param[in] options Tolerances, unrestricted-budget upper bound and explicit resource limits.
    @param[in] terminal_states Optional bit mask per database entry of permitted terminal assignments; bit
               @p 3*five+three with five in {OH, P} and three in {OH, P, cyclic-P}. Empty permits all six
               that agree with the termini stored in @p database; a mask replaces that check.
    @throws Exception::IllegalArgument for invalid options, terminal or linkage ('*') modifications,
            unsupported database termini, or a mask whose size differs from the database.
  */
  RNARouteCandidateGenerator(const std::vector<NASequence>& database,
                             const std::set<const Ribonucleotide*>& modifications,
                             const Options& options,
                             const std::vector<Byte>& terminal_states = {});

  /**
    @brief Search one spectrum with the same evidence/rules for targets and decoys.
    @param[in] spectrum Legacy processed MS2, retaining original peak identities.
    @param[in] precursor_masses Corrected neutral precursor hypotheses, in caller-defined order.
    @param[in] evidence Per-spectrum diagnostic permissions (all codes when filtering is disabled).
    @param[in] negative_mode Whether experimental ions are negative.
    @param[in] max_charge Largest charge magnitude; no FLASH MS2 processing is performed.
    @param[in] cancelled Optional cooperative cancellation after another spectrum made the file incomplete.
    @return Focused placements and auditable traversal counts, or explicit resource exhaustion.
  */
  Result generate(const MSSpectrum& spectrum,
                  const std::vector<double>& precursor_masses,
                  const RNAModificationDiagnostics::Evidence& evidence,
                  bool negative_mode,
                  Size max_charge,
                  const std::atomic<bool>* cancelled = nullptr) const;

private:
  struct TrieNode
  {
    std::array<Int, 4> next {-1, -1, -1, -1};
    std::vector<Size> descendants;
  };
  std::vector<NASequence> database_;
  std::set<const Ribonucleotide*> modifications_;
  std::array<std::vector<TrieNode>, 2> tries_;
  std::vector<double> minimum_masses_;
  std::vector<Byte> terminal_states_;
  Options options_;
};
} // namespace OpenMS
