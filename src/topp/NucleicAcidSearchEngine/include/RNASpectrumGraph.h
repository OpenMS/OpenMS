// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <array>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace OpenMS
{
/**
  @brief Compact, bounded RNA residue-mass DAG with tags queried on demand.

  The immutable RNA alphabet is prepared once and shared across spectra.
  Edges retain all matching residue alternatives through ranges in that table.
  RNATagCandidateFilter confirms series, charge and terminal origin before
  using a graph path to reject a candidate. Scoring spectra remain read-only.
*/
class RNASpectrumGraph
{
public:
  struct Options
  {
    double tolerance = 10.0;
    bool tolerance_ppm = true;
    Size min_tag_length = 2; ///< Minimum observed steps (edges); each spans up to max_step residues.
    Size max_tag_length = 4; ///< Maximum residues spanned by a tag.
    Size max_step = 2;
    Size max_peaks = 100;
    Size max_charge = 3;
    Size max_edges = 20000;
    Size max_query_edges = 20000;
    Size max_query_states = 100000;
  };

  struct Transition
  {
    double mass = 0.0;
    std::array<UInt, 2> residues {};
    UInt length = 0;
  };

  struct Alphabet
  {
    std::map<std::string, UInt> code_ids;
    std::vector<Transition> transitions;
  };

  /// Compact references preserve isobaric alternatives without copying label lists per path.
  struct Edge
  {
    UInt to = 0;
    UInt first_transition = 0;
    UInt past_transition = 0;
    UInt length = 0;
    bool operator==(const Edge&) const = default;
  };

  struct Layer
  {
    Size charge = 1;
    std::vector<std::vector<Edge>> outgoing;
  };

  struct Result
  {
    std::shared_ptr<const Alphabet> alphabet;
    Options options;
    std::vector<Size> peak_indices;
    std::vector<Layer> layers;
    std::string fallback_reason;
    Size peaks = 0;
    Size vertices = 0;
    Size edges = 0;
  };

  /**
    @brief Prepare the allowed RNA residue increments once for a batch.
    @param[in] alphabet Internal residues, including configured modifications.
    @param[in] options Fragment tolerances and graph/query bounds.
    @throws Exception::IllegalArgument If a tolerance, tag length or bound is invalid.
  */
  RNASpectrumGraph(const std::set<const Ribonucleotide*>& alphabet, const Options& options);

  /**
    @brief Build a mass DAG from a read-only processed MS2 spectrum.
    @param[in] spectrum Processed MS2; original indices are retained.
    @param[in] max_charge Largest fragment charge magnitude for this spectrum.
    @return Compact graph or an explicit enumeration fallback.
    @throws Exception::IllegalArgument If max_charge is zero.
  */
  Result build(const MSSpectrum& spectrum, Size max_charge) const;

  /**
    @brief Convenience entry point for a single spectrum.
    @param[in] spectrum Processed MS2.
    @param[in] alphabet Allowed RNA residues.
    @param[in] options Tolerances and resource bounds.
    @param[in] max_charge Largest fragment charge magnitude.
    @return Compact graph or enumeration fallback.
  */
  static Result build(const MSSpectrum& spectrum, const std::set<const Ribonucleotide*>& alphabet, const Options& options, Size max_charge);

private:
  std::shared_ptr<const Alphabet> alphabet_;
  Options options_;
  std::string fallback_reason_;
};
} // namespace OpenMS
