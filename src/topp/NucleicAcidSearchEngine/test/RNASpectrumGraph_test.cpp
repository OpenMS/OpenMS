// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#include "ModifiedNASequenceGenerator.h"
#include "RNASpectrumGraph.h"
#include "RNATagCandidateFilter.h"

#include <OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <algorithm>
#include <limits>
#include <random>

using namespace OpenMS;

namespace
{
MSSpectrum theoretical(const NASequence& sequence, const std::string& ion, Int charge = -1)
{
  NucleicAcidSpectrumGenerator generator;
  auto parameters = generator.getParameters();
  for (const auto& code : {"a", "b", "c", "d", "w", "x", "y", "z", "a-B"})
  {
    parameters.setValue("add_" + std::string(code) + "_ions", code == ion ? "true" : "false");
  }
  parameters.setValue("add_first_prefix_ion", "true");
  parameters.setValue("add_metainfo", "true");
  parameters.setValue("add_precursor_peaks", "false");
  generator.setParameters(parameters);
  MSSpectrum spectrum;
  generator.getSpectrum(spectrum, sequence, charge, charge);
  return spectrum;
}

MSSpectrum ladder(const MSSpectrum& spectrum, const std::set<std::string>& names)
{
  MSSpectrum result = spectrum;
  std::vector<Size> indices;
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    if (names.contains(spectrum.getStringDataArrays()[0][i])) { indices.push_back(i); }
  }
  result.select(indices);
  return result;
}

std::set<const Ribonucleotide*> alphabet()
{
  std::set<const Ribonucleotide*> result;
  for (const auto& code : {"A", "C", "G", "U", "mA", "Am"})
  {
    result.insert(RibonucleotideDB::getInstance()->getRibonucleotide(code));
  }
  return result;
}
RNATagCandidateFilter::Evidence checkedMatch(const NASequence& candidate,
                                             const RNASpectrumGraph::Result& graph,
                                             const MSSpectrum& observed,
                                             const MSSpectrum& model,
                                             const std::set<std::string>& ions,
                                             double tolerance,
                                             bool ppm)
{
  const auto expected = RNATagCandidateFilter::matchReference(candidate, graph, observed, model, ions, tolerance, ppm);
  const auto prepared = RNATagCandidateFilter::prepare(candidate, graph.alphabet, model, ions);
  const auto actual = RNATagCandidateFilter::match(prepared, graph, observed, tolerance, ppm);
  if (actual != expected || RNATagCandidateFilter::match(candidate, graph, observed, model, ions, tolerance, ppm) != expected)
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Prepared graph evidence differs from reference");
  }
  return actual;
}
} // namespace

START_TEST(RNASpectrumGraph, "$Id$")

START_SECTION((prefix and suffix tags preserve charge, orientation, and scoring spectra))
{
  const auto sequence = NASequence::fromString("ACGU");
  for (const auto& ion : {"a", "b", "c", "d", "w", "x", "y", "z"})
  {
    for (Int charge : {-2, -1, 1, 2})
    {
      const auto model = theoretical(sequence, ion, charge);
      const auto observed = model;
      const auto graph = RNASpectrumGraph::build(observed, alphabet(), {}, std::abs(charge));
      TEST_TRUE(graph.fallback_reason.empty())
      TEST_FALSE(graph.layers.empty())
      const auto evidence = checkedMatch(sequence, graph, observed, model, {ion}, 10, true);
      TEST_TRUE(evidence.supported)
      TEST_EQUAL(evidence.ion_series, ion)
      TEST_EQUAL(evidence.charge, std::abs(charge))
      TEST_EQUAL(observed, model)
      TEST_FALSE(RNATagCandidateFilter::isCompatible(NASequence::fromString("AAAA"), graph, {ion}))
    }
  }
}
END_SECTION

START_SECTION((mass differences alone cannot establish ion origin or terminal chemistry))
{
  const auto sequence = NASequence::fromString("ACGU");
  const auto model = theoretical(sequence, "b");
  auto shifted = model;
  for (auto& peak : shifted)
  {
    peak.setMZ(peak.getMZ() + 100.0);
  }
  const auto graph = RNASpectrumGraph::build(shifted, alphabet(), {}, 1);
  TEST_TRUE(RNATagCandidateFilter::isCompatible(sequence, graph, {"b"}))
  TEST_FALSE(checkedMatch(sequence, graph, shifted, model, {"b"}, 10, true).accepted())
  TEST_FALSE(checkedMatch(sequence, graph, model, theoretical(sequence, "y"), {"b"}, 10, true).accepted())
  TEST_FALSE(
    checkedMatch(NASequence::fromString("pACGU"), graph, model, theoretical(NASequence::fromString("pACGU"), "b"), {"b"}, 10, true).accepted())
}
END_SECTION

START_SECTION((one - and two - residue anchors can start after missing initial ions))
{
  const auto sequence = NASequence::fromString("ACGUA");
  const auto model = theoretical(sequence, "b");
  RNASpectrumGraph::Options options;
  options.min_tag_length = 1;
  auto observed = ladder(model, {"b2", "b3"});
  auto graph = RNASpectrumGraph::build(observed, alphabet(), options, 1);
  auto evidence = checkedMatch(sequence, graph, observed, model, {"b"}, 10, true);
  TEST_TRUE(evidence.supported)
  TEST_EQUAL(evidence.length, 1)
  TEST_EQUAL(evidence.start_ion, 2)
  observed = ladder(model, {"b2", "b4"});
  options.min_tag_length = 2;
  graph = RNASpectrumGraph::build(observed, alphabet(), options, 1);
  evidence = checkedMatch(sequence, graph, observed, model, {"b"}, 10, true);
  TEST_TRUE(evidence.supported)
  TEST_EQUAL(evidence.length, 2)
  TEST_EQUAL(evidence.start_ion, 2)
  options.max_step = 1;
  graph = RNASpectrumGraph::build(observed, alphabet(), options, 1);
  TEST_TRUE(graph.layers.empty())
  TEST_EQUAL(graph.fallback_reason, "no_tags")
}
END_SECTION

START_SECTION((endpoint tolerances and mixed one / two - residue steps preserve valid paths))
{
  const auto sequence = NASequence::fromString("ACGUA");
  const auto model = theoretical(sequence, "b");
  auto observed = ladder(model, {"b1", "b2", "b4"});
  for (Size i = 0; i < observed.size(); ++i)
  {
    observed[i].setMZ(observed[i].getMZ() * (1.0 + (i % 2 ? -9.9 : 9.9) * 1e-6));
  }
  RNASpectrumGraph::Options options;
  options.min_tag_length = 3;
  options.max_tag_length = 3;
  const RNASpectrumGraph builder(alphabet(), options);
  const auto graph = builder.build(observed, 1);
  TEST_TRUE(graph.fallback_reason.empty())
  TEST_TRUE(RNATagCandidateFilter::isCompatible(sequence, graph, {"b"}))
  const auto evidence = checkedMatch(sequence, graph, observed, model, {"b"}, 10, true);
  TEST_TRUE(evidence.supported)
  TEST_EQUAL(evidence.length, 3)
  TEST_EQUAL(evidence.start_ion, 1)
}
END_SECTION

START_SECTION((isobaric modifications remain alternatives and modification generation is reused))
{
  const auto sequence = NASequence::fromString("C[mA]GU");
  const auto model = theoretical(sequence, "b");
  const auto graph = RNASpectrumGraph::build(model, alphabet(), {}, 1);
  std::vector<NASequence> candidates;
  ModifiedNASequenceGenerator::applyVariableModifications(
    {RibonucleotideDB::getInstance()->getRibonucleotide("mA"), RibonucleotideDB::getInstance()->getRibonucleotide("Am")},
    NASequence::fromString("CAGU"), 1, candidates, true);
  Size supported = 0;
  for (const auto& candidate : candidates)
  {
    const auto evidence = checkedMatch(candidate, graph, model, theoretical(candidate, "b"), {"b"}, 10, true);
    if (evidence.supported) { ++supported; }
  }
  TEST_EQUAL(supported, 2)
}
END_SECTION

START_SECTION((graph and candidate - query bounds retain enumeration fallback))
{
  const auto sequence = NASequence::fromString("ACGUA");
  const auto model = theoretical(sequence, "b");
  RNASpectrumGraph::Options options;
  options.max_edges = 1;
  auto graph = RNASpectrumGraph::build(model, alphabet(), options, 1);
  TEST_EQUAL(graph.fallback_reason, "edge_limit")
  TEST_TRUE(graph.layers.empty())
  TEST_TRUE(checkedMatch(sequence, graph, model, model, {"b"}, 10, true).accepted())
  options = {};
  options.max_query_edges = 1;
  graph = RNASpectrumGraph::build(model, alphabet(), options, 1);
  TEST_TRUE(RNATagCandidateFilter::isCompatible(sequence, graph, {"b"}))
  auto evidence = checkedMatch(sequence, graph, model, model, {"b"}, 10, true);
  TEST_EQUAL(evidence.fallback_reason, "query_limit")
  TEST_TRUE(evidence.accepted())
  options = {};
  options.max_query_states = 1;
  graph = RNASpectrumGraph::build(model, alphabet(), options, 1);
  TEST_TRUE(RNATagCandidateFilter::isCompatible(sequence, graph, {"b"}))
  evidence = checkedMatch(sequence, graph, model, model, {"b"}, 10, true);
  TEST_EQUAL(evidence.fallback_reason, "query_state_limit")
  TEST_TRUE(evidence.accepted())
}
END_SECTION

START_SECTION((sparse spectra, unsupported chemistry, and a - B - only scoring retain enumeration))
{
  const auto sequence = NASequence::fromString("ACGU");
  const auto model = theoretical(sequence, "b");
  const auto graph = RNASpectrumGraph::build(model, alphabet(), {}, 1);
  TEST_EQUAL(checkedMatch(sequence, graph, model, model, {"a-B"}, 10, true).fallback_reason, "no_supported_ion_series")
  TEST_EQUAL(checkedMatch(NASequence::fromString("ACG"), graph, model, model, {"b"}, 10, true).fallback_reason, "short_candidate")
  TEST_EQUAL(checkedMatch(NASequence::fromString("A[Cm]GU"), graph, model, model, {"b"}, 10, true).fallback_reason, "unsupported_residue")
  auto sparse = ladder(model, {"b1"});
  const auto empty_graph = RNASpectrumGraph::build(sparse, alphabet(), {}, 1);
  TEST_EQUAL(empty_graph.fallback_reason, "no_tags")
  TEST_TRUE(checkedMatch(sequence, empty_graph, sparse, model, {"b"}, 10, true).accepted())
  auto no_annotations = model;
  no_annotations.getStringDataArrays().clear();
  TEST_EQUAL(checkedMatch(sequence, graph, model, no_annotations, {"b"}, 10, true).fallback_reason, "missing_theoretical_annotations")
  // A single-charge theoretical model cannot use an unrelated graph layer.
  const auto other_charge = theoretical(sequence, "b", -2);
  const auto evidence = checkedMatch(sequence, graph, model, other_charge, {"b"}, 10, true);
  TEST_EQUAL(evidence.fallback_reason, "no_tags_for_charge")
  TEST_TRUE(evidence.accepted())
}
END_SECTION

START_SECTION((input validation and traversal are deterministic))
{
  const auto sequence = NASequence::fromString("ACGU");
  const auto model = theoretical(sequence, "b");
  const auto first = RNASpectrumGraph::build(model, alphabet(), {}, 1);
  const auto second = RNASpectrumGraph::build(model, alphabet(), {}, 1);
  TEST_EQUAL(first.layers.size(), second.layers.size())
  TEST_EQUAL(first.edges, second.edges)
  TEST_EQUAL(first.vertices, second.vertices)
  TEST_TRUE(first.peak_indices == second.peak_indices)
  for (Size i = 0; i < first.layers.size(); ++i)
  {
    TEST_EQUAL(first.layers[i].charge, second.layers[i].charge)
    TEST_TRUE(first.layers[i].outgoing == second.layers[i].outgoing)
  }
  RNASpectrumGraph::Options options;
  options.tolerance = std::numeric_limits<double>::quiet_NaN();
  TEST_EXCEPTION(Exception::IllegalArgument, RNASpectrumGraph::build(model, alphabet(), options, 1))
  options = {};
  options.min_tag_length = 5;
  TEST_EXCEPTION(Exception::IllegalArgument, RNASpectrumGraph::build(model, alphabet(), options, 1))
}
END_SECTION


START_SECTION((partial residue domains retain every precursor and tag compatible completion))
{
  const auto source = NASequence::fromString("CAGUA");
  const auto desired = NASequence::fromString("C[mA]GUA");
  const std::set<const Ribonucleotide*> mods {RibonucleotideDB::getInstance()->getRibonucleotide("mA"),
                                              RibonucleotideDB::getInstance()->getRibonucleotide("Am"),
                                              RibonucleotideDB::getInstance()->getRibonucleotide("3'-p")};
  std::vector<NASequence> reference;
  ModifiedNASequenceGenerator::applyVariableModifications(mods, source, 2, reference);
  for (const auto& ion : {"b", "w"})
  {
    const auto model = theoretical(desired, ion, -2);
    const auto graph = RNASpectrumGraph::build(model, alphabet(), {}, 2);
    const double target = desired.getMonoWeight();
    ModifiedNASequenceGenerator::PlacementConstraints constraints;
    constraints.accepts = [&](double low, double high, const auto& alternatives) {
      return low <= target && target <= high && RNATagCandidateFilter::isCompatible(alternatives, graph, {ion}, 2);
    };
    std::vector<NASequence> actual, expected;
    const auto statistics = ModifiedNASequenceGenerator::applyVariableModificationsConstrained(mods, source, 2, actual, constraints);
    for (const auto& candidate : reference)
    {
      if (std::abs(candidate.getMonoWeight() - target) < 1e-7 && RNATagCandidateFilter::isCompatible(candidate, graph, {ion}, 2))
      {
        expected.push_back(candidate);
      }
    }
    std::sort(actual.begin(), actual.end());
    std::sort(expected.begin(), expected.end());
    TEST_TRUE(actual == expected)
    TEST_TRUE(expected.size() >= 2)
    TEST_TRUE(statistics.pruned_branches > 0)
    TEST_FALSE(statistics.enumeration_fallback)
  }
  ModifiedNASequenceGenerator::ResidueAlternatives alternatives;
  for (Size i = 0; i < source.size(); ++i)
  {
    alternatives.push_back({source[i]});
  }
  const auto model = theoretical(desired, "b");
  const auto graph = RNASpectrumGraph::build(model, alphabet(), {}, 1);
  alternatives[0].push_back(RibonucleotideDB::getInstance()->getRibonucleotide("Cm"));
  TEST_TRUE(RNATagCandidateFilter::isCompatible(alternatives, graph, {"b"}))
  alternatives[0].pop_back();
  alternatives[1].push_back(RibonucleotideDB::getInstance()->getRibonucleotide("mA"));
  TEST_TRUE(RNATagCandidateFilter::isCompatible(alternatives, graph, {"b"}))
  RNASpectrumGraph::Options bounded;
  bounded.max_query_edges = 1;
  const auto bounded_graph = RNASpectrumGraph::build(model, alphabet(), bounded, 1);
  TEST_TRUE(RNATagCandidateFilter::isCompatible(alternatives, bounded_graph, {"b"}))
}
END_SECTION


START_SECTION((partial domains preserve singleton query limit fallback))
{
  const auto* db = RibonucleotideDB::getInstance();
  auto alphabet = std::make_shared<RNASpectrumGraph::Alphabet>();
  alphabet->code_ids = {{"Am", 0}, {"mA", 1}, {"G", 2}};
  alphabet->transitions = {{343, {0, 0}, 1}, {343, {1, 0}, 1}, {345, {2, 0}, 1}};
  RNASpectrumGraph::Result graph;
  graph.alphabet = alphabet;
  graph.options.max_step = 1;
  graph.options.max_query_edges = 8;
  graph.peak_indices = {0, 1, 2, 3, 4};
  RNASpectrumGraph::Layer layer;
  layer.outgoing = {{{1, 0, 2, 1}}, {}, {{3, 0, 2, 1}}, {{4, 2, 3, 1}}, {}};
  graph.layers.push_back(layer);
  const auto candidate = NASequence::fromString("[mA][mA][mA][mA]");
  ModifiedNASequenceGenerator::ResidueAlternatives alternatives;
  for (Size i = 0; i < candidate.size(); ++i)
  {
    alternatives.push_back({candidate[i]});
  }
  alternatives[1].push_back(db->getRibonucleotide("Am"));
  // No completed compatible tag exists. The singleton hits its work limit;
  // a broader domain must therefore fail open too, despite earlier isobar matches.
  TEST_TRUE(RNATagCandidateFilter::isCompatible(candidate, graph, {"b"}))
  TEST_TRUE(RNATagCandidateFilter::isCompatible(alternatives, graph, {"b"}))
  graph.options.max_query_edges = 100;
  TEST_FALSE(RNATagCandidateFilter::isCompatible(candidate, graph, {"b"}))
  TEST_FALSE(RNATagCandidateFilter::isCompatible(alternatives, graph, {"b"}))
}
END_SECTION


START_SECTION((prepared sparse matching preserves witnesses, query limits and tolerance endpoints))
{
  const auto sequence = NASequence::fromString("ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU");
  for (const auto& ion : {"b", "w"})
  {
    const auto model = theoretical(sequence, ion, -2);
    for (Size maximum : {1, 2, 4, 32})
    {
      for (Size budget : {1, 8, 100, 20000})
      {
        RNASpectrumGraph::Options options;
        options.min_tag_length = 1;
        options.max_tag_length = maximum;
        options.max_query_edges = budget;
        const auto graph = RNASpectrumGraph::build(model, alphabet(), options, 2);
        checkedMatch(sequence, graph, model, model, {ion}, 10, true);
      }
    }
    for (bool ppm : {false, true})
    {
      RNASpectrumGraph::Options options;
      options.tolerance_ppm = ppm;
      options.tolerance = ppm ? 10.0 : 0.01;
      auto observed = model;
      for (Size i = 0; i < observed.size(); ++i)
      {
        const double delta = options.tolerance * (ppm ? model[i].getMZ() * 1e-6 : 1.0);
        observed[i].setMZ(model[i].getMZ() + (i % 2 ? -delta : delta));
      }
      const auto graph = RNASpectrumGraph::build(observed, alphabet(), options, 2);
      checkedMatch(sequence, graph, observed, model, {ion}, options.tolerance, ppm);
      auto reversed = model;
      std::reverse(reversed.begin(), reversed.end());
      for (auto& array : reversed.getStringDataArrays())
      {
        std::reverse(array.begin(), array.end());
      }
      for (auto& array : reversed.getIntegerDataArrays())
      {
        std::reverse(array.begin(), array.end());
      }
      checkedMatch(sequence, graph, observed, reversed, {ion}, options.tolerance, ppm);
    }
  }
}
END_SECTION

START_SECTION((prepared matching agrees with the reference on deterministic noisy and permuted candidates))
{
  std::mt19937 engine(91237);
  const std::string bases = "ACGU";
  for (Size iteration = 0; iteration < 200; ++iteration)
  {
    std::string text;
    for (Size i = 0; i < 8; ++i)
    {
      text += bases[engine() % bases.size()];
    }
    const auto source = NASequence::fromString(text);
    const std::string ion = iteration % 2 ? "b" : "w";
    const Int charge = iteration % 3 ? -1 : 2;
    const auto model = theoretical(source, ion, charge);
    auto observed = model;
    observed.getStringDataArrays().clear();
    observed.getIntegerDataArrays().clear();
    for (auto& peak : observed)
    {
      peak.setMZ(peak.getMZ() * (1.0 + (static_cast<Int>(engine() % 61) - 30) * 1e-6));
    }
    for (Size i = 0; i < 10; ++i)
    {
      observed.push_back(Peak1D(100.0 + engine() % 200000 / 100.0, 1.0));
    }
    observed.sortByPosition();
    RNASpectrumGraph::Options options;
    options.max_query_edges = iteration % 4 ? 20000 : 8;
    options.max_query_states = iteration % 7 ? 100000 : 1;
    const auto graph = RNASpectrumGraph::build(observed, alphabet(), options, std::abs(charge));
    checkedMatch(source, graph, observed, model, {ion}, 10, true);
    std::reverse(text.begin(), text.end());
    const auto candidate = NASequence::fromString(text);
    checkedMatch(candidate, graph, observed, theoretical(candidate, ion, charge), {ion}, 10, true);
  }
}
END_SECTION

END_TEST
