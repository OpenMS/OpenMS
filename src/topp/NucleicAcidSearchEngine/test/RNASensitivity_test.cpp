// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#include "ModifiedNASequenceGenerator.h"
#include "RNAModificationDiagnostics.h"
#include "RNARouteCandidateGenerator.h"
#include "RNARouteScorer.h"

#include <OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>

using namespace OpenMS;

namespace
{
MSSpectrum theory(const NASequence& sequence, Int charge = -1, const std::set<std::string>& ions = {"c", "y"})
{
  NucleicAcidSpectrumGenerator generator;
  auto parameters = generator.getParameters();
  for (const auto& code : {"a", "b", "c", "d", "w", "x", "y", "z", "a-B"})
  {
    parameters.setValue("add_" + std::string(code) + "_ions", ions.contains(code) ? "true" : "false");
  }
  parameters.setValue("add_first_prefix_ion", "true");
  parameters.setValue("add_metainfo", "true");
  parameters.setValue("add_precursor_peaks", "false");
  generator.setParameters(parameters);
  MSSpectrum model;
  generator.getSpectrum(model, sequence, charge < 0 ? -1 : 1, charge);
  return model;
}
} // namespace

START_TEST(RNASensitivity, "$Id$")
auto* db = RibonucleotideDB::getInstance();

START_SECTION((premodified database residues cannot expand the variable modification panel))
{
  const auto input = NASequence::fromString("[mA]UUUU");
  const auto* methyl = db->getRibonucleotide("mA");
  RNARouteCandidateGenerator::Options options;
  options.complete_unanchored_modifications = true;
  options.max_modifications = options.max_unanchored_modifications = 1;
  for (const auto& [base, placed] : std::vector<std::pair<std::string, std::string>> {{"GGCCA", "GGCC[mA]"}, {"AACGU", "A[mA]CGU"}})
  {
    const auto unmodified = NASequence::fromString(base), modified = NASequence::fromString(placed);
    const auto observed = theory(modified, -1, {"c"});
    for (Int policy = 0; policy < 3; ++policy)
    {
      const std::set<const Ribonucleotide*> variables = policy ? std::set<const Ribonucleotide*> {methyl} : std::set<const Ribonucleotide*> {};
      RNAModificationDiagnostics::Evidence evidence;
      if (policy == 2) { evidence.tested.insert("mA"); }
      const auto result
        = RNARouteCandidateGenerator({unmodified, input}, variables, options).generate(observed, {modified.getMonoWeight()}, evidence, true, 1);
      TEST_TRUE(result.complete)
      const bool found = std::any_of(result.candidates.begin(), result.candidates.end(), [&](const auto& hit) { return hit.sequence == modified; });
      TEST_EQUAL(found, policy == 1)
    }
  }
  // The original premodified input remains searchable without declaring it variable.
  const auto own = RNARouteCandidateGenerator({input}, {}, options).generate(theory(input, -1, {"c"}), {input.getMonoWeight()}, {}, true, 1);
  TEST_TRUE(own.complete)
  TEST_TRUE(std::any_of(own.candidates.begin(), own.candidates.end(), [&](const auto& hit) { return hit.sequence == input; }))
}
END_SECTION

START_SECTION((Da route anchors retain the observed tolerance after decharging))
{
  const auto candidate = NASequence::fromString("ACUGUACGp");
  RNARouteCandidateGenerator::Options options;
  options.tolerance = 0.02;
  options.tolerance_ppm = false;
  options.max_modifications = 0;
  options.min_route_nodes = 3;
  for (Int polarity : {-1, 1})
  {
    const auto model = theory(candidate, polarity * 2, {"c"});
    for (double error : {0.015, 0.025})
    {
      MSSpectrum observed;
      auto& charges = observed.getIntegerDataArrays().emplace_back();
      charges.setName("NASE:legacy_fragment_charge");
      for (Size i = 0; i < model.size(); ++i)
      {
        if (std::abs(model.getIntegerDataArrays()[0][i]) != 2) { continue; }
        observed.emplace_back(model[i].getMZ() + (observed.size() % 2 ? error : -error), 10);
        charges.push_back(2);
      }
      observed.sortByPosition();
      const auto result = RNARouteCandidateGenerator({candidate}, {}, options).generate(observed, {candidate.getMonoWeight()}, {}, polarity < 0, 2);
      TEST_TRUE(result.complete)
      const bool found = std::any_of(result.candidates.begin(), result.candidates.end(), [&](const auto& hit) { return hit.sequence == candidate; });
      TEST_EQUAL(found, error < options.tolerance)
    }
  }
}
END_SECTION

START_SECTION((raw diagnostic evidence is spectrum - specific and does not assert localization))
{
  const std::set<const Ribonucleotide*> mods {db->getRibonucleotide("mU"), db->getRibonucleotide("t6A")};
  RNAModificationDiagnostics diagnostics(mods, {{"t6A", {279.08418}}}, true, 20, true);
  MSSpectrum raw;
  raw.push_back(Peak1D(279.084625, 100));
  raw.push_back(Peak1D(319.033719, 50));
  const auto original = raw;
  const auto evidence = diagnostics.evaluate(raw);
  TEST_TRUE(evidence.allowed.contains("mU"))
  TEST_TRUE(evidence.allowed.contains("t6A"))
  TEST_TRUE(evidence.labile.contains("t6A"))
  TEST_TRUE(evidence.accepts(NASequence::fromString("A[t6A]C[mU]G")))
  TEST_EQUAL(raw, original)
  TEST_TRUE(diagnostics.evaluate(MSSpectrum()).allowed.empty())
  TEST_FALSE(diagnostics.evaluate(MSSpectrum()).accepts(NASequence::fromString("A[t6A]CG")))
  // A fixed or input modification outside this variable panel is not rejected.
  TEST_TRUE(diagnostics.evaluate(MSSpectrum()).accepts(NASequence::fromString("A[mC]G")))
  RNAModificationDiagnostics threshold(mods, {{"t6A", {279.08418}}}, true, 20, true, 0.75);
  TEST_FALSE(threshold.evaluate(raw).allowed.contains("mU"))
  TEST_TRUE(threshold.evaluate(raw).allowed.contains("t6A"))
}
END_SECTION

START_SECTION((multiply charged diagnostics require a raw isotope partner and respect precursor charge))
{
  const auto* residue = db->getRibonucleotide("mU");
  const std::set<const Ribonucleotide*> mods {residue};
  const NASequence single({residue}, nullptr, nullptr);
  for (Int polarity : {-1, 1})
  {
    RNAModificationDiagnostics legacy(mods, {}, polarity < 0, 20, true);
    RNAModificationDiagnostics charged(mods, {}, polarity < 0, 20, true, 0.0, 2);
    const double mz = single.getMonoWeight(NASequence::CIon, polarity * 2) / 2;
    MSSpectrum raw;
    raw.push_back(Peak1D(mz, 100));
    TEST_FALSE(charged.evaluate(raw).allowed.contains("mU"))
    RNAModificationDiagnostics wide(mods, {}, polarity < 0, 0.6, false, 0.0, 2);
    TEST_FALSE(wide.evaluate(raw).allowed.contains("mU"))
    raw.push_back(Peak1D(mz + Constants::C13C12_MASSDIFF_U / 2, 50));
    TEST_TRUE(wide.evaluate(raw).allowed.contains("mU"))
    TEST_FALSE(legacy.evaluate(raw).allowed.contains("mU"))
    const auto evidence = charged.evaluate(raw);
    TEST_TRUE(evidence.allowed.contains("mU"))
    TEST_EQUAL(evidence.matched_charge.at("mU"), 2)
    TEST_REAL_SIMILAR(evidence.matched_mz.at("mU"), mz)
    TEST_REAL_SIMILAR(evidence.isotope_mz.at("mU"), raw[1].getMZ())
    Precursor precursor;
    precursor.setCharge(polarity);
    raw.getPrecursors().push_back(precursor);
    TEST_FALSE(charged.evaluate(raw).allowed.contains("mU"))
    raw.getPrecursors()[0].setCharge(polarity * 2);
    TEST_TRUE(charged.evaluate(raw).allowed.contains("mU"))
  }
  // Fixed labile rules specify an observed m/z, not a neutral mass to expand.
  RNAModificationDiagnostics fixed(mods, {{"mU", {300.0}}}, true, 20, true, 0.0, 2);
  MSSpectrum half;
  half.push_back(Peak1D(150, 100));
  half.push_back(Peak1D(150 + Constants::C13C12_MASSDIFF_U / 2, 50));
  TEST_FALSE(fixed.evaluate(half).allowed.contains("mU"))
  half.push_back(Peak1D(300, 100));
  TEST_TRUE(fixed.evaluate(half).labile.contains("mU"))
}
END_SECTION

START_SECTION((route scoring requires anchored coherent ladders and unique observed nodes))
{
  const auto candidate = NASequence::fromString("ACUGUACGp");
  const auto model = theory(candidate, -2);
  auto observed = model;
  const auto query = RNATagCandidateFilter::prepare(candidate, nullptr, model, {"c", "y"});
  const auto result = RNARouteScorer::score(query, observed, true, {});
  TEST_TRUE(result.score > 0)
  TEST_TRUE(result.prefix_length >= 5)
  TEST_TRUE(result.suffix_length >= 5)
  TEST_EQUAL(result.annotations.size(), result.nodes)
  TEST_TRUE(result.score <= 2 * result.nodes)
  std::set<double> matched_mz;
  for (const auto& annotation : result.annotations)
  {
    matched_mz.insert(annotation.mz);
  }
  TEST_EQUAL(matched_mz.size(), result.nodes)
  TEST_REAL_SIMILAR(result.matched_intensity, double(result.nodes))
  TEST_EQUAL(model, observed)
  // Interior evidence without a one-/two-nucleotide origin must not create a route.
  std::vector<Size> interior;
  for (Size i = 0; i < model.size(); ++i)
  {
    const auto& name = model.getStringDataArrays()[0][i];
    if (name != "c1" && name != "c2" && name != "y1" && name != "y2") { interior.push_back(i); }
  }
  observed.select(interior);
  TEST_EQUAL(RNARouteScorer::score(query, observed, true, {}).score, 0)
  // Extra local background reduces normalized intensity support.
  observed = model;
  observed.getStringDataArrays().clear();
  observed.getIntegerDataArrays().clear();
  observed.push_back(Peak1D(model[0].getMZ() + 0.1, 100));
  observed.sortByPosition();
  const auto noisy = RNARouteScorer::score(query, observed, true, {});
  TEST_TRUE(noisy.score < result.score)
}
END_SECTION

START_SECTION((route scoring retains lighter prefixes needed by a later charge hypothesis))
{
  for (bool negative_mode : {true, false})
  {
    RNARouteScorer::Options options;
    options.min_nodes = 3;
    options.max_unknown_charge = 2;
    const auto sequence = NASequence::fromString("AAAAA");
    const auto model = theory(sequence, negative_mode ? -2 : 2, {"c"});
    const auto query = RNATagCandidateFilter::prepare(sequence, nullptr, model, {"c"});
    TEST_TRUE(query.valid_metadata)
    double first_mz = 0, third_mz = 0;
    for (const auto& group : query.groups)
    {
      if (group.charge != 1) { continue; }
      for (const auto& ion : group.ions)
      {
        if (ion.position == 2) { first_mz = ion.mz; }
        if (ion.position == 3) { third_mz = ion.mz; }
      }
    }
    TEST_TRUE(first_mz > 0)
    TEST_TRUE(third_mz > first_mz)
    const double delta = first_mz * options.tolerance * 1e-6;
    MSSpectrum observed;
    observed.emplace_back(first_mz, 10);
    observed.emplace_back(first_mz + 0.001, 2 / std::exp(-2 * std::pow(0.001 / delta, 2)));
    observed.emplace_back(third_mz, 3);
    // The poly-A c2/charge-1 and c4/charge-2 ions share peak 0. Only
    // peak 0 has ambiguous charge; the nearby alternative is singly charged.
    auto& charges = observed.getIntegerDataArrays().emplace_back();
    charges.setName("NASE:legacy_fragment_charge");
    for (Int charge : {0, 1, 1})
    {
      charges.push_back(charge);
    }
    const auto prepared = RNARouteScorer::prepareSpectrum(observed, negative_mode, options);
    const auto matched = RNARouteScorer::match(query, observed, prepared, negative_mode, options);
    TEST_EQUAL(matched.routes.size(), 1)
    ABORT_IF(matched.routes.empty())
    const auto& route = matched.routes[0];
    TEST_EQUAL(route.matches.size(), 3)
    ABORT_IF(route.matches.size() != 3)
    TEST_EQUAL(route.matches[0].peak, 1)
    TEST_EQUAL(route.matches[1].peak, 2)
    TEST_EQUAL(route.matches[2].peak, 0)
    TEST_EQUAL(route.matches[2].charge, 2)
    TEST_EQUAL(route.length, 4)
    TEST_REAL_SIMILAR(route.weighted_intensity, 15)
    const auto scored = RNARouteScorer::score(query, observed, prepared, negative_mode, options);
    TEST_EQUAL(scored.nodes, 3)
    TEST_EQUAL(scored.annotations.size(), 3)
  }
}
END_SECTION

START_SECTION((route scoring retains a shorter anchored prefix when a used peak recurs))
{
  RNARouteScorer::Options options;
  options.min_nodes = 4;
  options.max_step = 1;
  options.max_unknown_charge = 4;
  MSSpectrum observed;
  for (double mass : {300.0, 600.0, 900.0, 1500.0})
  {
    observed.emplace_back(mass - Constants::PROTON_MASS_U, 1);
  }
  RNATagCandidateFilter::PreparedQuery query;
  query.valid_metadata = true;
  query.sequence_length = 6;
  query.groups = {{'c', 1, {{1, observed[0].getMZ()}, {2, observed[1].getMZ()}, {3, observed[2].getMZ()}, {5, observed[3].getMZ()}}, true},
                  {'c', 4, {{4, observed[0].getMZ()}}, true}};
  const auto prepared = RNARouteScorer::prepareSpectrum(observed, true, options);
  const auto matched = RNARouteScorer::match(query, observed, prepared, true, options);
  TEST_EQUAL(matched.routes.size(), 1)
  ABORT_IF(matched.routes.empty())
  const auto& route = matched.routes[0];
  TEST_EQUAL(route.matches.size(), 4)
  ABORT_IF(route.matches.size() != 4)
  for (Size i = 0; i < route.matches.size(); ++i)
  {
    TEST_EQUAL(route.matches[i].position, i + 2)
  }
  TEST_EQUAL(route.matches[0].peak, 1)
  TEST_EQUAL(route.matches[1].peak, 2)
  TEST_EQUAL(route.matches[2].peak, 0)
  TEST_EQUAL(route.matches[3].peak, 3)
}
END_SECTION

START_SECTION((route scoring agrees with exhaustive paths for small charge-ambiguous spectra))
{
  for (unsigned mask = 0; mask < 16; ++mask)
  {
    for (const std::vector<double>& intensities : {std::vector<double> {1, 1, 1, 1}, {1, 2, 3, 4}, {10, 2, 3, 8}})
    {
      MSSpectrum observed;
      for (Size peak = 0; peak < 4; ++peak)
      {
        if (mask & (1u << peak)) { observed.emplace_back(300.0 * (peak + 1) - Constants::PROTON_MASS_U, intensities[peak]); }
      }
      RNATagCandidateFilter::PreparedQuery query;
      query.valid_metadata = true;
      query.sequence_length = 7;
      std::vector<RNARouteScorer::MatchedNode> nodes;
      for (Size charge = 1; charge <= 3; ++charge)
      {
        RNATagCandidateFilter::PreparedQuery::Group group;
        group.series = 'c';
        group.charge = charge;
        for (Size position = 1; position <= 6; ++position)
        {
          group.ions.push_back({position, 300.0 * position / charge - Constants::PROTON_MASS_U});
          Size index = 0;
          for (Size peak = 0; peak < 4; ++peak)
          {
            if (! (mask & (1u << peak))) { continue; }
            if ((peak + 1) * charge == position) { nodes.push_back({index, position, charge, intensities[peak], 300.0 * position, {}}); }
            ++index;
          }
        }
        query.groups.push_back(std::move(group));
      }
      std::sort(nodes.begin(), nodes.end(),
                [](const auto& a, const auto& b) { return std::tie(a.position, a.peak, a.charge) < std::tie(b.position, b.peak, b.charge); });
      for (Size step : {Size(1), Size(2)})
      {
        for (Size minimum = 1; minimum <= 6; ++minimum)
        {
          std::vector<Size> path, expected;
          double expected_weight = 0;
          std::set<Size> used;
          // Enumerate every complete path with the full peak history. This
          // oracle does not share the scorer's DP or expired-identity pruning.
          const auto visit = [&](const auto& self, double weight) -> void {
            if (path.size() >= minimum
                && (path.size() > expected.size()
                    || (path.size() == expected.size()
                        && (weight > expected_weight
                            || (weight == expected_weight
                                && std::lexicographical_compare(path.rbegin(), path.rend(), expected.rbegin(), expected.rend()))))))
            {
              expected = path;
              expected_weight = weight;
            }
            for (Size next = 0; next < nodes.size(); ++next)
            {
              const auto& node = nodes[next];
              if (used.contains(node.peak)) { continue; }
              if (path.empty())
              {
                if (node.position > 2) { continue; }
              }
              else
              {
                const auto& previous = nodes[path.back()];
                if (node.position <= previous.position || node.position - previous.position > step || node.neutral <= previous.neutral) { continue; }
              }
              path.push_back(next);
              used.insert(node.peak);
              self(self, weight + node.weighted);
              used.erase(node.peak);
              path.pop_back();
            }
          };
          visit(visit, 0);
          RNARouteScorer::Options options;
          options.tolerance = 1e-6;
          options.tolerance_ppm = false;
          options.max_unknown_charge = 3;
          options.max_step = step;
          options.min_nodes = minimum;
          const auto prepared = RNARouteScorer::prepareSpectrum(observed, true, options);
          const auto actual = RNARouteScorer::match(query, observed, prepared, true, options);
          TEST_EQUAL(actual.routes.size(), expected.empty() ? 0 : 1)
          if (expected.empty() || actual.routes.empty()) { continue; }
          const auto& route = actual.routes[0];
          TEST_EQUAL(route.matches.size(), expected.size())
          TEST_REAL_SIMILAR(route.weighted_intensity, expected_weight)
          if (route.matches.size() != expected.size()) { continue; }
          for (Size i = 0; i < expected.size(); ++i)
          {
            TEST_EQUAL(route.matches[i].peak, nodes[expected[i]].peak)
            TEST_EQUAL(route.matches[i].position, nodes[expected[i]].position)
            TEST_EQUAL(route.matches[i].charge, nodes[expected[i]].charge)
          }
        }
      }
    }
  }
}
END_SECTION

START_SECTION((route scoring expires irrelevant peak identities and fails on state exhaustion))
{
  MSSpectrum observed;
  RNATagCandidateFilter::PreparedQuery query;
  query.valid_metadata = true;
  query.sequence_length = 13;
  RNATagCandidateFilter::PreparedQuery::Group group;
  group.series = 'c';
  group.charge = 1;
  for (Size position = 1; position <= 12; ++position)
  {
    const double mz = 300.0 * position - Constants::PROTON_MASS_U;
    observed.emplace_back(mz, 2);
    observed.emplace_back(mz + 0.001, 1);
    group.ions.push_back({position, mz});
  }
  query.groups.push_back(std::move(group));
  RNARouteScorer::Options options;
  options.tolerance = 0.01;
  options.tolerance_ppm = false;
  options.min_nodes = 12;
  options.max_step = 1;
  options.max_unknown_charge = 1;
  // There are 2^12 valid histories, but expired peak identities permit one
  // state per endpoint. The exact search fits into 24 retained states.
  options.max_states = 24;
  const auto prepared = RNARouteScorer::prepareSpectrum(observed, true, options);
  const auto matched = RNARouteScorer::match(query, observed, prepared, true, options);
  TEST_EQUAL(matched.routes.size(), 1)
  ABORT_IF(matched.routes.empty())
  TEST_EQUAL(matched.routes[0].matches.size(), 12)
  TEST_REAL_SIMILAR(matched.routes[0].weighted_intensity, 24)
  options.max_states = 1;
  TEST_EXCEPTION(Exception::InvalidValue, RNARouteScorer::match(query, observed, prepared, true, options))
  options.max_states = 0;
  TEST_EXCEPTION(Exception::IllegalArgument, RNARouteScorer::match(query, observed, prepared, true, options))
}
END_SECTION

START_SECTION((mass accuracy calibration preserves matching and the default score))
{
  const auto candidate = NASequence::fromString("ACUGUACGp");
  const auto model = theory(candidate);
  auto observed = model;
  for (auto& peak : observed)
  {
    peak.setMZ(peak.getMZ() * (1 + 10e-6));
  }
  const auto query = RNATagCandidateFilter::prepare(candidate, nullptr, model, {"c", "y"});
  RNARouteScorer::Options options;
  const auto inherited = RNARouteScorer::score(query, observed, true, options);
  options.mass_accuracy_tolerance = options.tolerance;
  const auto explicit_default = RNARouteScorer::score(query, observed, true, options);
  TEST_EQUAL(inherited.score, explicit_default.score)
  options.mass_accuracy_tolerance = 5;
  const auto calibrated = RNARouteScorer::score(query, observed, true, options);
  TEST_EQUAL(calibrated.nodes, inherited.nodes)
  TEST_EQUAL(calibrated.prefix_length, inherited.prefix_length)
  TEST_EQUAL(calibrated.suffix_length, inherited.suffix_length)
  TEST_TRUE(calibrated.score < inherited.score)
  TEST_REAL_SIMILAR(calibrated.score, calibrated.nodes * (1 + calibrated.matched_intensity / calibrated.local_intensity))
  options.mass_accuracy_tolerance = std::numeric_limits<double>::denorm_min();
  TEST_TRUE(std::isfinite(RNARouteScorer::score(query, model, true, options).score))
  options.mass_accuracy_tolerance = -1;
  TEST_EXCEPTION(Exception::IllegalArgument, RNARouteScorer::score(query, observed, true, options))
}
END_SECTION

START_SECTION((base loss ladders score separately from intact a ions))
{
  const auto sequence = NASequence::fromString("ACUACUAGp");
  const auto observed = theory(sequence, -1, {"a-B"});
  const auto default_query = RNATagCandidateFilter::prepare(sequence, nullptr, observed, {"a-B"});
  TEST_EQUAL(RNARouteScorer::score(default_query, observed, true, {}).score, 0)
  const auto enabled_query = RNATagCandidateFilter::prepare(sequence, nullptr, observed, {"a-B"}, true);
  const auto scored = RNARouteScorer::score(enabled_query, observed, true, {});
  TEST_TRUE(scored.score > 0)
  TEST_TRUE(scored.nodes >= sequence.size() - 2)
  for (const auto& annotation : scored.annotations)
  {
    TEST_TRUE(annotation.annotation.ends_with("-B"))
  }
  TEST_TRUE(std::all_of(enabled_query.groups.begin(), enabled_query.groups.end(), [](const auto& group) { return group.series == 'A'; }))
}
END_SECTION

START_SECTION((the unanchored budget distinguishes intact and diagnostic only labile placements))
{
  const auto base = NASequence::fromString("ACAUAG");
  const auto target = NASequence::fromString("[t6A]C[t6A]U[t6A]G");
  const auto* modification = db->getRibonucleotide("t6A");
  RNAModificationDiagnostics::Evidence evidence;
  evidence.tested = evidence.allowed = evidence.labile = {"t6A"};
  RNARouteCandidateGenerator::Options options;
  options.max_modifications = 6;
  options.max_unanchored_modifications = 2;
  const auto limited
    = RNARouteCandidateGenerator({base}, {modification}, options).generate(theory(base), {target.getMonoWeight()}, evidence, true, 1);
  TEST_TRUE(limited.complete)
  TEST_TRUE(std::none_of(limited.candidates.begin(), limited.candidates.end(), [&](const auto& c) { return c.sequence == target; }))
  options.max_unanchored_modifications = 3;
  const auto complete
    = RNARouteCandidateGenerator({base}, {modification}, options).generate(theory(base), {target.getMonoWeight()}, evidence, true, 1);
  TEST_TRUE(complete.complete)
  TEST_TRUE(std::any_of(complete.candidates.begin(), complete.candidates.end(), [&](const auto& c) { return c.sequence == target; }))
  options.max_unanchored_modifications = 0;
  const auto intact
    = RNARouteCandidateGenerator({base}, {modification}, options).generate(theory(target), {target.getMonoWeight()}, evidence, true, 1);
  TEST_TRUE(intact.complete)
  TEST_TRUE(std::any_of(intact.candidates.begin(), intact.candidates.end(), [&](const auto& c) { return c.sequence == target; }))
  options.max_modifications = 2;
  const auto globally_limited
    = RNARouteCandidateGenerator({base}, {modification}, options).generate(theory(target), {target.getMonoWeight()}, evidence, true, 1);
  TEST_TRUE(globally_limited.complete)
  TEST_TRUE(std::none_of(globally_limited.candidates.begin(), globally_limited.candidates.end(), [&](const auto& c) { return c.sequence == target; }))
}
END_SECTION

START_SECTION((cancelled route searches fail explicitly and do not export partial candidates))
{
  const auto sequence = NASequence::fromString("ACUACUAG");
  RNARouteCandidateGenerator generator({sequence}, {}, {});
  std::atomic<bool> cancelled {true};
  const auto result = generator.generate(theory(sequence), {sequence.getMonoWeight()}, {}, true, 1, &cancelled);
  TEST_FALSE(result.complete)
  TEST_STRING_EQUAL(result.reason, "cancelled_after_incomplete_search")
  TEST_TRUE(result.candidates.empty())
  cancelled.store(false);
  const auto complete = generator.generate(theory(sequence), {sequence.getMonoWeight()}, {}, true, 1, &cancelled);
  TEST_TRUE(complete.complete)
  TEST_FALSE(complete.candidates.empty())
}
END_SECTION

START_SECTION((bounded generic completion indexes preserve every candidate and scaffold witness))
{
  const auto base = NASequence::fromString("ACUACUAG");
  const std::set<const Ribonucleotide*> mods {db->getRibonucleotide("mA"), db->getRibonucleotide("mU"), db->getRibonucleotide("t6A")};
  RNAModificationDiagnostics::Evidence evidence;
  evidence.tested = evidence.allowed = {"mA", "mU", "t6A"};
  evidence.labile = {"t6A"};
  std::vector<NASequence> exhaustive;
  ModifiedNASequenceGenerator::applyVariableModifications(mods, base, 2, exhaustive, true);
  RNARouteCandidateGenerator::Options options;
  options.complete_unanchored_modifications = true;
  options.max_unanchored_modifications = 2;
  options.max_modifications = 8;
  const RNARouteCandidateGenerator fast({base}, mods, options);
  options.index_precursor_mass = false;
  options.compress_equivalent_routes = false;
  options.compress_route_states = false;
  const RNARouteCandidateGenerator oracle({base}, mods, options);
  const auto keys = [](const auto& candidates) {
    std::vector<std::tuple<Size, Size, std::string, std::vector<std::string>>> result;
    for (const auto& candidate : candidates)
    {
      std::vector<std::string> labels;
      for (const auto& scaffold : candidate.support_scaffolds)
      {
        labels.push_back(scaffold->label);
      }
      std::sort(labels.begin(), labels.end());
      result.emplace_back(candidate.database_index, candidate.precursor_index, candidate.sequence.toString(), std::move(labels));
    }
    std::sort(result.begin(), result.end());
    return result;
  };
  for (auto target : exhaustive)
  {
    auto scaffold = target;
    for (Size pos = 0; pos < scaffold.size(); ++pos)
    {
      if (scaffold[pos]->getCode() == "t6A") { scaffold.set(pos, db->getRibonucleotide("A")); }
    }
    auto observed = theory(scaffold);
    std::vector<Size> early;
    for (Size i = 0; i < observed.size(); ++i)
    {
      const auto name = observed.getStringDataArrays()[0][i];
      if (name == "c1" || name == "c2" || name == "y1" || name == "y2") { early.push_back(i); }
    }
    observed.select(early);
    const auto actual = fast.generate(observed, {target.getMonoWeight()}, evidence, true, 1);
    const auto expected = oracle.generate(observed, {target.getMonoWeight()}, evidence, true, 1);
    TEST_TRUE(actual.complete)
    TEST_TRUE(expected.complete)
    TEST_TRUE(keys(actual.candidates) == keys(expected.candidates))
  }
}
END_SECTION

START_SECTION((bounded completion recovers modifications beyond the observed route))
{
  const auto base = NASequence::fromString("ACUACUAG");
  const auto target = NASequence::fromString("ACU[mA]C[mU]AG");
  auto observed = theory(base, -1, {"c"});
  std::vector<Size> keep;
  for (Size i = 0; i < observed.size(); ++i)
  {
    if (observed.getStringDataArrays()[0][i] == "c1" || observed.getStringDataArrays()[0][i] == "c2") { keep.push_back(i); }
  }
  observed.select(keep);
  const std::set<const Ribonucleotide*> mods {db->getRibonucleotide("mA"), db->getRibonucleotide("mU")};
  RNARouteCandidateGenerator::Options options;
  options.max_modifications = 8;
  options.complete_unanchored_modifications = true;
  options.max_unanchored_modifications = 1;
  const auto limited = RNARouteCandidateGenerator({base}, mods, options).generate(observed, {target.getMonoWeight()}, {}, true, 1);
  TEST_TRUE(limited.complete)
  TEST_TRUE(std::none_of(limited.candidates.begin(), limited.candidates.end(), [&](const auto& c) { return c.sequence == target; }))
  options.max_unanchored_modifications = 2;
  const auto completed = RNARouteCandidateGenerator({base}, mods, options).generate(observed, {target.getMonoWeight()}, {}, true, 1);
  TEST_TRUE(completed.complete)
  TEST_TRUE(std::any_of(completed.candidates.begin(), completed.candidates.end(), [&](const auto& c) { return c.sequence == target; }))
  options.index_precursor_mass = false;
  options.compress_route_states = false;
  const auto oracle = RNARouteCandidateGenerator({base}, mods, options).generate(observed, {target.getMonoWeight()}, {}, true, 1);
  std::set<std::string> actual, expected;
  for (const auto& c : completed.candidates)
  {
    actual.insert(c.sequence.toString());
  }
  for (const auto& c : oracle.candidates)
  {
    expected.insert(c.sequence.toString());
  }
  TEST_TRUE(actual == expected)
}
END_SECTION

START_SECTION((graph routes retain multiple modifications and all six physical terminal assignments))
{
  const auto base = NASequence::fromString("ACUACUAG");
  const auto modified = NASequence::fromString("[mA][mC][mU][mA]C[mU]AG");
  const std::set<const Ribonucleotide*> mods {db->getRibonucleotide("mA"), db->getRibonucleotide("mC"), db->getRibonucleotide("mU")};
  RNARouteCandidateGenerator::Options options;
  options.max_modifications = 8;
  options.max_completion_nodes = 1000000;
  RNARouteCandidateGenerator generator({base, NASequence::fromString("UACACUAG")}, mods, options);
  for (Size five = 0; five < 2; ++five)
  {
    for (Size three = 0; three < 3; ++three)
    {
      auto target = modified;
      target.setFivePrimeMod(five ? db->getRibonucleotide("5'-p") : nullptr);
      target.setThreePrimeMod(three == 1 ? db->getRibonucleotide("3'-p") : (three == 2 ? db->getRibonucleotide("3'-c") : nullptr));
      const auto spectrum = theory(target);
      const auto generated = generator.generate(spectrum, {target.getMonoWeight()}, {}, true, 1);
      TEST_TRUE(generated.complete)
      TEST_TRUE(generated.routes > 0)
      const auto found = std::find_if(generated.candidates.begin(), generated.candidates.end(),
                                      [&](const auto& candidate) { return candidate.database_index == 0 && candidate.sequence == target; });
      TEST_TRUE(found != generated.candidates.end())
      for (const auto& placed : generated.candidates)
      {
        TEST_TRUE(std::abs(placed.sequence.getMonoWeight() - target.getMonoWeight()) <= target.getMonoWeight() * 20e-6)
      }
    }
  }
  options.max_states = 1;
  const auto limited = RNARouteCandidateGenerator({base}, mods, options).generate(theory(modified), {modified.getMonoWeight()}, {}, true, 1);
  TEST_FALSE(limited.complete)
  TEST_TRUE(limited.candidates.empty())
}
END_SECTION

START_SECTION((terminal-state masks restrict assignments per database entry and linkage modifications are rejected))
{
  const auto base = NASequence::fromString("ACGUAGCUAG");
  auto target = base;
  target.setThreePrimeMod(db->getRibonucleotide("3'-p"));
  auto isobar = base;
  isobar.setFivePrimeMod(db->getRibonucleotide("5'-p"));
  const auto spectrum = theory(target);
  const auto sequences = [](const auto& result) {
    std::set<std::string> labels;
    for (const auto& candidate : result.candidates)
    {
      labels.insert(candidate.sequence.toString());
    }
    return labels;
  };
  // Without a mask, both isobaric phosphate assignments are generated.
  const auto all = RNARouteCandidateGenerator({base}, {}, {}).generate(spectrum, {target.getMonoWeight()}, {}, true, 1);
  TEST_TRUE(all.complete)
  TEST_TRUE(sequences(all) == (std::set<std::string> {target.toString(), isobar.toString()}))
  // Bit 3*five+three: only OH/P (an interior RNase T1-like product).
  const auto masked = RNARouteCandidateGenerator({base}, {}, {}, {Byte(1 << 1)}).generate(spectrum, {target.getMonoWeight()}, {}, true, 1);
  TEST_TRUE(masked.complete)
  TEST_TRUE(sequences(masked) == (std::set<std::string> {target.toString()}))
  const auto excluded = RNARouteCandidateGenerator({base}, {}, {}, {Byte(1 << 0)}).generate(spectrum, {target.getMonoWeight()}, {}, true, 1);
  TEST_TRUE(excluded.complete)
  TEST_TRUE(excluded.candidates.empty())
  TEST_EXCEPTION(Exception::IllegalArgument, RNARouteCandidateGenerator({base}, {}, {}, {Byte(1), Byte(1)}))
  // Phosphorothioate linkages have neighbor-dependent masses outside the graph alphabet.
  TEST_EXCEPTION(Exception::IllegalArgument, RNARouteCandidateGenerator({base}, {db->getRibonucleotide("A*")}, {}))
  TEST_EXCEPTION(Exception::IllegalArgument, RNARouteCandidateGenerator({NASequence::fromString("AC[G*]UAG")}, {}, {}))
}
END_SECTION

START_SECTION((diagnostic permissions constrain graph labels and preserve isobaric identities))
{
  const auto base = NASequence::fromString("ACACGU");
  const auto target = NASequence::fromString("A[mC]ACGUp");
  const std::set<const Ribonucleotide*> mods {db->getRibonucleotide("mC"), db->getRibonucleotide("Cm")};
  RNARouteCandidateGenerator::Options options;
  options.max_modifications = 1;
  RNARouteCandidateGenerator generator({base}, mods, options);
  auto evidence = RNAModificationDiagnostics::Evidence {};
  evidence.tested = {"mC", "Cm"};
  auto generated = generator.generate(theory(target), {target.getMonoWeight()}, evidence, true, 1);
  TEST_TRUE(generated.complete)
  TEST_TRUE(generated.candidates.empty())
  evidence.allowed = evidence.tested;
  generated = generator.generate(theory(target), {target.getMonoWeight()}, evidence, true, 1);
  TEST_TRUE(generated.complete)
  std::set<std::string> sequences;
  for (const auto& candidate : generated.candidates)
  {
    sequences.insert(candidate.sequence.toString());
  }
  TEST_TRUE(sequences.contains("A[mC]ACGUp"))
  TEST_TRUE(sequences.contains("A[Cm]ACGUp"))
}
END_SECTION

START_SECTION((focused placement matches an exhaustive small - space oracle and compressed route masks))
{
  const auto base = NASequence::fromString("ACUACG");
  const std::set<const Ribonucleotide*> mods {db->getRibonucleotide("mA"), db->getRibonucleotide("mC"), db->getRibonucleotide("mU")};
  std::vector<NASequence> exhaustive;
  ModifiedNASequenceGenerator::applyVariableModifications(mods, base, 2, exhaustive, true);
  RNARouteCandidateGenerator::Options options;
  options.max_modifications = 2;
  const RNARouteCandidateGenerator fast({base}, mods, options);
  options.compress_equivalent_routes = false;
  options.compress_route_states = false;
  options.index_precursor_mass = false;
  const RNARouteCandidateGenerator reference({base}, mods, options);
  for (auto candidate : exhaustive)
  {
    candidate.setThreePrimeMod(db->getRibonucleotide("3'-p"));
    for (Int charge : {-2, 2})
    {
      const auto observed = theory(candidate, charge);
      const auto a = fast.generate(observed, {candidate.getMonoWeight()}, {}, charge < 0, 2);
      const auto b = reference.generate(observed, {candidate.getMonoWeight()}, {}, charge < 0, 2);
      TEST_TRUE(a.complete)
      TEST_TRUE(b.complete)
      std::set<std::string> left, right;
      for (const auto& placed : a.candidates)
      {
        left.insert(placed.sequence.toString());
      }
      for (const auto& placed : b.candidates)
      {
        right.insert(placed.sequence.toString());
      }
      TEST_TRUE(left == right)
      TEST_TRUE(left.contains(candidate.toString()))
    }
  }
}
END_SECTION

START_SECTION((labile completion does not require intact modified backbone fragments))
{
  const auto base = NASequence::fromString("ACUGUACG");
  const auto target = NASequence::fromString("ACUGU[t6A]CGp");
  auto observed = theory(target);
  std::vector<Size> kept;
  for (Size i = 0; i < observed.size(); ++i)
  {
    const auto& name = observed.getStringDataArrays()[0][i];
    const Size position = std::stoul(name.substr(1));
    if ((name[0] == 'c' && position <= 5) || (name[0] == 'y' && position <= 2)) { kept.push_back(i); }
  }
  observed.select(kept);
  RNARouteCandidateGenerator generator({base}, {db->getRibonucleotide("t6A")}, {});
  RNAModificationDiagnostics::Evidence evidence;
  evidence.tested = {"t6A"};
  evidence.allowed = {"t6A"};
  evidence.labile = {"t6A"};
  const auto result = generator.generate(observed, {target.getMonoWeight()}, evidence, true, 1);
  TEST_TRUE(result.complete)
  TEST_TRUE(std::any_of(result.candidates.begin(), result.candidates.end(), [&](const auto& candidate) { return candidate.sequence == target; }))
  const auto scored = RNARouteScorer::score(RNATagCandidateFilter::prepare(target, nullptr, theory(target), {"c", "y"}), observed, true, {});
  TEST_TRUE(scored.score > 0)
  TEST_EQUAL(scored.prefix_length, 5)
  TEST_EQUAL(scored.suffix_length, 2)
  evidence.allowed.clear();
  evidence.labile.clear();
  TEST_TRUE(generator.generate(observed, {target.getMonoWeight()}, evidence, true, 1).candidates.empty())
}
END_SECTION

START_SECTION((explicit termini and higher legacy - assigned charges survive focused generation))
{
  const auto target = NASequence::fromString("pACUGUACGp");
  auto observed = theory(target, -4);
  std::vector<Size> charge_four;
  for (Size i = 0; i < observed.size(); ++i)
  {
    if (std::abs(observed.getIntegerDataArrays()[0][i]) == 4) { charge_four.push_back(i); }
  }
  observed.select(charge_four);
  MSSpectrum::IntegerDataArray legacy;
  legacy.setName("NASE:legacy_fragment_charge");
  legacy.assign(observed.size(), 4);
  observed.getIntegerDataArrays().push_back(legacy);
  const auto result = RNARouteCandidateGenerator({target}, {}, {}).generate(observed, {target.getMonoWeight()}, {}, true, 3);
  TEST_TRUE(result.complete)
  TEST_FALSE(result.candidates.empty())
  for (const auto& candidate : result.candidates)
  {
    TEST_TRUE(candidate.sequence.hasFivePrimeMod())
    TEST_TRUE(candidate.sequence.hasThreePrimeMod())
    TEST_STRING_EQUAL(candidate.sequence.getThreePrimeMod()->getCode(), "3'-p")
  }
}
END_SECTION

START_SECTION((mass - equivalent terminal hypotheses retain an exact score tie))
{
  const auto three = NASequence::fromString("ACUGUACGp");
  const auto five = NASequence::fromString("pACUGUACG");
  const auto observed = theory(three, -2);
  const std::set<std::string> ions {"a", "b", "c", "d", "w", "x", "y", "z"};
  const auto score = [&](const NASequence& sequence) {
    const auto model = theory(sequence, -2, ions);
    return RNARouteScorer::score(RNATagCandidateFilter::prepare(sequence, nullptr, model, ions), observed, true, {}).score;
  };
  TEST_TRUE(score(three) > 0)
  TEST_EQUAL(score(three), score(five))
}
END_SECTION

START_SECTION((precursor indexing preserves candidates at ppm and Da boundaries with negative labile shifts))
{
  const auto base = NASequence::fromString("ACUACG");
  auto loss = *db->getRibonucleotide("A");
  loss.setCode("lossA");
  loss.setFormula(EmpiricalFormula("C9H11N5O4"));
  loss.setMonoMass(loss.getFormula().getMonoWeight());
  auto target = base;
  target.set(0, &loss);
  target.setThreePrimeMod(db->getRibonucleotide("3'-p"));
  RNAModificationDiagnostics::Evidence evidence;
  evidence.tested = {"lossA"};
  evidence.allowed = {"lossA"};
  evidence.labile = {"lossA"};
  for (bool ppm : {true, false})
  {
    RNARouteCandidateGenerator::Options options;
    options.precursor_tolerance_ppm = ppm;
    options.precursor_tolerance = ppm ? 20 : 0.02;
    const double mass = target.getMonoWeight();
    const double tolerance = options.precursor_tolerance * (ppm ? mass * 1e-6 : 1);
    const std::vector<double> precursors {mass - tolerance + 1e-9, mass, mass + tolerance - 1e-9};
    const auto fast = RNARouteCandidateGenerator({base}, {&loss}, options).generate(theory(target), precursors, evidence, true, 1);
    options.index_precursor_mass = false;
    const auto reference = RNARouteCandidateGenerator({base}, {&loss}, options).generate(theory(target), precursors, evidence, true, 1);
    TEST_TRUE(fast.complete)
    TEST_TRUE(reference.complete)
    std::set<std::pair<Size, std::string>> left, right;
    for (const auto& candidate : fast.candidates)
    {
      left.emplace(candidate.precursor_index, candidate.sequence.toString());
    }
    for (const auto& candidate : reference.candidates)
    {
      right.emplace(candidate.precursor_index, candidate.sequence.toString());
    }
    TEST_TRUE(left == right)
    for (Size precursor = 0; precursor < precursors.size(); ++precursor)
    {
      TEST_TRUE(left.contains({precursor, target.toString()}))
    }
    const auto outside = RNARouteCandidateGenerator({base}, {&loss}, options).generate(theory(target), {mass + tolerance + 1e-5}, evidence, true, 1);
    TEST_TRUE(std::none_of(outside.candidates.begin(), outside.candidates.end(), [&](const auto& candidate) { return candidate.sequence == target; }))
  }
}
END_SECTION

START_SECTION((indexed completion preserves isobaric placements and scaffolds at precursor boundaries))
{
  const auto base = NASequence::fromString("AAAAAACG");
  auto loss = *db->getRibonucleotide("A");
  loss.setCode("lossA");
  loss.setFormula(EmpiricalFormula("C9H11N5O4"));
  loss.setMonoMass(loss.getFormula().getMonoWeight());
  auto isobar = *db->getRibonucleotide("A");
  isobar.setCode("sameA");
  isobar.setMonoMass(isobar.getFormula().getMonoWeight());
  const auto* methyl = db->getRibonucleotide("m6A");
  auto target = base;
  target.set(2, &loss);
  target.set(3, methyl);
  target.setThreePrimeMod(db->getRibonucleotide("3'-p"));
  RNAModificationDiagnostics::Evidence evidence;
  evidence.tested = evidence.allowed = evidence.labile = {"lossA", "m6A", "sameA"};
  for (bool generic : {false, true})
  {
    for (bool ppm : {false, true})
    {
      for (Size budget : {0, 1, 2, 3})
      {
        RNARouteCandidateGenerator::Options options;
        options.complete_unanchored_modifications = generic;
        options.max_modifications = 2;
        options.max_unanchored_modifications = budget;
        options.precursor_tolerance_ppm = ppm;
        options.precursor_tolerance = ppm ? 20 : 0.0001;
        const double mass = target.getMonoWeight();
        const double tolerance = options.precursor_tolerance * (ppm ? mass * 1e-6 : 1);
        const std::vector<double> precursors {mass - tolerance + 1e-9, mass, mass + tolerance - 1e-9};
        const auto fast = RNARouteCandidateGenerator({base}, {&loss, methyl, &isobar}, options).generate(theory(base), precursors, evidence, true, 1);
        options.index_precursor_mass = false;
        const auto reference
          = RNARouteCandidateGenerator({base}, {&loss, methyl, &isobar}, options).generate(theory(base), precursors, evidence, true, 1);
        const auto keys = [](const auto& candidates) {
          std::set<std::tuple<Size, std::string, std::vector<std::string>>> result;
          for (const auto& candidate : candidates)
          {
            std::vector<std::string> labels;
            for (const auto& scaffold : candidate.support_scaffolds)
            {
              labels.push_back(scaffold->label);
            }
            std::sort(labels.begin(), labels.end());
            result.emplace(candidate.precursor_index, candidate.sequence.toString(), labels);
          }
          return result;
        };
        TEST_TRUE(fast.complete)
        TEST_TRUE(reference.complete)
        TEST_TRUE(keys(fast.candidates) == keys(reference.candidates))
        for (Size precursor = 0; budget >= 2 && precursor < precursors.size(); ++precursor)
        {
          TEST_TRUE(std::any_of(fast.candidates.begin(), fast.candidates.end(),
                                [&](const auto& candidate) { return candidate.precursor_index == precursor && candidate.sequence == target; }))
        }
      }
    }
  }
}
END_SECTION

START_SECTION((a full - length route can place a terminal modification without scoring the precursor))
{
  const auto base = NASequence::fromString("AACUCG");
  const auto target = NASequence::fromString("AACUC[m7G]p");
  auto observed = theory(target, -1, {"d"});
  observed.getStringDataArrays().clear();
  observed.getIntegerDataArrays().clear();
  observed.emplace_back(target.getMonoWeight(NASequence::Full, -1), 0.7);
  observed.sortByPosition();
  RNARouteCandidateGenerator::Options options;
  const std::set<const Ribonucleotide*> modifications {db->getRibonucleotide("m7G")};
  const auto result = RNARouteCandidateGenerator({base}, modifications, options).generate(observed, {target.getMonoWeight()}, {}, true, 1);
  TEST_TRUE(result.complete)
  TEST_TRUE(std::any_of(result.candidates.begin(), result.candidates.end(), [&](const auto& candidate) { return candidate.sequence == target; }))
  options.allow_full_length_routes = false;
  const auto restricted = RNARouteCandidateGenerator({base}, modifications, options).generate(observed, {target.getMonoWeight()}, {}, true, 1);
  TEST_TRUE(restricted.complete)
  TEST_TRUE(
    std::none_of(restricted.candidates.begin(), restricted.candidates.end(), [&](const auto& candidate) { return candidate.sequence == target; }))
  const auto model = theory(target, -1, {"d"});
  const auto score = RNARouteScorer::score(RNATagCandidateFilter::prepare(target, nullptr, model, {"d"}), observed, true, {});
  TEST_EQUAL(score.nodes, target.size() - 1)
  for (const auto& annotation : score.annotations)
  {
    TEST_FALSE(annotation.annotation == "d6")
  }
}
END_SECTION

START_SECTION((attainable suffix masses preserve isobaric labile placements and tolerance boundaries))
{
  const auto base = NASequence::fromString("AAAAAACG");
  auto loss = *db->getRibonucleotide("A"), gain = loss;
  loss.setCode("lossA");
  loss.setFormula(EmpiricalFormula("C9H11N5O4"));
  gain.setCode("gainA");
  gain.setFormula(EmpiricalFormula("C11H15N5O4"));
  for (auto* residue : {&loss, &gain})
  {
    residue->setMonoMass(residue->getFormula().getMonoWeight());
  }
  RNAModificationDiagnostics::Evidence evidence;
  evidence.tested = evidence.allowed = evidence.labile = {"lossA", "gainA"};
  const std::set<const Ribonucleotide*> modifications {&loss, &gain};
  auto target = base;
  target.set(0, &loss);
  target.set(3, &gain);
  target.setThreePrimeMod(db->getRibonucleotide("3'-p"));
  const double mass = target.getMonoWeight();
  for (bool ppm : {true, false})
  {
    RNARouteCandidateGenerator::Options options;
    options.precursor_tolerance_ppm = ppm;
    options.precursor_tolerance = ppm ? 20 : 0.02;
    const double tolerance = options.precursor_tolerance * (ppm ? mass * 1e-6 : 1);
    const std::vector<double> precursors {mass - tolerance + 1e-9, mass, mass + tolerance - 1e-9};
    const auto fast = RNARouteCandidateGenerator({base}, modifications, options).generate(theory(base), precursors, evidence, true, 1);
    options.index_precursor_mass = false;
    const auto reference = RNARouteCandidateGenerator({base}, modifications, options).generate(theory(base), precursors, evidence, true, 1);
    TEST_TRUE(fast.complete)
    TEST_TRUE(reference.complete)
    std::set<std::pair<Size, std::string>> left, right;
    for (const auto& candidate : fast.candidates)
    {
      left.emplace(candidate.precursor_index, candidate.sequence.toString());
    }
    for (const auto& candidate : reference.candidates)
    {
      right.emplace(candidate.precursor_index, candidate.sequence.toString());
    }
    TEST_TRUE(left == right)
    TEST_TRUE(left.size() > 100)
    for (Size precursor = 0; precursor < precursors.size(); ++precursor)
    {
      TEST_TRUE(left.contains({precursor, target.toString()}))
    }
    TEST_TRUE(fast.completion_nodes < reference.completion_nodes)
  }
}
END_SECTION

START_SECTION((route scoring respects assigned charges and bounds only unknown hypotheses))
{
  const auto sequence = NASequence::fromString("ACUGUACGp");
  const std::set<std::string> ions {"c", "y"};
  const auto model = theory(sequence, -1, ions);
  const auto query = RNATagCandidateFilter::prepare(sequence, nullptr, model, ions);
  auto observed = model;
  MSSpectrum::IntegerDataArray charges;
  charges.setName("NASE:legacy_fragment_charge");
  charges.assign(observed.size(), 2);
  observed.getIntegerDataArrays().push_back(charges);
  TEST_EQUAL(RNARouteScorer::score(query, observed, true, {}).score, 0)
  observed.getIntegerDataArrays().back().assign(observed.size(), 1);
  const auto prepared = RNARouteScorer::prepareSpectrum(observed, true, {});
  TEST_TRUE(RNARouteScorer::score(query, observed, prepared, true, {}).score > 0)
  TEST_EQUAL(RNARouteScorer::score(query, observed, prepared, true, {}).score, RNARouteScorer::score(query, observed, true, {}).score)
  const auto high_model = theory(sequence, -4, ions);
  std::vector<Size> high_indices;
  for (Size i = 0; i < high_model.size(); ++i)
  {
    if (std::abs(high_model.getIntegerDataArrays()[0][i]) == 4) { high_indices.push_back(i); }
  }
  observed = high_model;
  observed.select(high_indices);
  const auto high_query = RNATagCandidateFilter::prepare(sequence, nullptr, high_model, ions);
  RNARouteScorer::Options options;
  options.max_unknown_charge = 1;
  TEST_EQUAL(RNARouteScorer::score(high_query, observed, true, options).score, 0)
  charges.assign(observed.size(), 4);
  observed.getIntegerDataArrays().push_back(charges);
  TEST_TRUE(RNARouteScorer::score(high_query, observed, true, options).score > 0)
}
END_SECTION

START_SECTION((local background uses actual neutral masses and preserves selected charge metadata))
{
  MSSpectrum observed;
  const double proton = Constants::PROTON_MASS_U;
  observed.emplace_back(100, 1);
  observed.emplace_back((100 + proton) / 2 - proton, 10);
  observed.emplace_back(100.01, 100);
  MSSpectrum::IntegerDataArray charges;
  charges.setName("NASE:legacy_fragment_charge");
  charges.assign({1, 2, 2});
  observed.getIntegerDataArrays().push_back(charges);
  observed.sortByPosition();
  RNARouteScorer::Options options;
  options.local_radius = 1;
  const auto prepared = RNARouteScorer::prepareSpectrum(observed, true, options);
  const auto& local = prepared.backgrounds[1].at(1);
  TEST_TRUE(std::find(local.begin(), local.end(), 0) != local.end())
  TEST_TRUE(std::find(local.begin(), local.end(), 1) != local.end())
  TEST_TRUE(std::find(local.begin(), local.end(), 2) == local.end())
  observed.select(RNASpectrumGraph::selectPeakIndices(observed, 2));
  TEST_EQUAL(observed.size(), 2)
  TEST_EQUAL((*RNASpectrumGraph::legacyCharges(observed))[0], 2)
  TEST_EQUAL((*RNASpectrumGraph::legacyCharges(observed))[1], 2)
}
END_SECTION

START_SECTION((one experimental peak cannot support multiple route rungs through charge hypotheses))
{
  MSSpectrum observed;
  observed.emplace_back(NASequence::fromString("A").getMonoWeight(NASequence::CIon, -1), 1);
  const auto target = NASequence::fromString("AAAAGp");
  const auto result = RNARouteCandidateGenerator({target}, {}, {}).generate(observed, {target.getMonoWeight()}, {}, true, 3);
  TEST_TRUE(result.complete)
  TEST_EQUAL(result.routes, 0)
  TEST_TRUE(result.candidates.empty())
}
END_SECTION

START_SECTION((graph transitions follow native fragment masses while precursor validation uses formulas))
{
  const auto base = NASequence::fromString("AAACG");
  auto modification = *db->getRibonucleotide("A");
  modification.setCode("staleA");
  modification.setFormula(EmpiricalFormula("C11H15N5O4"));
  // The native fragment generator uses cached mono masses; NASequence precursor
  // masses use formulas. Preserve both conventions even when they disagree.
  auto target = base;
  target.set(0, &modification);
  target.set(1, &modification);
  target.set(2, &modification);
  target.setThreePrimeMod(db->getRibonucleotide("3'-p"));
  const auto result = RNARouteCandidateGenerator({base}, {&modification}, {}).generate(theory(target), {target.getMonoWeight()}, {}, true, 1);
  TEST_TRUE(result.complete)
  TEST_TRUE(std::any_of(result.candidates.begin(), result.candidates.end(), [&](const auto& candidate) { return candidate.sequence == target; }))
}
END_SECTION

START_SECTION((labile completion retains and scores the canonical generating scaffold))
{
  auto base = NASequence::fromString("ACUGUACGp");
  auto target = base;
  target.set(5, db->getRibonucleotide("t6A"));
  RNAModificationDiagnostics::Evidence evidence;
  evidence.tested = evidence.allowed = evidence.labile = {"t6A"};
  const auto observed = theory(base);
  const auto result
    = RNARouteCandidateGenerator({base}, {db->getRibonucleotide("t6A")}, {}).generate(observed, {target.getMonoWeight()}, evidence, true, 1);
  TEST_TRUE(result.complete)
  const auto found
    = std::find_if(result.candidates.begin(), result.candidates.end(), [&](const auto& candidate) { return candidate.sequence == target; });
  TEST_TRUE(found != result.candidates.end())
  TEST_TRUE(
    std::any_of(found->support_scaffolds.begin(), found->support_scaffolds.end(), [&](const auto& support) { return support->sequence == base; }))
  const auto prepared = RNARouteScorer::prepareSpectrum(observed, true, {});
  const auto intact = RNARouteScorer::match(RNATagCandidateFilter::prepare(target, nullptr, theory(target), {"c", "y"}), observed, prepared, true, {},
                                            target.toString());
  const auto canonical
    = RNARouteScorer::match(RNATagCandidateFilter::prepare(base, nullptr, theory(base), {"c", "y"}), observed, prepared, true, {}, base.toString());
  const auto complete = RNARouteScorer::score(std::vector<const RNARouteScorer::PreparedRoutes*> {&intact, &canonical}, observed, prepared, true);
  const auto partial = RNARouteScorer::score(std::vector<const RNARouteScorer::PreparedRoutes*> {&intact}, observed, prepared, true);
  TEST_TRUE(complete.nodes > partial.nodes)
  TEST_STRING_EQUAL(complete.prefix_scaffold, base.toString())
  TEST_STRING_EQUAL(complete.suffix_scaffold, base.toString())
  TEST_EQUAL(complete.score, RNARouteScorer::score(RNATagCandidateFilter::prepare(base, nullptr, theory(base), {"c", "y"}), observed, true, {}).score)
  RNARouteCandidateGenerator::Options limited;
  limited.max_supports = 1;
  const auto exhausted
    = RNARouteCandidateGenerator({base}, {db->getRibonucleotide("t6A")}, limited).generate(observed, {target.getMonoWeight()}, evidence, true, 1);
  TEST_FALSE(exhausted.complete)
  TEST_TRUE(exhausted.candidates.empty())
}
END_SECTION

START_SECTION((overlapping intact and canonical labile routes retain independent side evidence))
{
  const auto base = NASequence::fromString("ACAAUCUG");
  auto target = base;
  target.set(0, db->getRibonucleotide("mA"));
  target.set(3, db->getRibonucleotide("t6A"));
  target.set(6, db->getRibonucleotide("mU"));
  target.setThreePrimeMod(db->getRibonucleotide("3'-p"));
  auto loss = target;
  loss.set(3, db->getRibonucleotide("A"));
  MSSpectrum observed;
  for (const auto& [sequence, series] : std::vector<std::pair<NASequence, char>> {{loss, 'c'}, {target, 'y'}})
  {
    const auto model = theory(sequence);
    for (Size i = 0; i < model.size(); ++i)
    {
      const auto& name = model.getStringDataArrays()[0][i];
      if (name[0] == series && std::stoul(name.substr(1)) <= 5) { observed.push_back(model[i]); }
    }
  }
  observed.sortByPosition();
  RNAModificationDiagnostics::Evidence evidence;
  evidence.tested = evidence.allowed = evidence.labile = {"t6A"};
  const std::set<const Ribonucleotide*> mods {db->getRibonucleotide("mA"), db->getRibonucleotide("mU"), db->getRibonucleotide("t6A")};
  RNARouteCandidateGenerator::Options options;
  const auto fast = RNARouteCandidateGenerator({base}, mods, options).generate(observed, {target.getMonoWeight()}, evidence, true, 1);
  options.index_precursor_mass = options.compress_equivalent_routes = options.compress_route_states = false;
  const auto reference = RNARouteCandidateGenerator({base}, mods, options).generate(observed, {target.getMonoWeight()}, evidence, true, 1);
  TEST_TRUE(fast.complete)
  TEST_TRUE(reference.complete)
  const auto keys = [](const auto& candidates) {
    std::map<std::string, std::set<std::string>> result;
    for (const auto& candidate : candidates)
    {
      auto& models = result[candidate.sequence.toString()];
      for (const auto& scaffold : candidate.support_scaffolds)
      {
        models.insert(scaffold->label);
      }
    }
    return result;
  };
  TEST_TRUE(keys(fast.candidates) == keys(reference.candidates))
  const auto found
    = std::find_if(fast.candidates.begin(), fast.candidates.end(), [&](const auto& candidate) { return candidate.sequence == target; });
  TEST_TRUE(found != fast.candidates.end())
  TEST_TRUE(
    std::any_of(found->support_scaffolds.begin(), found->support_scaffolds.end(), [&](const auto& scaffold) { return scaffold->sequence == loss; }))
  const auto prepared = RNARouteScorer::prepareSpectrum(observed, true, {});
  const auto intact = RNARouteScorer::match(RNATagCandidateFilter::prepare(target, nullptr, theory(target), {"c", "y"}), observed, prepared, true, {},
                                            target.toString());
  const auto canonical
    = RNARouteScorer::match(RNATagCandidateFilter::prepare(loss, nullptr, theory(loss), {"c", "y"}), observed, prepared, true, {}, loss.toString());
  const auto combined = RNARouteScorer::score(std::vector<const RNARouteScorer::PreparedRoutes*> {&intact, &canonical}, observed, prepared, true);
  TEST_TRUE(combined.score > RNARouteScorer::score(std::vector<const RNARouteScorer::PreparedRoutes*> {&intact}, observed, prepared, true).score)
  TEST_TRUE(combined.score > RNARouteScorer::score(std::vector<const RNARouteScorer::PreparedRoutes*> {&canonical}, observed, prepared, true).score)
  TEST_EQUAL(combined.prefix_length, 5)
  TEST_EQUAL(combined.suffix_length, 5)
  TEST_STRING_EQUAL(combined.prefix_scaffold, loss.toString())
  TEST_STRING_EQUAL(combined.suffix_scaffold, target.toString())
}
END_SECTION

END_TEST
