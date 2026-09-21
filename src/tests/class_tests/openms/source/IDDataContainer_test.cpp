// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <concepts>
#include <iterator>
#include <utility>

using namespace OpenMS;
using namespace OpenMS::IdentificationDataInternal;

// Some record types deliberately provide no value equality. Instantiating their
// containers (including DLL exports on MSVC) must not instantiate std::equal.
static_assert(! std::equality_comparable<InputFile>);
static_assert(! std::equality_comparable<InputFiles>);
static_assert(! std::equality_comparable<ParentGroup>);
static_assert(! std::equality_comparable<ParentGroups>);
static_assert(std::equality_comparable<AppliedProcessingSteps>);

START_TEST(IDDataContainer, "$Id$")

START_SECTION((equality compares record values))
{
  AppliedProcessingSteps steps;
  AppliedProcessingSteps other;
  TEST_TRUE(steps == other)
  steps.push_back(AppliedProcessingStep());
  TEST_FALSE(steps == other)
  other = steps;
  TEST_TRUE(steps == other)

  IdentificationData data;
  auto software = data.registerProcessingSoftware(ProcessingSoftware("test", "1.0"));
  auto step = data.registerProcessingStep(ProcessingStep(software));
  other.clear();
  other.push_back(AppliedProcessingStep(step));
  TEST_FALSE(steps == other)
}
END_SECTION

START_SECTION((ordered uniqueness, stable references, copy and move))
{
  InputFiles records;
  TEST_TRUE(records.begin() == records.end())
  auto a = records.insert(InputFile("a")).first;
  const auto* address = &*a;
  records.insert(InputFile("c"));
  records.insert(InputFile("b"));
  TEST_EQUAL(records.size(), 3)
  TEST_FALSE(records.insert(InputFile("a")).second)
  TEST_TRUE(&*records.find("a") == address)
  TEST_EQUAL(std::distance(records.begin(), records.end()), 3)
  TEST_EQUAL((--records.end())->name, "c")

  auto copied = records;
  TEST_TRUE(&*copied.find("a") != address)
  TEST_TRUE(records.modify(a, [](InputFile& value) { value.name = "d"; }))
  TEST_TRUE(&*a == address)
  TEST_EQUAL(a->name, "d")
  TEST_TRUE(a == --records.end())
  TEST_TRUE(copied.find("a") != copied.end())

  auto moved = std::move(records);
  TEST_TRUE(&*moved.find("d") == address)
  TEST_TRUE(++a == moved.end())
  records.insert(InputFile("reused"));
  TEST_EQUAL(records.size(), 1)
  TEST_EQUAL(moved.size(), 3)

  InputFiles swapped;
  auto b = moved.find("b");
  swapped.swap(moved);
  TEST_EQUAL(b->name, "b")
  TEST_EQUAL((++b)->name, "c")
  TEST_TRUE(moved.empty())
}
END_SECTION

START_SECTION((conflicting modifications erase only the modified record))
{
  InputFiles records;
  auto first = records.insert(InputFile("first")).first;
  auto second = records.insert(InputFile("second")).first;
  TEST_FALSE(records.modify(first, [](InputFile& value) { value.name = "second"; }))
  TEST_EQUAL(records.size(), 1)
  TEST_TRUE(records.begin() == second)
  TEST_TRUE(records.erase(second) == records.end())
}
END_SECTION

START_SECTION((composite - key ranges and processing - step order))
{
  IdentificationData data;
  auto file = data.registerInputFile(InputFile("file"));
  auto first = data.registerObservation(Observation("first", file));
  auto second = data.registerObservation(Observation("second", file));
  auto peptide = data.registerIdentifiedPeptide(IdentifiedPeptide(AASequence::fromString("PEPTIDE")));
  data.registerObservationMatch(ObservationMatch(peptide, first));
  data.registerObservationMatch(ObservationMatch(peptide, second));
  auto range = data.getObservationMatches().equal_range(first);
  TEST_EQUAL(std::distance(range.first, range.second), 1)
  TEST_TRUE(range.first->observation_ref == first)

  AppliedProcessingSteps steps;
  auto software = data.registerProcessingSoftware(ProcessingSoftware("test", "1.0"));
  auto step = data.registerProcessingStep(ProcessingStep(software));
  steps.push_back(AppliedProcessingStep(step));
  steps.push_back(AppliedProcessingStep());
  TEST_FALSE(steps.push_back(AppliedProcessingStep()).second)
  TEST_EQUAL(steps.size(), 2)
  TEST_TRUE(steps.get<1>().find(std::nullopt) != steps.get<1>().end())
  TEST_TRUE(steps.begin()->processing_step_opt == step)
  TEST_TRUE(steps.rbegin()->processing_step_opt == std::nullopt)
  const auto& const_steps = steps;
  auto ordered = const_steps.get<1>();
  TEST_EQUAL(std::distance(ordered.begin(), ordered.end()), 2)
  auto position = ordered.begin();
  TEST_TRUE(position->processing_step_opt == std::nullopt)
  TEST_TRUE((++position)->processing_step_opt == step)
  TEST_TRUE(++position == ordered.end())
  TEST_TRUE((--position)->processing_step_opt == step)
  TEST_TRUE((--position)->processing_step_opt == std::nullopt)
  auto step_range = steps.equal_range(step);
  TEST_EQUAL(std::distance(step_range.first, step_range.second), 1)
}
END_SECTION

END_TEST
