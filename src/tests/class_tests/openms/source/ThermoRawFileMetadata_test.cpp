// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/HANDLERS/ThermoRawFileMetadata.h>
///////////////////////////

#include <optional>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::Internal;

// Precursor reconstruction from synthetic scan metadata. Needs neither a RAW file nor the
// vendor bridge, so it runs on every platform including builds without WITH_THERMO_RAW.

namespace
{
ThermoReaction reaction(double mass, const std::string& activation = "CollisionInducedDissociation", double width = 2.0, double offset = 0.0)
{
  ThermoReaction result;
  result.precursor_mass = mass;
  result.activation = activation;
  result.isolation_width = width;
  result.isolation_offset = offset;
  result.collision_energy = 35.0;
  result.collision_energy_valid = true;
  return result;
}

ThermoScan scan(int number, int level, const std::string& filter, std::vector<ThermoReaction> reactions = {})
{
  ThermoScan result;
  result.scan_number = number;
  result.ms_level = level;
  result.filter = filter;
  result.reactions = std::move(reactions);
  result.native_id = "controllerType=0 controllerNumber=1 scan=" + std::to_string(number);
  return result;
}

void trailer(ThermoScan& scan, const std::string& label, const std::string& value)
{
  scan.trailer.push_back({label, value});
}
} // namespace

START_TEST(ThermoRawFileMetadata, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static std::optional<double> number(const std::string& value)))
{
  TEST_FALSE(ThermoRawFileMetadata::number("").has_value())        // missing values remain missing
  TEST_FALSE(ThermoRawFileMetadata::number("NaN").has_value())     // nonfinite values remain missing
  TEST_FALSE(ThermoRawFileMetadata::number("2,5").has_value())     // do not guess decimal separators
  TEST_FALSE(ThermoRawFileMetadata::number("12.3abc").has_value()) // trailing garbage is not a number
  TEST_TRUE(ThermoRawFileMetadata::number(" 0 ").has_value())      // zero is not missing
  TEST_EQUAL(ThermoRawFileMetadata::number(" 0 ").value(), 0.0)
  TEST_EQUAL(ThermoRawFileMetadata::number("499.5").value(), 499.5)
}
END_SECTION

START_SECTION((static double selectedIon(double target, std::optional<double> mono, double width)))
{
  TEST_EQUAL(ThermoRawFileMetadata::selectedIon(500, 499.5, 2), 499.5)       // retain plausible monoisotopic ion
  TEST_EQUAL(ThermoRawFileMetadata::selectedIon(500, 490, 2), 500)           // reject unrelated trailer monoisotopic mass
  TEST_EQUAL(ThermoRawFileMetadata::selectedIon(500, 497, 10), 497)          // wide-window selected ion rule
  TEST_EQUAL(ThermoRawFileMetadata::selectedIon(500, std::nullopt, 2), 500)  // missing monoisotopic mass
  TEST_EQUAL(ThermoRawFileMetadata::selectedIon(500, 0.0, 2), 500)           // unset (zero) monoisotopic mass
}
END_SECTION

START_SECTION((std::vector<ThermoPrecursor> precursors(const ThermoScan& scan)))
{
  ThermoRawFileMetadata parser;
  TEST_TRUE(parser.precursors(scan(10, 1, "FTMS + p NSI Full ms [200-2000]")).empty()) // MS1 has no precursors

  // MS2 with supplemental activation and trailer overrides for charge, monoisotopic m/z and width
  auto ms2 = scan(11, 2, "FTMS + c NSI Full ms2 500.00@etd35.00 500.00@hcd20.00 [100-1000]",
                  {reaction(500, "ElectronTransferDissociation", 2, 0.25), reaction(500, "HigherEnergyCollisionalDissociation")});
  trailer(ms2, "Master Scan Number:", "10");
  trailer(ms2, "Monoisotopic M/Z:", "499.5");
  trailer(ms2, "Charge State:", "3");
  trailer(ms2, "MS2 Isolation Width:", "4");
  const auto p2 = parser.precursors(ms2);
  TEST_EQUAL(p2.size(), 1) // supplemental activation is not an extra precursor
  ABORT_IF(p2.empty())
  TEST_EQUAL(p2[0].target_mz, 500)     // separate selected ion and isolation target
  TEST_EQUAL(p2[0].selected_mz, 499.5)
  TEST_TRUE(p2[0].width.has_value())   // trailer width overrides the reaction width
  TEST_EQUAL(p2[0].width.value(), 4)
  TEST_EQUAL(p2[0].lower_offset, 1.75) // preserve asymmetric window
  TEST_EQUAL(p2[0].upper_offset, 2.25)
  TEST_TRUE(p2[0].charge.has_value())
  TEST_EQUAL(p2[0].charge.value(), 3)
  TEST_EQUAL(p2[0].parent_scan, 10)
  TEST_STRING_EQUAL(p2[0].spectrum_ref, "controllerType=0 controllerNumber=1 scan=10")
  TEST_TRUE(p2[0].estimate_intensity)
  TEST_STRING_EQUAL(p2[0].activation.activation, "ElectronTransferDissociation")
  TEST_TRUE(p2[0].supplemental.has_value()) // preserve supplemental HCD
  TEST_STRING_EQUAL(p2[0].supplemental.value().activation, "HigherEnergyCollisionalDissociation")

  // MS3 with modern SPS trailer labels; lists its own selections and then the MS2 ancestor
  auto ms3 = scan(12, 3, "ITMS + c NSI Full ms3 500.00@etd35.00 500.00@hcd20.00 300.00@cid35.00 [80-600]",
                  {reaction(500, "ElectronTransferDissociation"), reaction(500, "HigherEnergyCollisionalDissociation"), reaction(300)});
  trailer(ms3, "Master Scan Number:", "11");
  trailer(ms3, "SPS Masses:", "300, 310,");
  trailer(ms3, "SPS Masses Continued:", "320");
  const auto p3 = parser.precursors(ms3);
  TEST_EQUAL(p3.size(), 4) // primary selection, two SPS selections, and the MS2 ancestor
  ABORT_IF(p3.size() != 4)
  TEST_EQUAL(p3[0].target_mz, 300) // use the current MS3 reaction, not reaction zero
  TEST_EQUAL(p3[0].parent_scan, 11)
  TEST_EQUAL(p3[1].target_mz, 310) // continued SPS list preserved
  TEST_EQUAL(p3[2].target_mz, 320)
  TEST_TRUE(p3[0].estimate_intensity) // only the primary SPS selection receives an intensity estimate
  TEST_FALSE(p3[1].estimate_intensity)
  TEST_FALSE(p3[2].estimate_intensity)
  TEST_EQUAL(p3[3].target_mz, 500) // ancestor references its own parent
  TEST_EQUAL(p3[3].parent_scan, 10)
  TEST_FALSE(p3[0].charge.has_value()) // missing charge remains missing
}
END_SECTION

START_SECTION([EXTRA] invalid master scan falls back to the filter hierarchy)
{
  ThermoRawFileMetadata parser;
  parser.precursors(scan(1, 1, "FTMS + p NSI Full ms [200-2000]"));
  const auto simple = parser.precursors(scan(2, 2, "ITMS + c NSI Full ms2 500.00@cid35.00 [100-1000]", {reaction(500)}));
  TEST_EQUAL(simple.size(), 1)
  ABORT_IF(simple.empty())
  TEST_EQUAL(simple[0].parent_scan, 1) // MS2 parent from the most recent MS1

  auto third = scan(3, 3, "ITMS + c NSI Full ms3 500.00@cid35.00 300.00@cid35.00 [80-600]", {reaction(500), reaction(300)});
  trailer(third, "Master Scan Number:", "999");
  trailer(third, "SPS Mass 1:", "300");
  trailer(third, "SPS Mass 2:", "310");
  const auto p = parser.precursors(third);
  TEST_EQUAL(p.size(), 3) // legacy SPS labels retained
  ABORT_IF(p.empty())
  TEST_EQUAL(p[0].parent_scan, 2) // invalid master falls back to the filter hierarchy
  TEST_EQUAL(p[0].target_mz, 300)
}
END_SECTION

START_SECTION([EXTRA] truncated file starts at the last reaction)
{
  ThermoRawFileMetadata parser;
  const auto orphan = scan(100, 3, "ITMS + c NSI Full ms3 500.00@cid35.00 300.00@cid35.00 [80-600]",
                           {reaction(500), reaction(300, "CollisionInducedDissociation", -1)});
  const auto op = parser.precursors(orphan);
  TEST_EQUAL(op.size(), 1)
  ABORT_IF(op.empty())
  TEST_EQUAL(op[0].target_mz, 300)
  TEST_EQUAL(op[0].parent_scan, 0)          // no fabricated parent
  TEST_TRUE(op[0].spectrum_ref.empty())
  TEST_FALSE(op[0].width.has_value())       // negative width is unknown
  TEST_EQUAL(op[0].lower_offset, 0.0)
  TEST_EQUAL(op[0].upper_offset, 0.0)
}
END_SECTION

START_SECTION([EXTRA] same-order master scan is a sibling)
{
  // Lumos/Tribrid files can point the master scan of an MS2 at a sibling MS2 of the same
  // precursor (e.g. HCD scan first, then EThcD). The sibling's reaction must not be consumed.
  ThermoRawFileMetadata parser;
  parser.precursors(scan(10, 1, "FTMS + p NSI Full ms [200-2000]"));
  auto sibling = scan(14, 2, "FTMS + c ESI d Full ms2 432.90@hcd30.00 [100-1300]", {reaction(432.9, "HigherEnergyCollisionalDissociation")});
  trailer(sibling, "Master Scan Number:", "10");
  const auto ps = parser.precursors(sibling);
  TEST_EQUAL(ps.size(), 1)
  ABORT_IF(ps.empty())
  TEST_EQUAL(ps[0].parent_scan, 10) // sibling descends from the MS1

  auto ethcd = scan(15, 2, "FTMS + c ESI d sa Full ms2 432.90@etd54.00 432.90@hcd30.00 [100-1300]",
                    {reaction(432.9, "ElectronTransferDissociation"), reaction(432.9, "HigherEnergyCollisionalDissociation")});
  trailer(ethcd, "Master Scan Number:", "14");
  const auto pe = parser.precursors(ethcd);
  TEST_EQUAL(pe.size(), 1)
  ABORT_IF(pe.empty())
  TEST_STRING_EQUAL(pe[0].activation.activation, "ElectronTransferDissociation") // same-order master does not consume the primary reaction
  TEST_TRUE(pe[0].supplemental.has_value()) // supplemental HCD retained
  TEST_STRING_EQUAL(pe[0].supplemental.value().activation, "HigherEnergyCollisionalDissociation")
  TEST_EQUAL(pe[0].parent_scan, 10) // resolve the sibling's own parent
  TEST_STRING_EQUAL(pe[0].spectrum_ref, "controllerType=0 controllerNumber=1 scan=10")
}
END_SECTION

START_SECTION([EXTRA] supplemental reaction retained for unexpected main activation types)
{
  ThermoRawFileMetadata parser;
  parser.precursors(scan(10, 1, "FTMS + p NSI Full ms [200-2000]"));
  auto unusual = scan(13, 2, "ITMS + c NSI Full ms2 500.00@cid35.00 500.00@hcd20.00 [100-1000]",
                      {reaction(500), reaction(500, "HigherEnergyCollisionalDissociation")});
  trailer(unusual, "Master Scan Number:", "10");
  const auto pu = parser.precursors(unusual);
  TEST_EQUAL(pu.size(), 1)
  ABORT_IF(pu.empty())
  TEST_STRING_EQUAL(pu[0].activation.activation, "CollisionInducedDissociation")
  TEST_TRUE(pu[0].supplemental.has_value())
  TEST_STRING_EQUAL(pu[0].supplemental.value().activation, "HigherEnergyCollisionalDissociation")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
