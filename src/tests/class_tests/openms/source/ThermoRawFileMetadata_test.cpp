// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Deliberately independent of libOpenMS and vendor binaries: precursor
// reconstruction can be tested on all platforms, including builds without
// WITH_THERMO_RAW.
#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include <OpenMS/FORMAT/HANDLERS/ThermoRawFileMetadata.h>
#include <iostream>
#include <stdexcept>

using Json = nlohmann::json;
using Metadata = OpenMS::Internal::ThermoRawFileMetadata;

namespace
{
void check(bool condition, const char* description)
{
  if (! condition) { throw std::runtime_error(description); }
}

Json reaction(double mass, const std::string& activation = "CollisionInducedDissociation", double width = 2.0, double offset = 0.0)
{
  return {{"precursor_mass", mass},     {"activation", activation}, {"isolation_width", width},
          {"isolation_offset", offset}, {"collision_energy", 35.0}, {"collision_energy_valid", true}};
}

Json scan(int number, int level, const std::string& filter, const Json& reactions = Json::array())
{
  return {{"scan_number", number},  {"ms_level", level},        {"filter", filter},
          {"reactions", reactions}, {"trailer", Json::array()}, {"native_id", "controllerType=0 controllerNumber=1 scan=" + std::to_string(number)}};
}

void trailer(Json& scan, const std::string& label, const std::string& value)
{ scan["trailer"].push_back({{"label", label}, {"value", value}}); }
} // namespace

int main()
{
  try
  {
    check(Metadata::number("").is_null(), "Missing values remain missing");
    check(Metadata::number("NaN").is_null(), "Nonfinite values remain missing");
    check(Metadata::number("2,5").is_null(), "Do not guess decimal separators");
    check(Metadata::number(" 0 ") == 0.0, "Zero is not missing");
    check(Metadata::selectedIon(500, 499.5, 2) == 499.5, "Retain plausible monoisotopic ion");
    check(Metadata::selectedIon(500, 490, 2) == 500, "Reject unrelated trailer monoisotopic mass");
    check(Metadata::selectedIon(500, 497, 10) == 497, "Wide-window selected ion rule");

    Metadata parser;
    check(parser.precursors(scan(10, 1, "FTMS + p NSI Full ms [200-2000]")).empty(), "MS1 has no precursors");
    auto ms2 = scan(11, 2, "FTMS + c NSI Full ms2 500.00@etd35.00 500.00@hcd20.00 [100-1000]",
                    Json::array({reaction(500, "ElectronTransferDissociation", 2, 0.25), reaction(500, "HigherEnergyCollisionalDissociation")}));
    trailer(ms2, "Master Scan Number:", "10");
    trailer(ms2, "Monoisotopic M/Z:", "499.5");
    trailer(ms2, "Charge State:", "3");
    trailer(ms2, "MS2 Isolation Width:", "4");
    auto p2 = parser.precursors(ms2);
    check(p2.size() == 1, "Supplemental activation is not an extra precursor");
    check(p2[0]["target_mz"] == 500 && p2[0]["selected_mz"] == 499.5, "Separate selected ion and isolation target");
    check(p2[0]["lower_offset"] == 1.75 && p2[0]["upper_offset"] == 2.25, "Preserve asymmetric window and trailer width override");
    check(p2[0]["charge"] == 3 && p2[0]["parent_scan"] == 10, "Preserve charge and parent scan");
    check(p2[0]["supplemental"]["activation"] == "HigherEnergyCollisionalDissociation", "Preserve supplemental HCD");

    auto ms3
      = scan(12, 3,
             "ITMS + c NSI Full ms3 500.00@etd35.00 500.00@hcd20.00 "
             "300.00@cid35.00 [80-600]",
             Json::array({reaction(500, "ElectronTransferDissociation"), reaction(500, "HigherEnergyCollisionalDissociation"), reaction(300)}));
    trailer(ms3, "Master Scan Number:", "11");
    trailer(ms3, "SPS Masses:", "300, 310,");
    trailer(ms3, "SPS Masses Continued:", "320");
    auto p3 = parser.precursors(ms3);
    check(p3.size() == 4, "Primary selection, two SPS selections, and MS2 ancestor");
    check(p3[0]["target_mz"] == 300 && p3[0]["parent_scan"] == 11, "Use current MS3 reaction, not reaction zero");
    check(p3[1]["target_mz"] == 310 && p3[2]["target_mz"] == 320, "Continued SPS list preserved");
    check(p3[0]["estimate_intensity"] == true && p3[1]["estimate_intensity"] == false,
          "Only the primary SPS selection receives an intensity estimate");
    check(p3[3]["target_mz"] == 500 && p3[3]["parent_scan"] == 10, "Ancestor references its own parent");
    check(p3[0]["charge"].is_null(), "Missing charge remains missing");

    Metadata fallback;
    fallback.precursors(scan(1, 1, "FTMS + p NSI Full ms [200-2000]"));
    auto simple = scan(2, 2, "ITMS + c NSI Full ms2 500.00@cid35.00 [100-1000]", Json::array({reaction(500)}));
    check(fallback.precursors(simple)[0]["parent_scan"] == 1, "MS2 parent from most recent MS1");
    auto third = scan(3, 3, "ITMS + c NSI Full ms3 500.00@cid35.00 300.00@cid35.00 [80-600]", Json::array({reaction(500), reaction(300)}));
    trailer(third, "Master Scan Number:", "999");
    trailer(third, "SPS Mass 1:", "300");
    trailer(third, "SPS Mass 2:", "310");
    auto p = fallback.precursors(third);
    check(p[0]["parent_scan"] == 2 && p.size() == 3, "Invalid master falls back to filter hierarchy; legacy SPS retained");

    Metadata isolated;
    auto orphan = scan(100, 3, "ITMS + c NSI Full ms3 500.00@cid35.00 300.00@cid35.00 [80-600]",
                       Json::array({reaction(500), reaction(300, "CollisionInducedDissociation", -1)}));
    auto op = isolated.precursors(orphan);
    check(op[0]["target_mz"] == 300, "Truncated file starts at last reaction");
    check(op[0]["spectrum_ref"] == "" && op[0]["width"].is_null(), "No fabricated parent or negative width");

    auto unusual = scan(13, 2, "ITMS + c NSI Full ms2 500.00@cid35.00 500.00@hcd20.00 [100-1000]",
                        Json::array({reaction(500), reaction(500, "HigherEnergyCollisionalDissociation")}));
    trailer(unusual, "Master Scan Number:", "10");
    check(parser.precursors(unusual)[0]["supplemental"].is_object(), "Retain supplemental reactions with unexpected main activation types");
    std::cout << "Thermo precursor metadata regression tests passed\n";
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
