// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch, Timo Sachsenberg, Chris Bielow, Jang Jang Jin$
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ElementDB.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <iostream>
#include <cmath>
#include <memory>

using namespace std;

namespace OpenMS
{
  ElementDB::ElementDB()
  {
    storeElements_();
  }

  ElementDB::~ElementDB()
  {
    clear_();
  }

  ElementDB* ElementDB::getInstance()
  {
    static ElementDB* db_ = new ElementDB;
    return db_;
  }

  const unordered_map<string, const Element*>& ElementDB::getNames() const
  {
    return names_;
  }

  const unordered_map<string, const Element*>& ElementDB::getSymbols() const
  {
    return symbols_;
  }

  const unordered_map<unsigned int, const Element*>& ElementDB::getAtomicNumbers() const
  {
    return atomic_numbers_;
  }

  const Element* ElementDB::getElement(const string& name) const
  {
    if (auto entry = symbols_.find(name); entry != symbols_.end())
    {
      return entry->second;
    }
    else
    {
      if (auto entry = names_.find(name); entry != names_.end())
      {
        return entry->second;
      }
    }
    return nullptr;
  }

  const Element* ElementDB::getElement(unsigned int atomic_number) const
  {
    if (auto entry = atomic_numbers_.find(atomic_number); entry != atomic_numbers_.end())
    {
      return entry->second;
    }
    return nullptr;
  }

  bool ElementDB::hasElement(const string& name) const
  {
    return (names_.count(name) == 1)|| (symbols_.count(name) == 1);
  }

  bool ElementDB::hasElement(unsigned int atomic_number) const
  {
    return atomic_numbers_.find(atomic_number) != atomic_numbers_.end();
  }

  void ElementDB::addElement(const std::string& name,
                             const std::string& symbol,
                             const unsigned int an,
                             const std::map<unsigned int, double>& abundance,
                             const std::map<unsigned int, double>& mass,
                             bool replace_existing)
  {
    if (hasElement(an) && !replace_existing)
    {      
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Element with atomic number ") + an + " already exists");
    }
    buildElement_(name, symbol, an, abundance, mass);
  }

  double ElementDB::calculateAvgWeight_(const map<unsigned int, double>& abundance, const map<unsigned int, double>& mass)
  {
    double avg = 0;
    // calculate weighted average
    for (const auto& it : abundance)
    {
      avg += mass.at(it.first) * abundance.at(it.first);
    }
    return avg;
  }

  double ElementDB::calculateMonoWeight_(const map<unsigned int, double>& abundance, const map<unsigned int, double>& mass)
  {
    double highest_abundance = -1.0;
    int highest_abundance_isotope = -1;

    // the monoisotopic weight is the *most abundant* isotope of an element
    for (const auto& it : abundance)
    {
      if (it.second > highest_abundance)
      {
        highest_abundance = it.second;
        highest_abundance_isotope = it.first;
      }
    }

    if (highest_abundance_isotope != -1) return mass.at(highest_abundance_isotope);
    else return 0.0;
  }

  
  template<class CONT, class KEY>
  void addIfUniqueOrThrow(CONT& container, const KEY& key, unique_ptr<const Element>& replacement)
  {
    auto elem = container.find(key);
    if (elem != container.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(key), "Already exists!");
    }
    container[key] = replacement.get();
  }

  void ElementDB::storeElements_()
  {	
    map<unsigned int, double> hydrogen_abundance = {{1u, 0.999885}, {2u, 0.000115}, {3u, 0.0}};
    map<unsigned int, double> hydrogen_mass = {{1u, 1.0078250319}, {2u, 2.01410178}, {3u, 3.01604927}};
    buildElement_("Hydrogen", "H", 1u, hydrogen_abundance, hydrogen_mass);


    map<unsigned int, double> helium_abundance = {{3u, 1.34e-06}, {4u, 0.9999986599999999}};
    map<unsigned int, double> helium_mass = {{3u, 3.0160293191}, {4u, 4.00260325415}};
    buildElement_("Helium", "He", 2u, helium_abundance, helium_mass);


    map<unsigned int, double> lithium_abundance = {{6u, 0.0759}, {7u, 0.9240999999999999}};
    map<unsigned int, double> lithium_mass = {{6u, 6.015122}, {7u, 7.016004}};
    buildElement_("Lithium", "Li", 3u, lithium_abundance, lithium_mass);


    map<unsigned int, double> beryllium_abundance = {{9u, 1.0}};
    map<unsigned int, double> beryllium_mass = {{9u, 9.0121822}};
    buildElement_("Beryllium", "Be", 4u, beryllium_abundance, beryllium_mass);


    map<unsigned int, double> bor_abundance = {{10u, 0.19899999999999998}, {11u, 0.8009999999999999}};
    map<unsigned int, double> bor_mass = {{10u, 10.012937000000001}, {11u, 11.009304999999999}};
    buildElement_("Boron", "B", 5u, bor_abundance, bor_mass);


    map<unsigned int, double> carbon_abundance = {{12u, 0.9893000000000001}, {13u, 0.010700000000000001}};
    map<unsigned int, double> carbon_mass = {{12u, 12.0}, {13u, 13.003355000000001}};
    buildElement_("Carbon", "C", 6u, carbon_abundance, carbon_mass);


    map<unsigned int, double> nitrogen_abundance = {{14u, 0.9963200000000001}, {15u, 0.00368}};
    map<unsigned int, double> nitrogen_mass = {{14u, 14.003074}, {15u, 15.000109}};
    buildElement_("Nitrogen", "N", 7u, nitrogen_abundance, nitrogen_mass);


    map<unsigned int, double> oxygen_abundance = {{16u, 0.9975700000000001}, {17u, 0.00037999999999999997}, {18u, 0.0020499999999999997}};
    map<unsigned int, double> oxygen_mass = {{16u, 15.994915000000001}, {17u, 16.999132}, {18u, 17.999168999999998}};
    buildElement_("Oxygen", "O", 8u, oxygen_abundance, oxygen_mass);


    map<unsigned int, double> fluorine_abundance = {{19u, 1.0}};
    map<unsigned int, double> fluorine_mass = {{19u, 18.99840322}};
    buildElement_("Fluorine", "F", 9u, fluorine_abundance, fluorine_mass);


    map<unsigned int, double> neon_abundance = {{20u, 0.9048},  {21u, 0.0027}, {22u, 0.0925}};
    map<unsigned int, double> neon_mass = {{20u,  19.99244018}, {21u, 20.9938467}, {22u, 21.9913851}};
    buildElement_("Neon", "Ne", 10u, neon_abundance, neon_mass);


    map<unsigned int, double> sodium_abundance = {{23u, 1.0}};
    map<unsigned int, double> sodium_mass = {{23u, 22.989769280899999}};
    buildElement_("Sodium", "Na", 11u, sodium_abundance, sodium_mass);


    map<unsigned int, double> magnesium_abundance = {{24u, 0.7898999999999999}, {25u, 0.1}, {26u, 0.1101}};
    map<unsigned int, double> magnesium_mass = {{24u, 23.985042}, {25u, 24.985837}, {26u, 25.982593000000001}};
    buildElement_("Magnesium", "Mg", 12u, magnesium_abundance, magnesium_mass);


    map<unsigned int, double> aluminium_abundance = {{27u, 1.0}};
    map<unsigned int, double> aluminium_mass = {{27u, 26.981538629999999}};
    buildElement_("Aluminium", "Al", 13u, aluminium_abundance, aluminium_mass);


    map<unsigned int, double> silicon_abundance = {{28u, 0.9220999999999999}, {29u, 0.0467}, {30u, 0.031}};
    map<unsigned int, double> silicon_mass = {{28u, 27.976926532499999}, {29u, 28.9764947}, {30u, 29.973770170000002}};
    buildElement_("Silicon", "Si", 14u, silicon_abundance, silicon_mass);


    map<unsigned int, double> phosphorus_abundance = {{31u, 1.0}};
    map<unsigned int, double> phosphorus_mass = {{31u, 30.973761490000001}};
    buildElement_("Phosphorus", "P", 15u, phosphorus_abundance, phosphorus_mass);


    map<unsigned int, double> sulfur_abundance = {{32u, 0.9493}, {33u, 0.0076}, {34u, 0.0429}, {36u, 0.0002}};
    map<unsigned int, double> sulfur_mass = {{32u, 31.972070729999999}, {33u, 32.971457999999998}, {34u, 33.967866999999998}, {36u, 35.967081}};
    buildElement_("Sulfur", "S", 16u, sulfur_abundance, sulfur_mass);


    map<unsigned int, double> chlorine_abundance = {{35u, 0.7576}, {37u, 0.24239999999999998}};
    map<unsigned int, double> chlorine_mass = {{35u, 34.968852679999998}, {37u, 36.965902589999999}};
    buildElement_("Chlorine", "Cl", 17u, chlorine_abundance, chlorine_mass);


    map<unsigned int, double> argon_abundance = {{36u, 0.003336}, {38u, 0.000629}, {40u, 0.996035}};
    map<unsigned int, double> argon_mass = {{36u, 35.967545106000003}, {38u, 37.9627324}, {40u, 39.9623831225}};
    buildElement_("Argon", "Ar", 18u, argon_abundance, argon_mass);


    map<unsigned int, double> potassium_abundance = {{39u, 0.932581}, {40u, 0.000117}, {41u, 0.067302}};
    map<unsigned int, double> potassium_mass = {{39u, 38.963706680000001}, {40u, 39.963998480000001}, {41u, 40.961825760000004}};
    buildElement_("Potassium", "K", 19u, potassium_abundance, potassium_mass);


    map<unsigned int, double> calcium_abundance = {{40u, 0.96941}, {42u, 0.00647}, {43u, 0.00135}, {44u, 0.02086}, {46u, 4e-05}, {48u, 0.00187}};
    map<unsigned int, double> calcium_mass = {{40u, 39.962590980000002}, {42u, 41.958618010000002}, {43u, 42.958766599999997}, {44u, 43.955481800000001}, {46u, 45.953692599999997}, {48u, 47.952534}};
    buildElement_("Calcium", "Ca", 20u, calcium_abundance, calcium_mass);


    map<unsigned int, double> scandium_abundance = {{45u, 1.0}};
    map<unsigned int, double> scandium_mass = {{45u, 44.955910000000003}};
    buildElement_("Scandium", "Sc", 21u, scandium_abundance, scandium_mass);


    map<unsigned int, double> titanium_abundance = {{46u, 0.0825}, {47u, 0.07440000000000001}, {48u, 0.7372}, {49u, 0.0541}, {50u, 0.0518}};
    map<unsigned int, double> titanium_mass = {{46u, 45.952631599999997}, {47u, 46.951763100000001}, {48u, 47.947946299999998}, {49u, 48.947870000000002}, {50u, 49.944791199999997}};
    buildElement_("Titanium", "Ti", 22u, titanium_abundance, titanium_mass);


    map<unsigned int, double> vanadium_abundance = {{50u, 0.0025}, {51u, 0.9975}};
    map<unsigned int, double> vanadium_mass = {{50u, 49.947158500000001}, {51u, 50.943959499999998}};
    buildElement_("Vanadium", "V", 23u, vanadium_abundance, vanadium_mass);


    map<unsigned int, double> chromium_abundance = {{50u, 0.043449999999999996}, {52u, 0.83789}, {53u, 0.09501}, {54u, 0.02365}};
    map<unsigned int, double> chromium_mass = {{50u, 49.946044200000003}, {52u, 51.940507500000003}, {53u, 52.940649399999998}, {54u, 53.938880400000002}};
    buildElement_("Chromium", "Cr", 24u, chromium_abundance, chromium_mass);


    map<unsigned int, double> manganese_abundance = {{55u, 1.0}};
    map<unsigned int, double> manganese_mass = {{55u, 54.938049999999997}};
    buildElement_("Manganese", "Mn", 25u, manganese_abundance, manganese_mass);


    map<unsigned int, double> ferrum_abundance = {{54u, 0.058449999999999995}, {56u, 0.91754}, {57u, 0.021191}, {58u, 0.002819}};
    map<unsigned int, double> ferrum_mass = {{54u, 53.939610500000001}, {56u, 55.934937499999997}, {57u, 56.935394000000002}, {58u, 57.933275600000002}};
    buildElement_("Ferrum", "Fe", 26u, ferrum_abundance, ferrum_mass);


    map<unsigned int, double> cobalt_abundance = {{59u, 1.0}};
    map<unsigned int, double> cobalt_mass = {{59u, 58.933194999999998}};
    buildElement_("Cobalt", "Co", 27u, cobalt_abundance, cobalt_mass);


    map<unsigned int, double> nickel_abundance = {{58u, 0.680169}, {60u, 0.262231}, {61u, 0.011399}, {62u, 0.036345}, {64u, 0.009256}};
    map<unsigned int, double> nickel_mass = {{58u, 57.935347999999998}, {60u, 59.930790999999999}, {61u, 60.931060000000002}, {62u, 61.928348999999997}, {64u, 63.927970000000002}};
    buildElement_("Nickel", "Ni", 28u, nickel_abundance, nickel_mass);


    map<unsigned int, double> copper_abundance = {{63u, 0.6917}, {65u, 0.30829999999999996}};
    map<unsigned int, double> copper_mass = {{63u, 62.929600999999998}, {65u, 64.927794000000006}};
    buildElement_("Copper", "Cu", 29u, copper_abundance, copper_mass);


    map<unsigned int, double> zinc_abundance = {{64u, 0.4863}, {66u, 0.27899999999999997}, {67u, 0.040999999999999995}, {68u, 0.1875}, {70u, 0.0062}};
    map<unsigned int, double> zinc_mass = {{64u, 63.929147}, {66u, 65.926036999999994}, {67u, 66.927131000000003}, {68u, 67.924847999999997}, {70u, 69.925325000000001}};
    buildElement_("Zinc", "Zn", 30u, zinc_abundance, zinc_mass);


    map<unsigned int, double> gallium_abundance = {{69u, 0.60108}, {71u, 0.39892000000000005}};
    map<unsigned int, double> gallium_mass = {{69u, 68.925573600000007}, {71u, 70.924701299999995}};
    buildElement_("Gallium", "Ga", 31u, gallium_abundance, gallium_mass);


    map<unsigned int, double> germanium_abundance = {{70u, 0.20379999999999998}, {72u, 0.2731}, {73u, 0.0776}, {74u, 0.36719999999999997}, {76u, 0.0776}};
    map<unsigned int, double> germanium_mass = {{70u, 69.924247399999999}, {72u, 71.922075800000002}, {73u, 72.9234589}, {74u, 73.921177799999995}, {76u, 75.921401}};
    buildElement_("Germanium", "Ge", 32u, germanium_abundance, germanium_mass);


    map<unsigned int, double> arsenic_abundance = {{75u, 1.0}};
    map<unsigned int, double> arsenic_mass = {{75u, 74.921596500000007}};
    buildElement_("Arsenic", "As", 33u, arsenic_abundance, arsenic_mass);


    map<unsigned int, double> selenium_abundance = {{74u, 0.00889}, {76u, 0.09366}, {77u, 0.07635}, {78u, 0.23772}, {80u, 0.49607}, {82u, 0.08731}};
    map<unsigned int, double> selenium_mass = {{74u, 73.922476399999994}, {76u, 75.919213600000006}, {77u, 76.919914000000006}, {78u, 77.917309099999997}, {80u, 79.916521299999999}, {82u, 81.916699399999999}};
    buildElement_("Selenium", "Se", 34u, selenium_abundance, selenium_mass);


    map<unsigned int, double> bromine_abundance = {{79u, 0.5069}, {81u, 0.49310000000000004}};
    map<unsigned int, double> bromine_mass = {{79u, 78.918337100000002}, {81u, 80.916290599999996}};
    buildElement_("Bromine", "Br", 35u, bromine_abundance, bromine_mass);


    map<unsigned int, double> krypton_abundance = {{78u, 0.0034999999999999996}, {80u, 0.0225}, {82u, 0.11599999999999999}, {83u, 0.115}, {84u, 0.57}, {86u, 0.17300000000000001}};
    map<unsigned int, double> krypton_mass = {{78u, 77.920400000000001}, {80u, 79.916380000000004}, {82u, 81.913482000000002}, {83u, 82.914135000000002}, {84u, 83.911507}, {86u, 85.910616000000005}};
    buildElement_("Krypton", "Kr", 36u, krypton_abundance, krypton_mass);


    map<unsigned int, double> rubidium_abundance = {{85u, 0.7217}};
    map<unsigned int, double> rubidium_mass = {{85u, 84.911789737999996}};
    buildElement_("Rubidium", "Rb", 37u, rubidium_abundance, rubidium_mass);


    map<unsigned int, double> strontium_abundance = {{84u, 0.005600000000000001}, {86u, 0.0986}, {87u, 0.07}, {88u, 0.8258}};
    map<unsigned int, double> strontium_mass = {{84u, 83.913425000000004}, {86u, 85.909260730900002}, {87u, 86.908877497000006}, {88u, 87.905612257100003}};
    buildElement_("Strontium", "Sr", 38u, strontium_abundance, strontium_mass);


    map<unsigned int, double> yttrium_abundance = {{89u, 1.0}};
    map<unsigned int, double> yttrium_mass = {{89u, 88.905850000000001}};
    buildElement_("Yttrium", "Y", 39u, yttrium_abundance, yttrium_mass);


    map<unsigned int, double> zirconium_abundance = {{90u, 0.5145000000000001}, {91u, 0.11220000000000001}, {92u, 0.17149999999999999}, {94u, 0.17379999999999998}, {96u, 0.0280}};
    map<unsigned int, double> zirconium_mass = {{90u, 89.9047044}, {91u, 90.905645800000002}, {92u, 91.905040799999995}, {94u, 93.906315199999995}, {96u, 95.9082776}};
    buildElement_("Zirconium", "Zr", 40u, zirconium_abundance, zirconium_mass);


    map<unsigned int, double> nibium_abundance = {{93u, 1.0}};
    map<unsigned int, double> nibium_mass = {{93u, 92.906378099999998}};
    buildElement_("Nibium", "Nb", 41u, nibium_abundance, nibium_mass);


    map<unsigned int, double> molybdenum_abundance = {{92u, 0.1484}, {94u, 0.0925}, {95u, 0.1592}, {96u, 0.1668}, {97u, 0.0955}, {98u, 0.2413}, {100u, 0.09630000000000001}};
    map<unsigned int, double> molybdenum_mass = {{92u, 91.906809999999993}, {94u, 93.905088000000006}, {95u, 94.905840999999995}, {96u, 95.904679000000002}, {97u, 96.906020999999996}, {98u, 97.905407999999994}, {100u, 99.907477}};
    buildElement_("Molybdenum", "Mo", 42u, molybdenum_abundance, molybdenum_mass);


    // Technitium(Tc) abundance is not known.


    map<unsigned int, double> ruthenium_abundance = {{96u, 0.0554}, {98u, 0.0187}, {99u, 0.1276}, {100u, 0.126}, {101u, 0.17059999999999997}, {102u, 0.3155}, {104u, 0.1862}};
    map<unsigned int, double> ruthenium_mass = {{96u, 95.907597999999993}, {98u, 97.905287000000001}, {99u, 98.9059393}, {100u, 99.904219499999996}, {101u, 100.905582100000004}, {102u, 101.904349300000007}, {104u, 103.905433000000002}};
    buildElement_("Ruthenium", "Ru", 44u, ruthenium_abundance, ruthenium_mass);


    map<unsigned int, double> rhodium_abundance = {{103u, 1.0}};
    map<unsigned int, double> rhodium_mass = {{103u, 102.905500000000004}};
    buildElement_("Rhodium", "Rh", 45u, rhodium_abundance, rhodium_mass);


    map<unsigned int, double> palladium_abundance = {{102u, 0.0102}, {104u, 0.1114}, {105u, 0.22329999999999997}, {106u, 0.2733}, {108u, 0.2646}, {110u, 0.11720000000000001}};
    map<unsigned int, double> palladium_mass = {{102u, 101.905608999999998}, {104u, 103.904036000000005}, {105u, 104.905085}, {106u, 105.903486000000001}, {108u, 107.903891999999999}, {110u, 109.905152999999999}};
    buildElement_("Palladium", "Pd", 46u, palladium_abundance, palladium_mass);


    map<unsigned int, double> silver_abundance = {{107u, 0.51839}, {109u, 0.48161000000000004}};
    map<unsigned int, double> silver_mass = {{107u, 106.905092999999994}, {109u, 108.904756000000006}};
    buildElement_("Silver", "Ag", 47u, silver_abundance, silver_mass);


    map<unsigned int, double> cadmium_abundance = {{106u, 0.0125}, {108u, 0.0089}, {110u, 0.1249}, {111u, 0.128}, {112u, 0.2413}, {113u, 0.1222}, {114u, 0.2873}, {116u, 0.07490000000000001}};
    map<unsigned int, double> cadmium_mass = {{106u, 105.906458000000001}, {108u, 107.904184000000001}, {110u, 109.903002099999995}, {111u, 110.904178099999996}, {112u, 111.902757800000003}, {113u, 112.904401699999994}, {114u, 113.903358499999996}, {116u, 115.904756000000006}};
    buildElement_("Cadmium", "Cd", 48u, cadmium_abundance, cadmium_mass);


    map<unsigned int, double> indium_abundance = {{113u, 0.0429}, {115u, 0.9571}};
    map<unsigned int, double> indium_mass = {{113u, 112.904060000000001}, {115u, 114.903878000000006}};
    buildElement_("Indium", "In", 49u, indium_abundance, indium_mass);


    map<unsigned int, double> tin_abundance = {{112u, 0.0097}, {114u, 0.0066}, {115u, 0.0034000000000000002}, {116u, 0.1454}, {117u, 0.0768}, {118u, 0.2422}, {119u, 0.0859}, {120u, 0.3258}, {122u, 0.0463}, {124u, 0.0579}};
    map<unsigned int, double> tin_mass = {{112u, 111.904818000000006}, {114u, 113.902777900000004}, {115u, 114.903341999999995}, {116u, 115.901741000000001}, {117u, 116.902951999999999}, {118u, 117.901602999999994}, {119u, 118.903307999999996}, {120u, 119.902194699999996}, {122u, 121.903439000000006}, {124u, 123.905273899999997}};
    buildElement_("Tin", "Sn", 50u, tin_abundance, tin_mass);


    map<unsigned int, double> antimony_abundance = {{121u, 0.5721}, {123u, 0.4279}};
    map<unsigned int, double> antimony_mass = {{121u, 120.903815699999996}, {123u, 122.904213999999996}};
    buildElement_("Antimony", "Sb", 51u, antimony_abundance, antimony_mass);


    map<unsigned int, double> tellurium_abundance = {{120u, 0.0009}, {122u, 0.0255}, {124u, 0.047400000000000005}, {125u, 0.0707}, {126u, 0.1884}, {128u, 0.31739999999999996}, {130u, 0.3408}};
    map<unsigned int, double> tellurium_mass = {{120u, 119.904020000000003}, {122u, 121.9030439}, {124u, 123.902817900000002}, {125u, 124.904430700000006}, {126u, 125.903311700000003}, {128u, 127.904463100000001}, {130u, 129.906224400000014}};
    buildElement_("Tellurium", "Te", 52u, tellurium_abundance, tellurium_mass);


    map<unsigned int, double> iodine_abundance = {{127u, 1.0}};
    map<unsigned int, double> iodine_mass = {{127u, 126.904472999999996}};
    buildElement_("Iodine", "I", 53u, iodine_abundance, iodine_mass);


    map<unsigned int, double> xenon_abundance = {{128u, 0.0191}, {129u, 0.264}, {130u, 0.040999999999999995}, {131u, 0.212}, {132u, 0.26899999999999996}, {134u, 0.10400000000000001}, {136u, 0.08900000000000001}};
    map<unsigned int, double> xenon_mass = {{128u, 127.903531000000001}, {129u, 128.904779999999988}, {130u, 129.903509000000014}, {131u, 130.90507199999999}, {132u, 131.904144000000002}, {134u, 133.905394999999999}, {136u, 135.90721400000001}};
    buildElement_("Xenon", "Xe", 54u, xenon_abundance, xenon_mass);


    map<unsigned int, double> caesium_abundance = {{133u, 1.0}};
    map<unsigned int, double> caesium_mass = {{133u, 132.905451932999995}};
    buildElement_("Caesium", "Cs", 55u, caesium_abundance, caesium_mass);


    map<unsigned int, double> barium_abundance = {{132u, 0.00101}, {134u, 0.024169999999999997}, {135u, 0.06591999999999999}, {136u, 0.07854}, {137u, 0.11231999999999999}, {138u, 0.71698}};
    map<unsigned int, double> barium_mass = {{132u, 131.9050613}, {134u, 133.904508399999997}, {135u, 134.905688599999991}, {136u, 135.904575899999998}, {137u, 136.905827399999993}, {138u, 137.905247199999991}};
    buildElement_("Barium", "Ba", 56u, barium_abundance, barium_mass);


    map<unsigned int, double> lanthanum_abundance = {{138u, 0.00089}, {139u, 0.99911}};
    map<unsigned int, double> lanthanum_mass = {{138u, 137.907112000000012}, {139u, 138.906353300000006}};
    buildElement_("Lanthanum", "La", 57u, lanthanum_abundance, lanthanum_mass);


    map<unsigned int, double> cerium_abundance = {{136u, 0.00185}, {138u, 0.00251}, {140u, 0.8845000000000001}, {142u, 0.11114}};
    map<unsigned int, double> cerium_mass = {{136u, 135.907172000000003}, {138u, 137.905991}, {140u, 139.905438699999991}, {142u, 141.909244000000001}};
    buildElement_("Cerium", "Ce", 58u, cerium_abundance, cerium_mass);


    map<unsigned int, double> praseodymium_abundance = {{141u, 1.0}};
    map<unsigned int, double> praseodymium_mass = {{141u, 140.907646999999997}};
    buildElement_("Praseodymium", "Pr", 59u, praseodymium_abundance, praseodymium_mass);


    map<unsigned int, double> neodymium_abundance = {{142u, 0.272}, {143u, 0.122}, {144u, 0.23800000000000002}, {145u, 0.083}, {146u, 0.172}, {148u, 0.057999999999999996}, {150u, 0.055999999999999994}};
    map<unsigned int, double> neodymium_mass = {{142u, 141.907723299999987}, {143u, 142.909814299999994}, {144u, 143.910087299999987}, {145u, 144.912573600000002}, {146u, 145.913116900000006}, {148u, 147.916892999999988}, {150u, 149.920891000000012}};
    buildElement_("Neodymium", "Nd", 60u, neodymium_abundance, neodymium_mass);


    // Promethium(Pm) abundance is not known.


    map<unsigned int, double> samarium_abundance = {{144u, 0.0308}, {147u, 0.15}, {148u, 0.1125}, {149u, 0.1382}, {150u, 0.0737}, {152u, 0.26739999999999997}, {154u, 0.2274}};
    map<unsigned int, double> samarium_mass = {{144u, 143.911999000000009}, {147u, 146.9148979}, {148u, 147.914822700000002}, {149u, 148.917184700000007}, {150u, 149.917275499999988}, {152u, 151.919732399999987}, {154u, 153.92220929999999}};
    buildElement_("Samarium", "Sm", 62u, samarium_abundance, samarium_mass);


    map<unsigned int, double> europium_abundance = {{151u, 0.4781}, {153u, 0.5219}};
    map<unsigned int, double> europium_mass = {{151u, 150.919857}, {153u, 152.921237}};
    buildElement_("Europium", "Eu", 63u, europium_abundance, europium_mass);


    map<unsigned int, double> gadolinium_abundance = {{152u, 0.002}, {154u, 0.0218}, {155u, 0.14800000000000002}, {156u, 0.2047}, {157u, 0.1565}, {158u, 0.2484}, {160u, 0.2186}};
    map<unsigned int, double> gadolinium_mass = {{152u, 151.919791000000004}, {154u, 153.920865600000013}, {155u, 154.92262199999999}, {156u, 155.922122699999989}, {157u, 156.923960099999988}, {158u, 157.924103900000006}, {160u, 159.927054099999992}};
    buildElement_("Gadolinium", "Gd", 64u, gadolinium_abundance, gadolinium_mass);


    map<unsigned int, double> terbium_abundance = {{159u, 1.0}};
    map<unsigned int, double> terbium_mass = {{159u, 158.925354}};
    buildElement_("Terbium", "Tb", 65u, terbium_abundance, terbium_mass);


    map<unsigned int, double> dysprosium_abundance = {{156u, 0.00056}, {158u, 0.00095}, {160u, 0.02329}, {161u, 0.18889}, {162u, 0.25475}, {163u, 0.24896}, {164u, 0.28260}};
    map<unsigned int, double> dysprosium_mass = {{156u, 155.924284}, {158u, 157.92441}, {160u,  159.925203}, {161u, 160.926939}, {162u, 161.926804}, {163u, 162.928737}, {164u,  163.929181}};
    buildElement_("Dysprosium", "Dy", 66u, dysprosium_abundance, dysprosium_mass);


    map<unsigned int, double> holmium_abundance = {{165u, 1.0}};
    map<unsigned int, double> holmium_mass = {{165u,  164.930328}};
    buildElement_("Holmium", "Ho", 67u, holmium_abundance, holmium_mass);


    map<unsigned int, double> erbium_abundance = {{162u, 0.00056}, {164u, 0.01601}, {166u, 0.33503}, {167u, 0.22869}, {168u, 0.26978}, {170u, 0.14910}};
    map<unsigned int, double> erbium_mass = {{162u, 161.928787}, {164u, 163.929207}, {166u, 165.930299}, {167u, 166.932054}, {168u, 167.932376}, {170u, 169.93547}};
    buildElement_("Erbium", "Er", 68u, erbium_abundance, erbium_mass);


    map<unsigned int, double> thulium_abundance = {{169u, 1.0}};
    map<unsigned int, double> thulium_mass = {{169u,  168.934218}};
    buildElement_("Thulium", "Tm", 69u, thulium_abundance, thulium_mass);


    map<unsigned int, double> ytterbium_abundance = {{168u, 0.00126}, {170u, 0.03023}, {171u, 0.14216}, {172u, 0.21754}, {173u, 0.16098}, {174u, 0.31896}, {176u, 0.12887}};
    map<unsigned int, double> ytterbium_mass = {{168u,  167.933889}, {170u, 169.93476725}, {171u, 170.93633152}, {172u, 171.93638666}, {173u, 172.93821622}, {174u, 173.93886755}, {176u, 175.9425747}};
    buildElement_("Ytterbium", "Yb", 70u, ytterbium_abundance, ytterbium_mass);


    map<unsigned int, double> lutetium_abundance = {{175u,  0.97401}, {176u, 0.02599}};
    map<unsigned int, double> lutetium_mass = {{175u, 174.940777}, {176u, 175.942692}};
    buildElement_("Lutetium", "Lu", 71u, lutetium_abundance, lutetium_mass);


    map<unsigned int, double> hafnium_abundance = {{176u, 0.0526}, {177u, 0.18600000000000003}, {178u, 0.2728}, {179u, 0.1362}, {180u, 0.3508}};
    map<unsigned int, double> hafnium_mass = {{176u, 175.941408599999988}, {177u, 176.943220700000012}, {178u, 177.943698799999993}, {179u, 178.945816100000002}, {180u, 179.946550000000002}};
    buildElement_("Hafnium", "Hf", 72u, hafnium_abundance, hafnium_mass);


    map<unsigned int, double> tantalum_abundance = {{180u, 0.0001176}, {181u, 0.99988}};
    map<unsigned int, double> tantalum_mass = {{180u, 179.94747}, {181u, 180.947995800000001}};
    buildElement_("Tantalum", "Ta", 73u, tantalum_abundance, tantalum_mass);


    map<unsigned int, double> tungsten_abundance = {{180u, 0.0012}, {182u, 0.265}, {183u, 0.1431}, {184u, 0.3064}, {186u, 0.2843}};
    map<unsigned int, double> tungsten_mass = {{180u, 179.946704000000011}, {182u, 181.948204199999992}, {183u, 182.950222999999994}, {184u, 183.950930999999997}, {186u, 185.954364099999992}};
    buildElement_("Tungsten", "W", 74u, tungsten_abundance, tungsten_mass);


    map<unsigned int, double> rhenium_abundance = {{185u, 0.374}, {187u, 0.626}};
    map<unsigned int, double> rhenium_mass = {{185u, 184.952955000000003}, {187u, 186.95575310000001}};
    buildElement_("Rhenium", "Re", 75u, rhenium_abundance, rhenium_mass);


    map<unsigned int, double> osmium_abundance = {{184u, 0.0002}, {186u, 0.0159}, {187u, 0.0196}, {188u, 0.1324}, {189u, 0.1615}, {190u, 0.2626}, {192u, 0.4078}};
    map<unsigned int, double> osmium_mass = {{184u, 183.952493}, {186u, 185.953838}, {187u, 186.955750}, {188u, 187.955837}, {189u, 188.958146}, {190u, 189.958446}, {192u, 191.96148}};
    buildElement_("Osmium", "Os", 76u, osmium_abundance, osmium_mass);


    map<unsigned int, double> iridium_abundance = {{191u, 0.3723}, {193u, 0.6277}};
    map<unsigned int, double> iridium_mass = {{191u, 190.960591}, {193u, 192.962924}};
    buildElement_("Iridium", "Ir", 77u, rhenium_abundance, rhenium_mass);


    // Pt-190 is radioactive but with a very long half-life. Since its natural occurence is very low, we neglect it by default (m=189.959930 abund.frac.=0.00014)
    // TODO re-evaluate inclusion?
    map<unsigned int, double> platinum_abundance = {{192u, 0.00782}, {194u, 0.32966999999999996}, {195u, 0.33832}, {196u, 0.25242000000000003}, {198u, 0.07163}};
    map<unsigned int, double> platinum_mass = {{192u, 191.961038000000002}, {194u, 193.962680299999988}, {195u, 194.964791100000014}, {196u, 195.964951500000012}, {198u, 197.967893000000004}};
    buildElement_("Platinum", "Pt", 78u, platinum_abundance, platinum_mass);


    map<unsigned int, double> gold_abundance = {{197u, 1.0}};
    map<unsigned int, double> gold_mass = {{197u, 196.96655100000001}};
    buildElement_("Gold", "Au", 79u, gold_abundance, gold_mass);


    map<unsigned int, double> mercury_abundance = {{196u, 0.0015}, {198u, 0.09970000000000001}, {199u, 0.16870000000000002}, {200u, 0.231}, {201u, 0.1318}, {202u, 0.2986}, {204u, 0.0687}};
    map<unsigned int, double> mercury_mass = {{196u, 195.965833000000004}, {198u, 197.966768999999999}, {199u, 198.968279899999999}, {200u, 199.968325999999991}, {201u, 200.970302299999986}, {202u, 201.970642999999996}, {204u, 203.973493899999994}};
    buildElement_("Mercury", "Hg", 80u, mercury_abundance, mercury_mass);


    map<unsigned int, double> thallium_abundance = {{203u, 0.2952}, {205u, 0.7048000000000001}};
    map<unsigned int, double> thallium_mass = {{203u, 202.972344200000009}, {205u, 204.97442749999999}};
    buildElement_("Thallium", "Tl", 81u, thallium_abundance, thallium_mass);


    map<unsigned int, double> lead_abundance = {{204u, 0.013999999999999999}, {206u, 0.24100000000000002}, {207u, 0.221}, {208u, 0.524}};
    map<unsigned int, double> lead_mass = {{204u, 203.973043600000011}, {206u, 205.974465299999991}, {207u, 206.975896900000009}, {208u, 207.976653800000008}};
    buildElement_("Lead", "Pb", 82u, lead_abundance, lead_mass);


    map<unsigned int, double> bismuth_abundance = {{209u, 1.0}};
    map<unsigned int, double> bismuth_mass = {{209u, 208.980398699999995}};
    buildElement_("Bismuth", "Bi", 83u, bismuth_abundance, bismuth_mass);


    // Polonium (Pb) abundance is not known.

    // Astatine(At) abundance is not known.

    // Radon(Rn) abundance is not known.

    // Radium(Ra) abundance is not known.

    map<unsigned int, double> thorium_abundance = {{230u, 0.0002}, {232u, 0.9998}};
    map<unsigned int, double> thorium_mass = {{230u, 230.033133800000002}, {232u, 232.038055299999996}};
    buildElement_("Thorium", "Th", 90u, thorium_abundance, thorium_mass);


    map<unsigned int, double> protactinium_abundance = {{231u, 1.0}};
    map<unsigned int, double> protactinium_mass = {{231u, 231.03588}};
    buildElement_("Protactinium", "Pa", 91u, protactinium_abundance, protactinium_mass);


    map<unsigned int, double> uranium_abundance = {{234u,  0.000054}, {235u, 0.007204}, {238u, 0.992742}};
    map<unsigned int, double> uranium_mass = {{234u,  234.040950}, {235u,  235.043928}, {238u,   238.05079}};
    buildElement_("Uranium", "U", 92u, uranium_abundance, uranium_mass);

    // special case for deuterium and tritium: add symbol alias
    const Element* deuterium = getElement("(2)H");
    symbols_["D"] = deuterium;
    const Element* tritium = getElement("(3)H");
    symbols_["T"] = tritium;

    // Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts and Og Abundances are not known.

  }


  void ElementDB::buildElement_(const string& name, const string& symbol, const unsigned int an, const map<unsigned int, double>& abundance, const map<unsigned int, double>& mass)
  {
    IsotopeDistribution isotopes = parseIsotopeDistribution_(abundance, mass);
    double avg_weight = calculateAvgWeight_(abundance, mass);
    double mono_weight = calculateMonoWeight_(abundance, mass);

    addElementToMaps_(name, symbol, an, make_unique<const Element>(name, symbol, an, avg_weight, mono_weight, isotopes));
    storeIsotopes_(name, symbol, an, mass, isotopes);
  }

  void overwrite(const Element* old, unique_ptr<const Element>& new_e)
  {
    if (old->getSymbol() != new_e->getSymbol())
    { // -- this would invalidate the lookup, since e_ptr->getSymbols().at("O")->getSymbol() == 'P'
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, new_e->getSymbol(),
                                    "Replacing element with name " + old->getName() + " and symbol " + old->getSymbol() + " has different new symbol: " + new_e->getSymbol());
    }
    if (old->getName() != new_e->getName())
    { // -- this would invalidate the lookup, since e_ptr->getName().at("Oxygen")->getName() == 'Something'
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, new_e->getSymbol(), "Replacing element with name " + old->getName() + " has different new name: " + new_e->getName());
    }
    if (old->getAtomicNumber() != new_e->getAtomicNumber())
    { // -- this would invalidate the lookup, since e_ptr->getAtomicNumbers().at(12)->getAtomicNumber() == 14
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, new_e->getSymbol(),
                                    "Replacing element with atomic number " + String(old->getAtomicNumber()) + " has different new atomic number: " + String(new_e->getAtomicNumber()));
    }
    // ... overwrite
    *(const_cast<Element*>(old)) = *new_e;
  }

  void ElementDB::addElementToMaps_(const string& name, const string& symbol, const unsigned int an, unique_ptr<const Element> e)
  {
    // overwrite existing element if it already exists
    // find() has to be protected here in a parallel context
    if (atomic_numbers_.find(an) != atomic_numbers_.end())
    {
      // in order to ensure that existing elements are still valid and memory
      // addresses do not change, we have to modify the Element in place
      // instead of replacing it.
      overwrite(atomic_numbers_[an], e);
      // do not release 'e' here; it needs to be deleted when it goes out if scope
    }
    else
    {
      addIfUniqueOrThrow(names_, name, e);
      addIfUniqueOrThrow(symbols_, symbol, e);
      addIfUniqueOrThrow(atomic_numbers_, an, e);
      e.release(); // allocation will be cleaned up by ~ElementDB now
    }
  }

  void ElementDB::storeIsotopes_(const string& name, const string& symbol, const unsigned int an, const map<unsigned int, double>& mass, const IsotopeDistribution& isotopes)
  {
    for (const auto& isotope : isotopes)
    {
      double atomic_mass = isotope.getMZ();
      unsigned int mass_number = std::round(atomic_mass);
      string iso_name = "(" + std::to_string(mass_number) + ")" + name;
      string iso_symbol = "(" + std::to_string(mass_number) + ")" + symbol;

      // set avg and mono to same value for isotopes (old hack...)
      double iso_avg_weight = mass.at(mass_number);
      double iso_mono_weight = iso_avg_weight;
      IsotopeDistribution iso_isotopes;
      IsotopeDistribution::ContainerType iso_container;
      iso_container.push_back(Peak1D(atomic_mass, 1.0));
      iso_isotopes.set(iso_container);  

      auto iso_element = make_unique<const Element>(iso_name, iso_symbol, an, iso_avg_weight, iso_mono_weight, iso_isotopes);
      if (auto has_elem = names_.find(iso_name); has_elem != names_.end())
      { // already exists: overwrite (affects all maps, since they all point to the same thing)
        overwrite(has_elem->second, iso_element);
        // do not release 'iso_element' here; it needs to be deleted when it goes out if scope
      }
      else
      {
        addIfUniqueOrThrow(names_, iso_name, iso_element);
        addIfUniqueOrThrow(symbols_, iso_symbol, iso_element);
        iso_element.release(); // allocation will be cleaned up by ~ElementDB now
      }
    } 
  }

  IsotopeDistribution ElementDB::parseIsotopeDistribution_(const map<unsigned int, double>& abundance, const map<unsigned int, double>& mass)
  {
    IsotopeDistribution::ContainerType dist;
    
    for (map<unsigned int, double>::const_iterator it = abundance.begin(); it != abundance.end(); ++it)
    { 
      dist.push_back(Peak1D(mass.at(it->first) , abundance.at(it->first)));
    }

    IsotopeDistribution iso_dist;
    iso_dist.set(dist);
    
    return iso_dist;
  }

  void ElementDB::clear_()
  {
    // names_ has the union of all Element*, deleting this is sufficient to avoid mem leaks
    for (auto it = names_.begin(); it != names_.end(); ++it)
    {
      delete it->second;
    }
    names_.clear();
    symbols_.clear();
    atomic_numbers_.clear();
  }

} // namespace OpenMS
