// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch, Timo Sachsenberg, Chris Bielow, Jang Jang Jin‚$
// --------------------------------------------------------------------------
//
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/SYSTEM/File.h>

#include <iostream>

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

  const ElementDB* ElementDB::getInstance()
  {
    static ElementDB* db_ = new ElementDB;
    return db_;
  }

  const map<string, const Element*>& ElementDB::getNames() const
  {
    return names_;
  }

  const map<string, const Element*>& ElementDB::getSymbols() const
  {
    return symbols_;
  }

  const map<unsigned int, const Element*>& ElementDB::getAtomicNumbers() const
  {
    return atomic_numbers_;
  }

  const Element* ElementDB::getElement(const string& name) const
  {
    if (auto entry = names_.find(name); entry != names_.end())
    {
      return entry->second;
    }
    else
    {
      if (auto entry = symbols_.find(name); entry != symbols_.end())
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
    return atomic_numbers_.count(atomic_number) == 1;
  }

  double ElementDB::calculateAvgWeight_(const map<unsigned int, double>& abundance, const map<unsigned int, double>& mass)
  {
    double avg = 0;
    // calculate weighted average
    for (map<unsigned int, double>::const_iterator it = abundance.begin(); it != abundance.end(); ++it)
    {
      avg += mass.at(it->first) * abundance.at(it->first);
    }
    return avg;
  }

  double ElementDB::calculateMonoWeight_(const map<unsigned int, double>& mass)
  {
    double smallest_weight = 1e10;

    for (map<unsigned int, double>::const_iterator it = mass.begin(); it != mass.end(); ++it)
    {
      if (it->second < smallest_weight)
      {
        smallest_weight = it->second;
      }
    }

    return smallest_weight;
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
    buildElement_("Bor", "B", 5u, bor_abundance, bor_mass);


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


    map<unsigned int, double> argon_abundance = {{36u, 0.003336}, {38u, 0.000629}, {40u, 0.996035}};
    map<unsigned int, double> argon_mass = {{36u, 35.967545106000003}, {38u, 37.9627324}, {40u, 39.9623831225}};
    buildElement_("Argon", "Ar", 18u, argon_abundance, argon_mass);


    map<unsigned int, double> titanium_abundance = {{46u, 0.0825}, {47u, 0.07440000000000001}, {48u, 0.7372}, {49u, 0.0541}, {50u, 0.0518}};
    map<unsigned int, double> titanium_mass = {{46u, 45.952631599999997}, {47u, 46.951763100000001}, {48u, 47.947946299999998}, {49u, 48.947870000000002}, {50u, 49.944791199999997}};
    buildElement_("Titanium", "Ti", 22u, titanium_abundance, titanium_mass);


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


    map<unsigned int, double> potassium_abundance = {{39u, 0.932581}, {40u, 0.000117}, {41u, 0.067302}};
    map<unsigned int, double> potassium_mass = {{39u, 38.963706680000001}, {40u, 39.963998480000001}, {41u, 40.961825760000004}};
    buildElement_("Potassium", "K", 19u, potassium_abundance, potassium_mass);


    map<unsigned int, double> calcium_abundance = {{40u, 0.96941}, {42u, 0.00647}, {43u, 0.00135}, {44u, 0.02086}, {46u, 4e-05}, {48u, 0.00187}};
    map<unsigned int, double> calcium_mass = {{40u, 39.962590980000002}, {42u, 41.958618010000002}, {43u, 42.958766599999997}, {44u, 43.955481800000001}, {46u, 45.953692599999997}, {48u, 47.952534}};
    buildElement_("Calcium", "Ca", 20u, calcium_abundance, calcium_mass);


    map<unsigned int, double> scandium_abundance = {{45u, 1.0}};
    map<unsigned int, double> scandium_mass = {{45u, 44.955910000000003}};
    buildElement_("Scandium", "Sc", 21u, scandium_abundance, scandium_mass);


    map<unsigned int, double> vanadium_abundance = {{50u, 0.0025}, {51u, 0.9975}};
    map<unsigned int, double> vanadium_mass = {{50u, 49.947158500000001}, {51u, 50.943959499999998}};
    buildElement_("Vanadium", "V", 23u, vanadium_abundance, vanadium_mass);


    map<unsigned int, double> chromium_abundance = {{50u, 0.043449999999999996}, {52u, 0.83789}, {53u, 0.09501}, {54u, 0.02365}};
    map<unsigned int, double> chromium_mass = {{50u, 49.946044200000003}, {52u, 51.940507500000003}, {53u, 52.940649399999998}, {54u, 53.938880400000002}};
    buildElement_("Chromium", "Cr", 24u, chromium_abundance, chromium_mass);


    map<unsigned int, double> tellurium_abundance = {{120u, 0.0009}, {122u, 0.0255}, {124u, 0.047400000000000005}, {125u, 0.0707}, {126u, 0.1884}, {128u, 0.31739999999999996}, {130u, 0.3408}};
    map<unsigned int, double> tellurium_mass = {{120u, 119.904020000000003}, {122u, 121.9030439}, {124u, 123.902817900000002}, {125u, 124.904430700000006}, {126u, 125.903311700000003}, {128u, 127.904463100000001}, {130u, 129.906224400000014}};
    buildElement_("Tellurium", "Te", 52u, tellurium_abundance, tellurium_mass);


    map<unsigned int, double> barium_abundance = {{132u, 0.00101}, {134u, 0.024169999999999997}, {135u, 0.06591999999999999}, {136u, 0.07854}, {137u, 0.11231999999999999}, {138u, 0.71698}};
    map<unsigned int, double> barium_mass = {{132u, 131.9050613}, {134u, 133.904508399999997}, {135u, 134.905688599999991}, {136u, 135.904575899999998}, {137u, 136.905827399999993}, {138u, 137.905247199999991}};
    buildElement_("Barium", "Ba", 56u, barium_abundance, barium_mass);


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


    map<unsigned int, double> rubidium_abundance = {{85u, 0.7217}};
    map<unsigned int, double> rubidium_mass = {{85u, 84.911789737999996}};
    buildElement_("Rubidium", "Rb", 37u, rubidium_abundance, rubidium_mass);


    map<unsigned int, double> strontium_abundance = {{84u, 0.005600000000000001}, {86u, 0.0986}, {87u, 0.07}, {88u, 0.8258}};
    map<unsigned int, double> strontium_mass = {{84u, 83.913425000000004}, {86u, 85.909260730900002}, {87u, 86.908877497000006}, {88u, 87.905612257100003}};
    buildElement_("Strontium", "Sr", 38u, strontium_abundance, strontium_mass);


    map<unsigned int, double> yttrium_abundance = {{89u, 1.0}};
    map<unsigned int, double> yttrium_mass = {{89u, 88.905850000000001}};
    buildElement_("Yttrium", "Y", 39u, yttrium_abundance, yttrium_mass);


    map<unsigned int, double> zirconium_abundance = {{90u, 0.5145000000000001}, {91u, 0.11220000000000001}, {92u, 0.17149999999999999}, {94u, 0.17379999999999998}};
    map<unsigned int, double> zirconium_mass = {{90u, 89.9047044}, {91u, 90.905645800000002}, {92u, 91.905040799999995}, {94u, 93.906315199999995}};
    buildElement_("Zirconium", "Zr", 40u, zirconium_abundance, zirconium_mass);


    map<unsigned int, double> nibium_abundance = {{93u, 1.0}};
    map<unsigned int, double> nibium_mass = {{93u, 92.906378099999998}};
    buildElement_("Nibium", "Nb", 41u, nibium_abundance, nibium_mass);


    map<unsigned int, double> ruthenium_abundance = {{96u, 0.0554}, {98u, 0.0187}, {99u, 0.1276}, {100u, 0.126}, {101u, 0.17059999999999997}, {102u, 0.3155}, {104u, 0.1862}};
    map<unsigned int, double> ruthenium_mass = {{96u, 95.907597999999993}, {98u, 97.905287000000001}, {99u, 98.9059393}, {100u, 99.904219499999996}, {101u, 100.905582100000004}, {102u, 101.904349300000007}, {104u, 103.905433000000002}};
    buildElement_("Ruthenium", "Ru", 44u, ruthenium_abundance, ruthenium_mass);


    map<unsigned int, double> tin_abundance = {{112u, 0.0097}, {114u, 0.0066}, {115u, 0.0034000000000000002}, {116u, 0.1454}, {117u, 0.0768}, {118u, 0.2422}, {119u, 0.0859}, {120u, 0.3258}, {122u, 0.0463}, {124u, 0.0579}};
    map<unsigned int, double> tin_mass = {{112u, 111.904818000000006}, {114u, 113.902777900000004}, {115u, 114.903341999999995}, {116u, 115.901741000000001}, {117u, 116.902951999999999}, {118u, 117.901602999999994}, {119u, 118.903307999999996}, {120u, 119.902194699999996}, {122u, 121.903439000000006}, {124u, 123.905273899999997}};
    buildElement_("Tin", "Sn", 50u, tin_abundance, tin_mass);


    map<unsigned int, double> antimony_abundance = {{121u, 0.5721}, {123u, 0.4279}};
    map<unsigned int, double> antimony_mass = {{121u, 120.903815699999996}, {123u, 122.904213999999996}};
    buildElement_("Antimony", "Sb", 51u, antimony_abundance, antimony_mass);


    map<unsigned int, double> selenium_abundance = {{74u, 0.00889}, {76u, 0.09366}, {77u, 0.07635}, {78u, 0.23772}, {80u, 0.49607}, {82u, 0.08731}};
    map<unsigned int, double> selenium_mass = {{74u, 73.922476399999994}, {76u, 75.919213600000006}, {77u, 76.919914000000006}, {78u, 77.917309099999997}, {80u, 79.916521299999999}, {82u, 81.916699399999999}};
    buildElement_("Selenium", "Se", 34u, selenium_abundance, selenium_mass);


    map<unsigned int, double> bromine_abundance = {{79u, 0.5069}, {81u, 0.49310000000000004}};
    map<unsigned int, double> bromine_mass = {{79u, 78.918337100000002}, {81u, 80.916290599999996}};
    buildElement_("Bromine", "Br", 35u, bromine_abundance, bromine_mass);


    map<unsigned int, double> krypton_abundance = {{78u, 0.0034999999999999996}, {80u, 0.0225}, {82u, 0.11599999999999999}, {83u, 0.115}, {84u, 0.57}, {86u, 0.17300000000000001}};
    map<unsigned int, double> krypton_mass = {{78u, 77.920400000000001}, {80u, 79.916380000000004}, {82u, 81.913482000000002}, {83u, 82.914135000000002}, {84u, 83.911507}, {86u, 85.910616000000005}};
    buildElement_("Krypton", "Kr", 36u, krypton_abundance, krypton_mass);


    map<unsigned int, double> molybdenum_abundance = {{92u, 0.1484}, {94u, 0.0925}, {95u, 0.1592}, {96u, 0.1668}, {97u, 0.0955}, {98u, 0.2413}, {100u, 0.09630000000000001}};
    map<unsigned int, double> molybdenum_mass = {{92u, 91.906809999999993}, {94u, 93.905088000000006}, {95u, 94.905840999999995}, {96u, 95.904679000000002}, {97u, 96.906020999999996}, {98u, 97.905407999999994}, {100u, 99.907477}};
    buildElement_("Molybdenum", "Mo", 42u, molybdenum_abundance, molybdenum_mass);


    map<unsigned int, double> technitium_abundance = {{97u, 0.0}, {98u, 0.0}, {99u, 0.0}};
    map<unsigned int, double> technitium_mass = {{97u, 96.906363999999996}, {98u, 97.907214999999994}, {99u, 98.906254000000004}};
    buildElement_("Technitium", "Tc", 43u, technitium_abundance, technitium_mass);


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


    map<unsigned int, double> iodine_abundance = {{127u, 1.0}};
    map<unsigned int, double> iodine_mass = {{127u, 126.904472999999996}};
    buildElement_("Iodine", "I", 53u, iodine_abundance, iodine_mass);


    map<unsigned int, double> xenon_abundance = {{128u, 0.0191}, {129u, 0.264}, {130u, 0.040999999999999995}, {131u, 0.212}, {132u, 0.26899999999999996}, {134u, 0.10400000000000001}, {136u, 0.08900000000000001}};
    map<unsigned int, double> xenon_mass = {{128u, 127.903531000000001}, {129u, 128.904779999999988}, {130u, 129.903509000000014}, {131u, 130.90507199999999}, {132u, 131.904144000000002}, {134u, 133.905394999999999}, {136u, 135.90721400000001}};
    buildElement_("Xenon", "Xe", 54u, xenon_abundance, xenon_mass);


    map<unsigned int, double> caesium_abundance = {{133u, 1.0}};
    map<unsigned int, double> caesium_mass = {{133u, 132.905451932999995}};
    buildElement_("Caesium", "Cs", 55u, caesium_abundance, caesium_mass);


    map<unsigned int, double> cerium_abundance = {{136u, 0.00185}, {138u, 0.00251}, {140u, 0.8845000000000001}, {142u, 0.11114}};
    map<unsigned int, double> cerium_mass = {{136u, 135.907172000000003}, {138u, 137.905991}, {140u, 139.905438699999991}, {142u, 141.909244000000001}};
    buildElement_("Cerium", "Ce", 58u, cerium_abundance, cerium_mass);


    map<unsigned int, double> praseodymium_abundance = {{141u, 1.0}};
    map<unsigned int, double> praseodymium_mass = {{141u, 140.907646999999997}};
    buildElement_("Praseodymium", "Pr", 59u, praseodymium_abundance, praseodymium_mass);


    map<unsigned int, double> gadolinium_abundance = {{152u, 0.002}, {154u, 0.0218}, {155u, 0.14800000000000002}, {156u, 0.2047}, {157u, 0.1565}, {158u, 0.2484}, {160u, 0.2186}};
    map<unsigned int, double> gadolinium_mass = {{152u, 151.919791000000004}, {154u, 153.920865600000013}, {155u, 154.92262199999999}, {156u, 155.922122699999989}, {157u, 156.923960099999988}, {158u, 157.924103900000006}, {160u, 159.927054099999992}};
    buildElement_("Gadolinium", "Gd", 64u, gadolinium_abundance, gadolinium_mass);


    map<unsigned int, double> hafnium_abundance = {{176u, 0.0526}, {177u, 0.18600000000000003}, {178u, 0.2728}, {179u, 0.1362}, {180u, 0.3508}};
    map<unsigned int, double> hafnium_mass = {{176u, 175.941408599999988}, {177u, 176.943220700000012}, {178u, 177.943698799999993}, {179u, 178.945816100000002}, {180u, 179.946550000000002}};
    buildElement_("Hafnium", "Hf", 72u, hafnium_abundance, hafnium_mass);


    map<unsigned int, double> tantalum_abundance = {{181u, 1.0}};
    map<unsigned int, double> tantalum_mass = {{181u, 180.947995800000001}};
    buildElement_("Tantalum", "Ta", 73u, tantalum_abundance, tantalum_mass);


    map<unsigned int, double> platinum_abundance = {{192u, 0.00782}, {194u, 0.32966999999999996}, {195u, 0.33832}, {196u, 0.25242000000000003}, {198u, 0.07163}};
    map<unsigned int, double> platinum_mass = {{192u, 191.961038000000002}, {194u, 193.962680299999988}, {195u, 194.964791100000014}, {196u, 195.964951500000012}, {198u, 197.967893000000004}};
    buildElement_("Platinum", "Pt", 78u, platinum_abundance, platinum_mass);


    map<unsigned int, double> tungsten_abundance = {{180u, 0.0012}, {182u, 0.265}, {183u, 0.1431}, {184u, 0.3064}, {186u, 0.2843}};
    map<unsigned int, double> tungsten_mass = {{180u, 179.946704000000011}, {182u, 181.948204199999992}, {183u, 182.950222999999994}, {184u, 183.950930999999997}, {186u, 185.954364099999992}};
    buildElement_("Tungsten", "W", 74u, tungsten_abundance, tungsten_mass);


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


    map<unsigned int, double> rhenium_abundance = {{185u, 0.374}, {187u, 0.626}};
    map<unsigned int, double> rhenium_mass = {{185u, 184.952955000000003}, {187u, 186.95575310000001}};
    buildElement_("Rhenium", "Re", 75u, rhenium_abundance, rhenium_mass);


    map<unsigned int, double> neodymium_abundance = {{142u, 0.272}, {143u, 0.122}, {144u, 0.23800000000000002}, {145u, 0.083}, {146u, 0.172}, {148u, 0.057999999999999996}, {150u, 0.055999999999999994}};
    map<unsigned int, double> neodymium_mass = {{142u, 141.907723299999987}, {143u, 142.909814299999994}, {144u, 143.910087299999987}, {145u, 144.912573600000002}, {146u, 145.913116900000006}, {148u, 147.916892999999988}, {150u, 149.920891000000012}};
    buildElement_("Neodymium", "Nd", 60u, neodymium_abundance, neodymium_mass);


    map<unsigned int, double> thorium_abundance = {{230u, 0.0002}, {232u, 0.9998}};
    map<unsigned int, double> thorium_mass = {{230u, 230.033133800000002}, {232u, 232.038055299999996}};
    buildElement_("Thorium", "Th", 90u, thorium_abundance, thorium_mass);


    map<unsigned int, double> lanthanum_abundance = {{138u, 0.00089}, {139u, 0.99911}};
    map<unsigned int, double> lanthanum_mass = {{138u, 137.907112000000012}, {139u, 138.906353300000006}};
    buildElement_("Lanthanum", "La", 57u, lanthanum_abundance, lanthanum_mass);


    map<unsigned int, double> samarium_abundance = {{144u, 0.0308}, {147u, 0.15}, {148u, 0.1125}, {149u, 0.1382}, {150u, 0.0737}, {152u, 0.26739999999999997}, {154u, 0.2274}};
    map<unsigned int, double> samarium_mass = {{144u, 143.911999000000009}, {147u, 146.9148979}, {148u, 147.914822700000002}, {149u, 148.917184700000007}, {150u, 149.917275499999988}, {152u, 151.919732399999987}, {154u, 153.92220929999999}};
    buildElement_("Samarium", "Sm", 62u, samarium_abundance, samarium_mass);
  }

  void ElementDB::buildElement_(const string& name, const string& symbol, const unsigned int an, const map<unsigned int, double>& abundance, const map<unsigned int, double>& mass)
  {
    IsotopeDistribution isotopes = parseIsotopeDistribution_(abundance, mass);
    double avg_weight = calculateAvgWeight_(abundance, mass);
    double mono_weight = calculateMonoWeight_(mass);

    Element* e = new Element(name, symbol, an, avg_weight, mono_weight, isotopes);
    addElementToMaps_(name, symbol, an, e);
    storeIsotopes_(name, symbol, an, mass, isotopes);
  }

  void ElementDB::addElementToMaps_(const string& name, const string& symbol, const unsigned int an, const Element* e)
  {
    names_[name] = e;
    symbols_[symbol] = e;
    atomic_numbers_[an] = e;
  }

  void ElementDB::storeIsotopes_(const string& name, const string& symbol, const unsigned int an, const map<unsigned int, double>& mass, const IsotopeDistribution& isotopes)
  {
    for (const auto& isotope : isotopes)
    {
      double atomic_mass = isotope.getMZ();
      unsigned int mass_number = round(atomic_mass);
      string iso_name = "(" + std::to_string(mass_number) + ")" + name;
      string iso_symbol = "(" + std::to_string(mass_number) + ")" + symbol;

      // set avg and mono to same value for isotopes (old hack...)
      double iso_avg_weight = mass.at(mass_number);
      double iso_mono_weight = iso_avg_weight;
      IsotopeDistribution iso_isotopes;
      IsotopeDistribution::ContainerType iso_container;
      iso_container.push_back(Peak1D(atomic_mass, 1.0));
      iso_isotopes.set(iso_container);  

      Element* iso_e = new Element(iso_name, iso_symbol, an, iso_avg_weight, iso_mono_weight, iso_isotopes);
      names_[iso_name] = iso_e;
      symbols_[iso_symbol] = iso_e;
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
    map<string, const Element*>::iterator it = names_.begin();
    for (; it != names_.end(); ++it)
    {
      delete it->second;
    }
    names_.clear();
    symbols_.clear();
    atomic_numbers_.clear();
  }

} // namespace OpenMS
