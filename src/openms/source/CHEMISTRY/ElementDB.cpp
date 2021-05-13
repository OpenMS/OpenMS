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
// $Authors: Andreas Bertsch, Timo Sachsenberg, Chris Bielow $
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
    if (names_.count(name) == 1)
    {
      return names_.at(name);
    }
    else
    {
      if (symbols_.count(name) == 1)
      {
        return symbols_.at(name);
      }
    }
    return nullptr;
  }

  const Element* ElementDB::getElement(unsigned int atomic_number) const
  {
    if (atomic_numbers_.count(atomic_number) == 1)
    {
      return atomic_numbers_.at(atomic_number);
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

    map<unsigned int, double> hydrogen_abundance = {{(unsigned int) 1, 0.999885}, {(unsigned int) 2, 0.000115}, {(unsigned int) 3, 0.0}, };
    map<unsigned int, double> hydrogen_mass = {{(unsigned int) 1, 1.0078250319}, {(unsigned int) 2, 2.01410178}, {(unsigned int) 3, 3.01604927}};
    buildElement_("Hydrogen", "H", (unsigned int) 1, hydrogen_abundance, hydrogen_mass);


    map<unsigned int, double> helium_abundance = {{(unsigned int) 3, 1.34e-06}, {(unsigned int) 4, 0.9999986599999999}, };
    map<unsigned int, double> helium_mass = {{(unsigned int) 3, 3.0160293191}, {(unsigned int) 4, 4.00260325415}};
    buildElement_("Helium", "He", (unsigned int) 2, helium_abundance, helium_mass);


    map<unsigned int, double> lithium_abundance = {{(unsigned int) 6, 0.0759}, {(unsigned int) 7, 0.9240999999999999}, };
    map<unsigned int, double> lithium_mass = {{(unsigned int) 6, 6.015122}, {(unsigned int) 7, 7.016004}};
    buildElement_("Lithium", "Li", (unsigned int) 3, lithium_abundance, lithium_mass);


    map<unsigned int, double> beryllium_abundance = {{(unsigned int) 9, 1.0}, };
    map<unsigned int, double> beryllium_mass = {{(unsigned int) 9, 9.0121822}};
    buildElement_("Beryllium", "Be", (unsigned int) 4, beryllium_abundance, beryllium_mass);


    map<unsigned int, double> bor_abundance = {{(unsigned int) 10, 0.19899999999999998}, {(unsigned int) 11, 0.8009999999999999}, };
    map<unsigned int, double> bor_mass = {{(unsigned int) 10, 10.012937000000001}, {(unsigned int) 11, 11.009304999999999}};
    buildElement_("Bor", "B", (unsigned int) 5, bor_abundance, bor_mass);


    map<unsigned int, double> carbon_abundance = {{(unsigned int) 12, 0.9893000000000001}, {(unsigned int) 13, 0.010700000000000001}, };
    map<unsigned int, double> carbon_mass = {{(unsigned int) 12, 12.0}, {(unsigned int) 13, 13.003355000000001}};
    buildElement_("Carbon", "C", (unsigned int) 6, carbon_abundance, carbon_mass);


    map<unsigned int, double> nitrogen_abundance = {{(unsigned int) 14, 0.9963200000000001}, {(unsigned int) 15, 0.00368}, };
    map<unsigned int, double> nitrogen_mass = {{(unsigned int) 14, 14.003074}, {(unsigned int) 15, 15.000109}};
    buildElement_("Nitrogen", "N", (unsigned int) 7, nitrogen_abundance, nitrogen_mass);


    map<unsigned int, double> oxygen_abundance = {{(unsigned int) 16, 0.9975700000000001}, {(unsigned int) 17, 0.00037999999999999997}, {(unsigned int) 18, 0.0020499999999999997}, };
    map<unsigned int, double> oxygen_mass = {{(unsigned int) 16, 15.994915000000001}, {(unsigned int) 17, 16.999132}, {(unsigned int) 18, 17.999168999999998}};
    buildElement_("Oxygen", "O", (unsigned int) 8, oxygen_abundance, oxygen_mass);


    map<unsigned int, double> fluorine_abundance = {{(unsigned int) 19, 1.0}, };
    map<unsigned int, double> fluorine_mass = {{(unsigned int) 19, 18.99840322}};
    buildElement_("Fluorine", "F", (unsigned int) 9, fluorine_abundance, fluorine_mass);


    map<unsigned int, double> argon_abundance = {{(unsigned int) 36, 0.003336}, {(unsigned int) 38, 0.000629}, {(unsigned int) 40, 0.996035}, };
    map<unsigned int, double> argon_mass = {{(unsigned int) 36, 35.967545106000003}, {(unsigned int) 38, 37.9627324}, {(unsigned int) 40, 39.9623831225}};
    buildElement_("Argon", "Ar", (unsigned int) 18, argon_abundance, argon_mass);


    map<unsigned int, double> titanium_abundance = {{(unsigned int) 46, 0.0825}, {(unsigned int) 47, 0.07440000000000001}, {(unsigned int) 48, 0.7372}, {(unsigned int) 49, 0.0541}, {(unsigned int) 50, 0.0518}, };
    map<unsigned int, double> titanium_mass = {{(unsigned int) 46, 45.952631599999997}, {(unsigned int) 47, 46.951763100000001}, {(unsigned int) 48, 47.947946299999998}, {(unsigned int) 49, 48.947870000000002}, {(unsigned int) 50, 49.944791199999997}};
    buildElement_("Titanium", "Ti", (unsigned int) 22, titanium_abundance, titanium_mass);


    map<unsigned int, double> sodium_abundance = {{(unsigned int) 23, 1.0}, };
    map<unsigned int, double> sodium_mass = {{(unsigned int) 23, 22.989769280899999}};
    buildElement_("Sodium", "Na", (unsigned int) 11, sodium_abundance, sodium_mass);


    map<unsigned int, double> magnesium_abundance = {{(unsigned int) 24, 0.7898999999999999}, {(unsigned int) 25, 0.1}, {(unsigned int) 26, 0.1101}, };
    map<unsigned int, double> magnesium_mass = {{(unsigned int) 24, 23.985042}, {(unsigned int) 25, 24.985837}, {(unsigned int) 26, 25.982593000000001}};
    buildElement_("Magnesium", "Mg", (unsigned int) 12, magnesium_abundance, magnesium_mass);


    map<unsigned int, double> aluminium_abundance = {{(unsigned int) 27, 1.0}, };
    map<unsigned int, double> aluminium_mass = {{(unsigned int) 27, 26.981538629999999}};
    buildElement_("Aluminium", "Al", (unsigned int) 13, aluminium_abundance, aluminium_mass);


    map<unsigned int, double> silicon_abundance = {{(unsigned int) 28, 0.9220999999999999}, {(unsigned int) 29, 0.0467}, {(unsigned int) 30, 0.031}, };
    map<unsigned int, double> silicon_mass = {{(unsigned int) 28, 27.976926532499999}, {(unsigned int) 29, 28.9764947}, {(unsigned int) 30, 29.973770170000002}};
    buildElement_("Silicon", "Si", (unsigned int) 14, silicon_abundance, silicon_mass);


    map<unsigned int, double> phosphorus_abundance = {{(unsigned int) 31, 1.0}, };
    map<unsigned int, double> phosphorus_mass = {{(unsigned int) 31, 30.973761490000001}};
    buildElement_("Phosphorus", "P", (unsigned int) 15, phosphorus_abundance, phosphorus_mass);


    map<unsigned int, double> sulfur_abundance = {{(unsigned int) 32, 0.9493}, {(unsigned int) 33, 0.0076}, {(unsigned int) 34, 0.0429}, {(unsigned int) 36, 0.0002}, };
    map<unsigned int, double> sulfur_mass = {{(unsigned int) 32, 31.972070729999999}, {(unsigned int) 33, 32.971457999999998}, {(unsigned int) 34, 33.967866999999998}, {(unsigned int) 36, 35.967081}};
    buildElement_("Sulfur", "S", (unsigned int) 16, sulfur_abundance, sulfur_mass);


    map<unsigned int, double> chlorine_abundance = {{(unsigned int) 35, 0.7576}, {(unsigned int) 37, 0.24239999999999998}, };
    map<unsigned int, double> chlorine_mass = {{(unsigned int) 35, 34.968852679999998}, {(unsigned int) 37, 36.965902589999999}};
    buildElement_("Chlorine", "Cl", (unsigned int) 17, chlorine_abundance, chlorine_mass);


    map<unsigned int, double> potassium_abundance = {{(unsigned int) 39, 0.932581}, {(unsigned int) 40, 0.000117}, {(unsigned int) 41, 0.067302}, };
    map<unsigned int, double> potassium_mass = {{(unsigned int) 39, 38.963706680000001}, {(unsigned int) 40, 39.963998480000001}, {(unsigned int) 41, 40.961825760000004}};
    buildElement_("Potassium", "K", (unsigned int) 19, potassium_abundance, potassium_mass);


    map<unsigned int, double> calcium_abundance = {{(unsigned int) 40, 0.96941}, {(unsigned int) 42, 0.00647}, {(unsigned int) 43, 0.00135}, {(unsigned int) 44, 0.02086}, {(unsigned int) 46, 4e-05}, {(unsigned int) 48, 0.00187}, };
    map<unsigned int, double> calcium_mass = {{(unsigned int) 40, 39.962590980000002}, {(unsigned int) 42, 41.958618010000002}, {(unsigned int) 43, 42.958766599999997}, {(unsigned int) 44, 43.955481800000001}, {(unsigned int) 46, 45.953692599999997}, {(unsigned int) 48, 47.952534}};
    buildElement_("Calcium", "Ca", (unsigned int) 20, calcium_abundance, calcium_mass);


    map<unsigned int, double> scandium_abundance = {{(unsigned int) 45, 1.0}, };
    map<unsigned int, double> scandium_mass = {{(unsigned int) 45, 44.955910000000003}};
    buildElement_("Scandium", "Sc", (unsigned int) 21, scandium_abundance, scandium_mass);


    map<unsigned int, double> vanadium_abundance = {{(unsigned int) 50, 0.0025}, {(unsigned int) 51, 0.9975}, };
    map<unsigned int, double> vanadium_mass = {{(unsigned int) 50, 49.947158500000001}, {(unsigned int) 51, 50.943959499999998}};
    buildElement_("Vanadium", "V", (unsigned int) 23, vanadium_abundance, vanadium_mass);


    map<unsigned int, double> chromium_abundance = {{(unsigned int) 50, 0.043449999999999996}, {(unsigned int) 52, 0.83789}, {(unsigned int) 53, 0.09501}, {(unsigned int) 54, 0.02365}, };
    map<unsigned int, double> chromium_mass = {{(unsigned int) 50, 49.946044200000003}, {(unsigned int) 52, 51.940507500000003}, {(unsigned int) 53, 52.940649399999998}, {(unsigned int) 54, 53.938880400000002}};
    buildElement_("Chromium", "Cr", (unsigned int) 24, chromium_abundance, chromium_mass);


    map<unsigned int, double> tellurium_abundance = {{(unsigned int) 120, 0.0009}, {(unsigned int) 122, 0.0255}, {(unsigned int) 124, 0.047400000000000005}, {(unsigned int) 125, 0.0707}, {(unsigned int) 126, 0.1884}, {(unsigned int) 128, 0.31739999999999996}, {(unsigned int) 130, 0.3408}, };
    map<unsigned int, double> tellurium_mass = {{(unsigned int) 120, 119.904020000000003}, {(unsigned int) 122, 121.9030439}, {(unsigned int) 124, 123.902817900000002}, {(unsigned int) 125, 124.904430700000006}, {(unsigned int) 126, 125.903311700000003}, {(unsigned int) 128, 127.904463100000001}, {(unsigned int) 130, 129.906224400000014}};
    buildElement_("Tellurium", "Te", (unsigned int) 52, tellurium_abundance, tellurium_mass);


    map<unsigned int, double> barium_abundance = {{(unsigned int) 132, 0.00101}, {(unsigned int) 134, 0.024169999999999997}, {(unsigned int) 135, 0.06591999999999999}, {(unsigned int) 136, 0.07854}, {(unsigned int) 137, 0.11231999999999999}, {(unsigned int) 138, 0.71698}, };
    map<unsigned int, double> barium_mass = {{(unsigned int) 132, 131.9050613}, {(unsigned int) 134, 133.904508399999997}, {(unsigned int) 135, 134.905688599999991}, {(unsigned int) 136, 135.904575899999998}, {(unsigned int) 137, 136.905827399999993}, {(unsigned int) 138, 137.905247199999991}};
    buildElement_("Barium", "Ba", (unsigned int) 56, barium_abundance, barium_mass);


    map<unsigned int, double> manganese_abundance = {{(unsigned int) 55, 1.0}, };
    map<unsigned int, double> manganese_mass = {{(unsigned int) 55, 54.938049999999997}};
    buildElement_("Manganese", "Mn", (unsigned int) 25, manganese_abundance, manganese_mass);


    map<unsigned int, double> ferrum_abundance = {{(unsigned int) 54, 0.058449999999999995}, {(unsigned int) 56, 0.91754}, {(unsigned int) 57, 0.021191}, {(unsigned int) 58, 0.002819}, };
    map<unsigned int, double> ferrum_mass = {{(unsigned int) 54, 53.939610500000001}, {(unsigned int) 56, 55.934937499999997}, {(unsigned int) 57, 56.935394000000002}, {(unsigned int) 58, 57.933275600000002}};
    buildElement_("Ferrum", "Fe", (unsigned int) 26, ferrum_abundance, ferrum_mass);


    map<unsigned int, double> cobalt_abundance = {{(unsigned int) 59, 1.0}, };
    map<unsigned int, double> cobalt_mass = {{(unsigned int) 59, 58.933194999999998}};
    buildElement_("Cobalt", "Co", (unsigned int) 27, cobalt_abundance, cobalt_mass);


    map<unsigned int, double> nickel_abundance = {{(unsigned int) 58, 0.680169}, {(unsigned int) 60, 0.262231}, {(unsigned int) 61, 0.011399}, {(unsigned int) 62, 0.036345}, {(unsigned int) 64, 0.009256}, };
    map<unsigned int, double> nickel_mass = {{(unsigned int) 58, 57.935347999999998}, {(unsigned int) 60, 59.930790999999999}, {(unsigned int) 61, 60.931060000000002}, {(unsigned int) 62, 61.928348999999997}, {(unsigned int) 64, 63.927970000000002}};
    buildElement_("Nickel", "Ni", (unsigned int) 28, nickel_abundance, nickel_mass);


    map<unsigned int, double> copper_abundance = {{(unsigned int) 63, 0.6917}, {(unsigned int) 65, 0.30829999999999996}, };
    map<unsigned int, double> copper_mass = {{(unsigned int) 63, 62.929600999999998}, {(unsigned int) 65, 64.927794000000006}};
    buildElement_("Copper", "Cu", (unsigned int) 29, copper_abundance, copper_mass);


    map<unsigned int, double> zinc_abundance = {{(unsigned int) 64, 0.4863}, {(unsigned int) 66, 0.27899999999999997}, {(unsigned int) 67, 0.040999999999999995}, {(unsigned int) 68, 0.1875}, {(unsigned int) 70, 0.0062}, };
    map<unsigned int, double> zinc_mass = {{(unsigned int) 64, 63.929147}, {(unsigned int) 66, 65.926036999999994}, {(unsigned int) 67, 66.927131000000003}, {(unsigned int) 68, 67.924847999999997}, {(unsigned int) 70, 69.925325000000001}};
    buildElement_("Zinc", "Zn", (unsigned int) 30, zinc_abundance, zinc_mass);


    map<unsigned int, double> gallium_abundance = {{(unsigned int) 69, 0.60108}, {(unsigned int) 71, 0.39892000000000005}, };
    map<unsigned int, double> gallium_mass = {{(unsigned int) 69, 68.925573600000007}, {(unsigned int) 71, 70.924701299999995}};
    buildElement_("Gallium", "Ga", (unsigned int) 31, gallium_abundance, gallium_mass);


    map<unsigned int, double> germanium_abundance = {{(unsigned int) 70, 0.20379999999999998}, {(unsigned int) 72, 0.2731}, {(unsigned int) 73, 0.0776}, {(unsigned int) 74, 0.36719999999999997}, };
    map<unsigned int, double> germanium_mass = {{(unsigned int) 70, 69.924247399999999}, {(unsigned int) 72, 71.922075800000002}, {(unsigned int) 73, 72.9234589}, {(unsigned int) 74, 73.921177799999995}};
    buildElement_("Germanium", "Ge", (unsigned int) 32, germanium_abundance, germanium_mass);


    map<unsigned int, double> arsenic_abundance = {{(unsigned int) 75, 1.0}, };
    map<unsigned int, double> arsenic_mass = {{(unsigned int) 75, 74.921596500000007}};
    buildElement_("Arsenic", "As", (unsigned int) 33, arsenic_abundance, arsenic_mass);


    map<unsigned int, double> rubidium_abundance = {{(unsigned int) 85, 0.7217}, };
    map<unsigned int, double> rubidium_mass = {{(unsigned int) 85, 84.911789737999996}};
    buildElement_("Rubidium", "Rb", (unsigned int) 37, rubidium_abundance, rubidium_mass);


    map<unsigned int, double> strontium_abundance = {{(unsigned int) 84, 0.005600000000000001}, {(unsigned int) 86, 0.0986}, {(unsigned int) 87, 0.07}, {(unsigned int) 88, 0.8258}, };
    map<unsigned int, double> strontium_mass = {{(unsigned int) 84, 83.913425000000004}, {(unsigned int) 86, 85.909260730900002}, {(unsigned int) 87, 86.908877497000006}, {(unsigned int) 88, 87.905612257100003}};
    buildElement_("Strontium", "Sr", (unsigned int) 38, strontium_abundance, strontium_mass);


    map<unsigned int, double> yttrium_abundance = {{(unsigned int) 89, 1.0}, };
    map<unsigned int, double> yttrium_mass = {{(unsigned int) 89, 88.905850000000001}};
    buildElement_("Yttrium", "Y", (unsigned int) 39, yttrium_abundance, yttrium_mass);


    map<unsigned int, double> zirconium_abundance = {{(unsigned int) 90, 0.5145000000000001}, {(unsigned int) 91, 0.11220000000000001}, {(unsigned int) 92, 0.17149999999999999}, {(unsigned int) 94, 0.17379999999999998}, };
    map<unsigned int, double> zirconium_mass = {{(unsigned int) 90, 89.9047044}, {(unsigned int) 91, 90.905645800000002}, {(unsigned int) 92, 91.905040799999995}, {(unsigned int) 94, 93.906315199999995}};
    buildElement_("Zirconium", "Zr", (unsigned int) 40, zirconium_abundance, zirconium_mass);


    map<unsigned int, double> nibium_abundance = {{(unsigned int) 93, 1.0}, };
    map<unsigned int, double> nibium_mass = {{(unsigned int) 93, 92.906378099999998}};
    buildElement_("Nibium", "Nb", (unsigned int) 41, nibium_abundance, nibium_mass);


    map<unsigned int, double> ruthenium_abundance = {{(unsigned int) 96, 0.0554}, {(unsigned int) 98, 0.0187}, {(unsigned int) 99, 0.1276}, {(unsigned int) 100, 0.126}, {(unsigned int) 101, 0.17059999999999997}, {(unsigned int) 102, 0.3155}, {(unsigned int) 104, 0.1862}, };
    map<unsigned int, double> ruthenium_mass = {{(unsigned int) 96, 95.907597999999993}, {(unsigned int) 98, 97.905287000000001}, {(unsigned int) 99, 98.9059393}, {(unsigned int) 100, 99.904219499999996}, {(unsigned int) 101, 100.905582100000004}, {(unsigned int) 102, 101.904349300000007}, {(unsigned int) 104, 103.905433000000002}};
    buildElement_("Ruthenium", "Ru", (unsigned int) 44, ruthenium_abundance, ruthenium_mass);


    map<unsigned int, double> tin_abundance = {{(unsigned int) 112, 0.0097}, {(unsigned int) 114, 0.0066}, {(unsigned int) 115, 0.0034000000000000002}, {(unsigned int) 116, 0.1454}, {(unsigned int) 117, 0.0768}, {(unsigned int) 118, 0.2422}, {(unsigned int) 119, 0.0859}, {(unsigned int) 120, 0.3258}, {(unsigned int) 122, 0.0463}, {(unsigned int) 124, 0.0579}, };
    map<unsigned int, double> tin_mass = {{(unsigned int) 112, 111.904818000000006}, {(unsigned int) 114, 113.902777900000004}, {(unsigned int) 115, 114.903341999999995}, {(unsigned int) 116, 115.901741000000001}, {(unsigned int) 117, 116.902951999999999}, {(unsigned int) 118, 117.901602999999994}, {(unsigned int) 119, 118.903307999999996}, {(unsigned int) 120, 119.902194699999996}, {(unsigned int) 122, 121.903439000000006}, {(unsigned int) 124, 123.905273899999997}};
    buildElement_("Tin", "Sn", (unsigned int) 50, tin_abundance, tin_mass);


    map<unsigned int, double> antimony_abundance = {{(unsigned int) 121, 0.5721}, {(unsigned int) 123, 0.4279}, };
    map<unsigned int, double> antimony_mass = {{(unsigned int) 121, 120.903815699999996}, {(unsigned int) 123, 122.904213999999996}};
    buildElement_("Antimony", "Sb", (unsigned int) 51, antimony_abundance, antimony_mass);


    map<unsigned int, double> selenium_abundance = {{(unsigned int) 74, 0.00889}, {(unsigned int) 76, 0.09366}, {(unsigned int) 77, 0.07635}, {(unsigned int) 78, 0.23772}, {(unsigned int) 80, 0.49607}, {(unsigned int) 82, 0.08731}, };
    map<unsigned int, double> selenium_mass = {{(unsigned int) 74, 73.922476399999994}, {(unsigned int) 76, 75.919213600000006}, {(unsigned int) 77, 76.919914000000006}, {(unsigned int) 78, 77.917309099999997}, {(unsigned int) 80, 79.916521299999999}, {(unsigned int) 82, 81.916699399999999}};
    buildElement_("Selenium", "Se", (unsigned int) 34, selenium_abundance, selenium_mass);


    map<unsigned int, double> bromine_abundance = {{(unsigned int) 79, 0.5069}, {(unsigned int) 81, 0.49310000000000004}, };
    map<unsigned int, double> bromine_mass = {{(unsigned int) 79, 78.918337100000002}, {(unsigned int) 81, 80.916290599999996}};
    buildElement_("Bromine", "Br", (unsigned int) 35, bromine_abundance, bromine_mass);


    map<unsigned int, double> krypton_abundance = {{(unsigned int) 78, 0.0034999999999999996}, {(unsigned int) 80, 0.0225}, {(unsigned int) 82, 0.11599999999999999}, {(unsigned int) 83, 0.115}, {(unsigned int) 84, 0.57}, {(unsigned int) 86, 0.17300000000000001}, };
    map<unsigned int, double> krypton_mass = {{(unsigned int) 78, 77.920400000000001}, {(unsigned int) 80, 79.916380000000004}, {(unsigned int) 82, 81.913482000000002}, {(unsigned int) 83, 82.914135000000002}, {(unsigned int) 84, 83.911507}, {(unsigned int) 86, 85.910616000000005}};
    buildElement_("Krypton", "Kr", (unsigned int) 36, krypton_abundance, krypton_mass);


    map<unsigned int, double> molybdenum_abundance = {{(unsigned int) 92, 0.1484}, {(unsigned int) 94, 0.0925}, {(unsigned int) 95, 0.1592}, {(unsigned int) 96, 0.1668}, {(unsigned int) 97, 0.0955}, {(unsigned int) 98, 0.2413}, {(unsigned int) 100, 0.09630000000000001}, };
    map<unsigned int, double> molybdenum_mass = {{(unsigned int) 92, 91.906809999999993}, {(unsigned int) 94, 93.905088000000006}, {(unsigned int) 95, 94.905840999999995}, {(unsigned int) 96, 95.904679000000002}, {(unsigned int) 97, 96.906020999999996}, {(unsigned int) 98, 97.905407999999994}, {(unsigned int) 100, 99.907477}};
    buildElement_("Molybdenum", "Mo", (unsigned int) 42, molybdenum_abundance, molybdenum_mass);


    map<unsigned int, double> technitium_abundance = {{(unsigned int) 97, 0.0}, {(unsigned int) 98, 0.0}, {(unsigned int) 99, 0.0}, };
    map<unsigned int, double> technitium_mass = {{(unsigned int) 97, 96.906363999999996}, {(unsigned int) 98, 97.907214999999994}, {(unsigned int) 99, 98.906254000000004}};
    buildElement_("Technitium", "Tc", (unsigned int) 43, technitium_abundance, technitium_mass);


    map<unsigned int, double> rhodium_abundance = {{(unsigned int) 103, 1.0}, };
    map<unsigned int, double> rhodium_mass = {{(unsigned int) 103, 102.905500000000004}};
    buildElement_("Rhodium", "Rh", (unsigned int) 45, rhodium_abundance, rhodium_mass);


    map<unsigned int, double> palladium_abundance = {{(unsigned int) 102, 0.0102}, {(unsigned int) 104, 0.1114}, {(unsigned int) 105, 0.22329999999999997}, {(unsigned int) 106, 0.2733}, {(unsigned int) 108, 0.2646}, {(unsigned int) 110, 0.11720000000000001}, };
    map<unsigned int, double> palladium_mass = {{(unsigned int) 102, 101.905608999999998}, {(unsigned int) 104, 103.904036000000005}, {(unsigned int) 105, 104.905085}, {(unsigned int) 106, 105.903486000000001}, {(unsigned int) 108, 107.903891999999999}, {(unsigned int) 110, 109.905152999999999}};
    buildElement_("Palladium", "Pd", (unsigned int) 46, palladium_abundance, palladium_mass);


    map<unsigned int, double> silver_abundance = {{(unsigned int) 107, 0.51839}, {(unsigned int) 109, 0.48161000000000004}, };
    map<unsigned int, double> silver_mass = {{(unsigned int) 107, 106.905092999999994}, {(unsigned int) 109, 108.904756000000006}};
    buildElement_("Silver", "Ag", (unsigned int) 47, silver_abundance, silver_mass);


    map<unsigned int, double> cadmium_abundance = {{(unsigned int) 106, 0.0125}, {(unsigned int) 108, 0.0089}, {(unsigned int) 110, 0.1249}, {(unsigned int) 111, 0.128}, {(unsigned int) 112, 0.2413}, {(unsigned int) 113, 0.1222}, {(unsigned int) 114, 0.2873}, {(unsigned int) 116, 0.07490000000000001}, };
    map<unsigned int, double> cadmium_mass = {{(unsigned int) 106, 105.906458000000001}, {(unsigned int) 108, 107.904184000000001}, {(unsigned int) 110, 109.903002099999995}, {(unsigned int) 111, 110.904178099999996}, {(unsigned int) 112, 111.902757800000003}, {(unsigned int) 113, 112.904401699999994}, {(unsigned int) 114, 113.903358499999996}, {(unsigned int) 116, 115.904756000000006}};
    buildElement_("Cadmium", "Cd", (unsigned int) 48, cadmium_abundance, cadmium_mass);


    map<unsigned int, double> indium_abundance = {{(unsigned int) 113, 0.0429}, {(unsigned int) 115, 0.9571}, };
    map<unsigned int, double> indium_mass = {{(unsigned int) 113, 112.904060000000001}, {(unsigned int) 115, 114.903878000000006}};
    buildElement_("Indium", "In", (unsigned int) 49, indium_abundance, indium_mass);


    map<unsigned int, double> iodine_abundance = {{(unsigned int) 127, 1.0}, };
    map<unsigned int, double> iodine_mass = {{(unsigned int) 127, 126.904472999999996}};
    buildElement_("Iodine", "I", (unsigned int) 53, iodine_abundance, iodine_mass);


    map<unsigned int, double> xenon_abundance = {{(unsigned int) 128, 0.0191}, {(unsigned int) 129, 0.264}, {(unsigned int) 130, 0.040999999999999995}, {(unsigned int) 131, 0.212}, {(unsigned int) 132, 0.26899999999999996}, {(unsigned int) 134, 0.10400000000000001}, {(unsigned int) 136, 0.08900000000000001}, };
    map<unsigned int, double> xenon_mass = {{(unsigned int) 128, 127.903531000000001}, {(unsigned int) 129, 128.904779999999988}, {(unsigned int) 130, 129.903509000000014}, {(unsigned int) 131, 130.90507199999999}, {(unsigned int) 132, 131.904144000000002}, {(unsigned int) 134, 133.905394999999999}, {(unsigned int) 136, 135.90721400000001}};
    buildElement_("Xenon", "Xe", (unsigned int) 54, xenon_abundance, xenon_mass);


    map<unsigned int, double> caesium_abundance = {{(unsigned int) 133, 1.0}, };
    map<unsigned int, double> caesium_mass = {{(unsigned int) 133, 132.905451932999995}};
    buildElement_("Caesium", "Cs", (unsigned int) 55, caesium_abundance, caesium_mass);


    map<unsigned int, double> cerium_abundance = {{(unsigned int) 136, 0.00185}, {(unsigned int) 138, 0.00251}, {(unsigned int) 140, 0.8845000000000001}, {(unsigned int) 142, 0.11114}, };
    map<unsigned int, double> cerium_mass = {{(unsigned int) 136, 135.907172000000003}, {(unsigned int) 138, 137.905991}, {(unsigned int) 140, 139.905438699999991}, {(unsigned int) 142, 141.909244000000001}};
    buildElement_("Cerium", "Ce", (unsigned int) 58, cerium_abundance, cerium_mass);


    map<unsigned int, double> praseodymium_abundance = {{(unsigned int) 141, 1.0}, };
    map<unsigned int, double> praseodymium_mass = {{(unsigned int) 141, 140.907646999999997}};
    buildElement_("Praseodymium", "Pr", (unsigned int) 59, praseodymium_abundance, praseodymium_mass);


    map<unsigned int, double> gadolinium_abundance = {{(unsigned int) 152, 0.002}, {(unsigned int) 154, 0.0218}, {(unsigned int) 155, 0.14800000000000002}, {(unsigned int) 156, 0.2047}, {(unsigned int) 157, 0.1565}, {(unsigned int) 158, 0.2484}, {(unsigned int) 160, 0.2186}, };
    map<unsigned int, double> gadolinium_mass = {{(unsigned int) 152, 151.919791000000004}, {(unsigned int) 154, 153.920865600000013}, {(unsigned int) 155, 154.92262199999999}, {(unsigned int) 156, 155.922122699999989}, {(unsigned int) 157, 156.923960099999988}, {(unsigned int) 158, 157.924103900000006}, {(unsigned int) 160, 159.927054099999992}};
    buildElement_("Gadolinium", "Gd", (unsigned int) 64, gadolinium_abundance, gadolinium_mass);


    map<unsigned int, double> hafnium_abundance = {{(unsigned int) 176, 0.0526}, {(unsigned int) 177, 0.18600000000000003}, {(unsigned int) 178, 0.2728}, {(unsigned int) 179, 0.1362}, {(unsigned int) 180, 0.3508}, };
    map<unsigned int, double> hafnium_mass = {{(unsigned int) 176, 175.941408599999988}, {(unsigned int) 177, 176.943220700000012}, {(unsigned int) 178, 177.943698799999993}, {(unsigned int) 179, 178.945816100000002}, {(unsigned int) 180, 179.946550000000002}};
    buildElement_("Hafnium", "Hf", (unsigned int) 72, hafnium_abundance, hafnium_mass);


    map<unsigned int, double> tantalum_abundance = {{(unsigned int) 181, 1.0}, };
    map<unsigned int, double> tantalum_mass = {{(unsigned int) 181, 180.947995800000001}};
    buildElement_("Tantalum", "Ta", (unsigned int) 73, tantalum_abundance, tantalum_mass);


    map<unsigned int, double> platinum_abundance = {{(unsigned int) 192, 0.00782}, {(unsigned int) 194, 0.32966999999999996}, {(unsigned int) 195, 0.33832}, {(unsigned int) 196, 0.25242000000000003}, {(unsigned int) 198, 0.07163}, };
    map<unsigned int, double> platinum_mass = {{(unsigned int) 192, 191.961038000000002}, {(unsigned int) 194, 193.962680299999988}, {(unsigned int) 195, 194.964791100000014}, {(unsigned int) 196, 195.964951500000012}, {(unsigned int) 198, 197.967893000000004}};
    buildElement_("Platinum", "Pt", (unsigned int) 78, platinum_abundance, platinum_mass);


    map<unsigned int, double> tungsten_abundance = {{(unsigned int) 180, 0.0012}, {(unsigned int) 182, 0.265}, {(unsigned int) 183, 0.1431}, {(unsigned int) 184, 0.3064}, {(unsigned int) 186, 0.2843}, };
    map<unsigned int, double> tungsten_mass = {{(unsigned int) 180, 179.946704000000011}, {(unsigned int) 182, 181.948204199999992}, {(unsigned int) 183, 182.950222999999994}, {(unsigned int) 184, 183.950930999999997}, {(unsigned int) 186, 185.954364099999992}};
    buildElement_("Tungsten", "W", (unsigned int) 74, tungsten_abundance, tungsten_mass);


    map<unsigned int, double> gold_abundance = {{(unsigned int) 197, 1.0}, };
    map<unsigned int, double> gold_mass = {{(unsigned int) 197, 196.96655100000001}};
    buildElement_("Gold", "Au", (unsigned int) 79, gold_abundance, gold_mass);


    map<unsigned int, double> mercury_abundance = {{(unsigned int) 196, 0.0015}, {(unsigned int) 198, 0.09970000000000001}, {(unsigned int) 199, 0.16870000000000002}, {(unsigned int) 200, 0.231}, {(unsigned int) 201, 0.1318}, {(unsigned int) 202, 0.2986}, {(unsigned int) 204, 0.0687}, };
    map<unsigned int, double> mercury_mass = {{(unsigned int) 196, 195.965833000000004}, {(unsigned int) 198, 197.966768999999999}, {(unsigned int) 199, 198.968279899999999}, {(unsigned int) 200, 199.968325999999991}, {(unsigned int) 201, 200.970302299999986}, {(unsigned int) 202, 201.970642999999996}, {(unsigned int) 204, 203.973493899999994}};
    buildElement_("Mercury", "Hg", (unsigned int) 80, mercury_abundance, mercury_mass);


    map<unsigned int, double> thallium_abundance = {{(unsigned int) 203, 0.2952}, {(unsigned int) 205, 0.7048000000000001}, };
    map<unsigned int, double> thallium_mass = {{(unsigned int) 203, 202.972344200000009}, {(unsigned int) 205, 204.97442749999999}};
    buildElement_("Thallium", "Tl", (unsigned int) 81, thallium_abundance, thallium_mass);


    map<unsigned int, double> lead_abundance = {{(unsigned int) 204, 0.013999999999999999}, {(unsigned int) 206, 0.24100000000000002}, {(unsigned int) 207, 0.221}, {(unsigned int) 208, 0.524}, };
    map<unsigned int, double> lead_mass = {{(unsigned int) 204, 203.973043600000011}, {(unsigned int) 206, 205.974465299999991}, {(unsigned int) 207, 206.975896900000009}, {(unsigned int) 208, 207.976653800000008}};
    buildElement_("Lead", "Pb", (unsigned int) 82, lead_abundance, lead_mass);


    map<unsigned int, double> bismuth_abundance = {{(unsigned int) 209, 1.0}, };
    map<unsigned int, double> bismuth_mass = {{(unsigned int) 209, 208.980398699999995}};
    buildElement_("Bismuth", "Bi", (unsigned int) 83, bismuth_abundance, bismuth_mass);


    map<unsigned int, double> rhenium_abundance = {{(unsigned int) 185, 0.374}, {(unsigned int) 187, 0.626}, };
    map<unsigned int, double> rhenium_mass = {{(unsigned int) 185, 184.952955000000003}, {(unsigned int) 187, 186.95575310000001}};
    buildElement_("Rhenium", "Re", (unsigned int) 75, rhenium_abundance, rhenium_mass);


    map<unsigned int, double> neodymium_abundance = {{(unsigned int) 142, 0.272}, {(unsigned int) 143, 0.122}, {(unsigned int) 144, 0.23800000000000002}, {(unsigned int) 145, 0.083}, {(unsigned int) 146, 0.172}, {(unsigned int) 148, 0.057999999999999996}, {(unsigned int) 150, 0.055999999999999994}, };
    map<unsigned int, double> neodymium_mass = {{(unsigned int) 142, 141.907723299999987}, {(unsigned int) 143, 142.909814299999994}, {(unsigned int) 144, 143.910087299999987}, {(unsigned int) 145, 144.912573600000002}, {(unsigned int) 146, 145.913116900000006}, {(unsigned int) 148, 147.916892999999988}, {(unsigned int) 150, 149.920891000000012}};
    buildElement_("Neodymium", "Nd", (unsigned int) 60, neodymium_abundance, neodymium_mass);


    map<unsigned int, double> thorium_abundance = {{(unsigned int) 230, 0.0002}, {(unsigned int) 232, 0.9998}, };
    map<unsigned int, double> thorium_mass = {{(unsigned int) 230, 230.033133800000002}, {(unsigned int) 232, 232.038055299999996}};
    buildElement_("Thorium", "Th", (unsigned int) 90, thorium_abundance, thorium_mass);


    map<unsigned int, double> lanthanum_abundance = {{(unsigned int) 138, 0.00089}, {(unsigned int) 139, 0.99911}, };
    map<unsigned int, double> lanthanum_mass = {{(unsigned int) 138, 137.907112000000012}, {(unsigned int) 139, 138.906353300000006}};
    buildElement_("Lanthanum", "La", (unsigned int) 57, lanthanum_abundance, lanthanum_mass);


    map<unsigned int, double> samarium_abundance = {{(unsigned int) 144, 0.0308}, {(unsigned int) 147, 0.15}, {(unsigned int) 148, 0.1125}, {(unsigned int) 149, 0.1382}, {(unsigned int) 150, 0.0737}, {(unsigned int) 152, 0.26739999999999997}, {(unsigned int) 154, 0.2274}, };
    map<unsigned int, double> samarium_mass = {{(unsigned int) 144, 143.911999000000009}, {(unsigned int) 147, 146.9148979}, {(unsigned int) 148, 147.914822700000002}, {(unsigned int) 149, 148.917184700000007}, {(unsigned int) 150, 149.917275499999988}, {(unsigned int) 152, 151.919732399999987}, {(unsigned int) 154, 153.92220929999999}};
    buildElement_("Samarium", "Sm", (unsigned int) 62, samarium_abundance, samarium_mass);

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
    
    vector<unsigned int> keys;
    for (map<unsigned int, double>::const_iterator it = abundance.begin(); it != abundance.end(); ++it)
    {
      keys.push_back(it->first);
    }

    // calculate weighted average
    for (vector<unsigned int>::iterator it = keys.begin(); it != keys.end(); ++it)
    {
      dist.push_back(Peak1D(mass.at(*it) , abundance.at(*it)));
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
