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
    storeElements();
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

  const Map<String, const Element*>& ElementDB::getNames() const
  {
    return names_;
  }

  const Map<String, const Element*>& ElementDB::getSymbols() const
  {
    return symbols_;
  }

  const Map<UInt, const Element*>& ElementDB::getAtomicNumbers() const
  {
    return atomic_numbers_;
  }

  const Element* ElementDB::getElement(const String& name) const
  {
    if (names_.has(name))
    {
      return names_[name];
    }
    else
    {
      if (symbols_.has(name))
      {
        return symbols_[name];
      }
    }
    return nullptr;
  }

  const Element* ElementDB::getElement(UInt atomic_number) const
  {
    if (atomic_numbers_.has(atomic_number))
    {
      return atomic_numbers_[atomic_number];
    }
    return nullptr;
  }

  bool ElementDB::hasElement(const String& name) const
  {
    return names_.has(name) || symbols_.has(name);
  }

  bool ElementDB::hasElement(UInt atomic_number) const
  {
    return atomic_numbers_.has(atomic_number);
  }

  double ElementDB::calculateAvgWeight_(const Map<UInt, double>& Z_to_abundance, const Map<UInt, double>& Z_to_mass)
  {
    double avg = 0;
    // calculate weighted average
    for (Map<UInt, double>::const_iterator it = Z_to_abundance.begin(); it != Z_to_abundance.end(); ++it)
    {
      avg += Z_to_mass[it->first] * Z_to_abundance[it->first];
    }
    return avg;
  }

  double ElementDB::calculateMonoWeight_(const Map<UInt, double>& Z_to_mass)
  {
    double smallest_weight = 1e10;

    for (Map<UInt, double>::const_iterator it = Z_to_mass.begin(); it != Z_to_mass.end(); ++it)
    {
      if (it->second < smallest_weight)
      {
        smallest_weight = it->second;
      }
    }

    return smallest_weight;
  }

  void ElementDB::storeElements()
  {	

    Map<UInt, double> hydrogen_abundance;
    hydrogen_abundance[UInt(1)] = 0.999885;
    hydrogen_abundance[UInt(2)] = 1.15e-04;
    hydrogen_abundance[UInt(3)] = 0.0;

    Map<UInt, double> hydrogen_mass;
    hydrogen_mass[UInt(1)] = 1.0078250319;
    hydrogen_mass[UInt(2)] = 2.01410178;
    hydrogen_mass[UInt(3)] = 3.01604927;

    IsotopeDistribution hydrogen_isotopes = parseIsotopeDistribution_(hydrogen_abundance, hydrogen_mass);
    double hydrogen_avg_weight = calculateAvgWeight_(hydrogen_abundance, hydrogen_mass);
    double hydrogen_mono_weight = calculateMonoWeight_(hydrogen_mass);
    String hydrogen_name = "Hydrogen";
    String hydrogen_symbol = "H";
    UInt hydrogen_an = UInt(1);

    Element* hydrogen = new Element(hydrogen_name, hydrogen_symbol, hydrogen_an, hydrogen_avg_weight, hydrogen_mono_weight, hydrogen_isotopes);
    addElementToMaps(hydrogen_name, hydrogen_symbol, hydrogen_an, hydrogen);
    storeIsotopes(hydrogen_name, hydrogen_symbol, hydrogen_an, hydrogen_mass, hydrogen_isotopes);


    Map<UInt, double> helium_abundance;
    helium_abundance[UInt(3)] = 1.34e-06;
    helium_abundance[UInt(4)] = 0.99999866;

    Map<UInt, double> helium_mass;
    helium_mass[UInt(3)] = 3.0160293191;
    helium_mass[UInt(4)] = 4.00260325415;

    IsotopeDistribution helium_isotopes = parseIsotopeDistribution_(helium_abundance, helium_mass);
    double helium_avg_weight = calculateAvgWeight_(helium_abundance, helium_mass);
    double helium_mono_weight = calculateMonoWeight_(helium_mass);
    String helium_name = "Helium";
    String helium_symbol = "He";
    UInt helium_an = UInt(2);

    Element* helium = new Element(helium_name, helium_symbol, helium_an, helium_avg_weight, helium_mono_weight, helium_isotopes);
    addElementToMaps(helium_name, helium_symbol, helium_an, helium);
    storeIsotopes(helium_name, helium_symbol, helium_an, helium_mass, helium_isotopes);


    Map<UInt, double> lithium_abundance;
    lithium_abundance[UInt(6)] = 0.0759;
    lithium_abundance[UInt(7)] = 0.9241;

    Map<UInt, double> lithium_mass;
    lithium_mass[UInt(6)] = 6.015122;
    lithium_mass[UInt(7)] = 7.016004;

    IsotopeDistribution lithium_isotopes = parseIsotopeDistribution_(lithium_abundance, lithium_mass);
    double lithium_avg_weight = calculateAvgWeight_(lithium_abundance, lithium_mass);
    double lithium_mono_weight = calculateMonoWeight_(lithium_mass);
    String lithium_name = "Lithium";
    String lithium_symbol = "Li";
    UInt lithium_an = UInt(3);

    Element* lithium = new Element(lithium_name, lithium_symbol, lithium_an, lithium_avg_weight, lithium_mono_weight, lithium_isotopes);
    addElementToMaps(lithium_name, lithium_symbol, lithium_an, lithium);
    storeIsotopes(lithium_name, lithium_symbol, lithium_an, lithium_mass, lithium_isotopes);


    Map<UInt, double> beryllium_abundance;
    beryllium_abundance[UInt(9)] =  1.0;

    Map<UInt, double> beryllium_mass;
    beryllium_mass[UInt(9)] =  9.0121822;

    IsotopeDistribution beryllium_isotopes = parseIsotopeDistribution_(beryllium_abundance, beryllium_mass);
    double beryllium_avg_weight = calculateAvgWeight_(beryllium_abundance, beryllium_mass);
    double beryllium_mono_weight = calculateMonoWeight_(beryllium_mass);
    String beryllium_name = "Beryllium";
    String beryllium_symbol = "Be";
    UInt beryllium_an = UInt(4);

    Element* beryllium = new Element(beryllium_name, beryllium_symbol, beryllium_an, beryllium_avg_weight, beryllium_mono_weight, beryllium_isotopes);
    addElementToMaps(beryllium_name, beryllium_symbol, beryllium_an, beryllium);
    storeIsotopes(beryllium_name, beryllium_symbol, beryllium_an, beryllium_mass, beryllium_isotopes);

    Map<UInt, double> bor_abundance;
    bor_abundance[UInt(10)] =  0.199;
    bor_abundance[UInt(11)] =  0.801;

    Map<UInt, double> bor_mass;
    bor_mass[UInt(10)] =  10.012937000000001;
    bor_mass[UInt(11)] =  11.009304999999999;

    IsotopeDistribution bor_isotopes = parseIsotopeDistribution_(bor_abundance, bor_mass);
    double bor_avg_weight = calculateAvgWeight_(bor_abundance, bor_mass);
    double bor_mono_weight = calculateMonoWeight_(bor_mass);
    String bor_name = "Bor";
    String bor_symbol = "B";
    UInt bor_an = UInt(5);

    Element* bor = new Element(bor_name, bor_symbol, bor_an, bor_avg_weight, bor_mono_weight, bor_isotopes);
    addElementToMaps(bor_name, bor_symbol, bor_an, bor);
    storeIsotopes(bor_name, bor_symbol, bor_an, bor_mass, bor_isotopes);


    Map<UInt, double> carbon_abundance;
    carbon_abundance[UInt(12)] =  0.9893;
    carbon_abundance[UInt(13)] =  0.0107;

    Map<UInt, double> carbon_mass;
    carbon_mass[UInt(12)] =  12.0;
    carbon_mass[UInt(13)] =  13.003355000000001;

    IsotopeDistribution carbon_isotopes = parseIsotopeDistribution_(carbon_abundance, carbon_mass);
    double carbon_avg_weight = calculateAvgWeight_(carbon_abundance, carbon_mass);
    double carbon_mono_weight = calculateMonoWeight_(carbon_mass);
    String carbon_name = "Carbon";
    String carbon_symbol = "C";
    UInt carbon_an = UInt(6);

    Element* carbon = new Element(carbon_name, carbon_symbol, carbon_an, carbon_avg_weight, carbon_mono_weight, carbon_isotopes);
    addElementToMaps(carbon_name, carbon_symbol, carbon_an, carbon);
    storeIsotopes(carbon_name, carbon_symbol, carbon_an, carbon_mass, carbon_isotopes);


    Map<UInt, double> nitrogen_abundance;
    nitrogen_abundance[UInt(14)] =  0.99632;
    nitrogen_abundance[UInt(15)] =  3.68e-03;

    Map<UInt, double> nitrogen_mass;
    nitrogen_mass[UInt(14)] =  14.003074;
    nitrogen_mass[UInt(15)] =  15.000109;

    IsotopeDistribution nitrogen_isotopes = parseIsotopeDistribution_(nitrogen_abundance, nitrogen_mass);
    double nitrogen_avg_weight = calculateAvgWeight_(nitrogen_abundance, nitrogen_mass);
    double nitrogen_mono_weight = calculateMonoWeight_(nitrogen_mass);
    String nitrogen_name = "Nitrogen";
    String nitrogen_symbol = "N";
    UInt nitrogen_an = UInt(7);

    Element* nitrogen = new Element(nitrogen_name, nitrogen_symbol, nitrogen_an, nitrogen_avg_weight, nitrogen_mono_weight, nitrogen_isotopes);
    addElementToMaps(nitrogen_name, nitrogen_symbol, nitrogen_an, nitrogen);
    storeIsotopes(nitrogen_name, nitrogen_symbol, nitrogen_an, nitrogen_mass, nitrogen_isotopes);


    Map<UInt, double> oxygen_abundance;
    oxygen_abundance[UInt(16)] =  0.99757;
    oxygen_abundance[UInt(17)] =  3.8e-04;
    oxygen_abundance[UInt(18)] =  2.05e-03;

    Map<UInt, double> oxygen_mass;
    oxygen_mass[UInt(16)] =  15.994915000000001;
    oxygen_mass[UInt(17)] =  16.999132;
    oxygen_mass[UInt(18)] =  17.999168999999998;

    IsotopeDistribution oxygen_isotopes = parseIsotopeDistribution_(oxygen_abundance, oxygen_mass);
    double oxygen_avg_weight = calculateAvgWeight_(oxygen_abundance, oxygen_mass);
    double oxygen_mono_weight = calculateMonoWeight_(oxygen_mass);
    String oxygen_name = "Oxygen";
    String oxygen_symbol = "O";
    UInt oxygen_an = UInt(8);

    Element* oxygen = new Element(oxygen_name, oxygen_symbol, oxygen_an, oxygen_avg_weight, oxygen_mono_weight, oxygen_isotopes);
    addElementToMaps(oxygen_name, oxygen_symbol, oxygen_an, oxygen);
    storeIsotopes(oxygen_name, oxygen_symbol, oxygen_an, oxygen_mass, oxygen_isotopes);


    Map<UInt, double> fluorine_abundance;
    fluorine_abundance[UInt(19)] =  1.0;

    Map<UInt, double> fluorine_mass;
    fluorine_mass[UInt(19)] =  18.99840322;

    IsotopeDistribution fluorine_isotopes = parseIsotopeDistribution_(fluorine_abundance, fluorine_mass);
    double fluorine_avg_weight = calculateAvgWeight_(fluorine_abundance, fluorine_mass);
    double fluorine_mono_weight = calculateMonoWeight_(fluorine_mass);
    String fluorine_name = "Fluorine";
    String fluorine_symbol = "F";
    UInt fluorine_an = UInt(9);

    Element* fluorine = new Element(fluorine_name, fluorine_symbol, fluorine_an, fluorine_avg_weight, fluorine_mono_weight, fluorine_isotopes);
    addElementToMaps(fluorine_name, fluorine_symbol, fluorine_an, fluorine);
    storeIsotopes(fluorine_name, fluorine_symbol, fluorine_an, fluorine_mass, fluorine_isotopes);


    Map<UInt, double> argon_abundance;
    argon_abundance[UInt(36)] =  3.336e-03;
    argon_abundance[UInt(38)] =  6.29e-04;
    argon_abundance[UInt(40)] =  0.996035;

    Map<UInt, double> argon_mass;
    argon_mass[UInt(36)] =  35.967545106000003;
    argon_mass[UInt(38)] =  37.9627324;
    argon_mass[UInt(40)] =  39.9623831225;

    IsotopeDistribution argon_isotopes = parseIsotopeDistribution_(argon_abundance, argon_mass);
    double argon_avg_weight = calculateAvgWeight_(argon_abundance, argon_mass);
    double argon_mono_weight = calculateMonoWeight_(argon_mass);
    String argon_name = "Argon";
    String argon_symbol = "Ar";
    UInt argon_an = UInt(18);

    Element* argon = new Element(argon_name, argon_symbol, argon_an, argon_avg_weight, argon_mono_weight, argon_isotopes);
    addElementToMaps(argon_name, argon_symbol, argon_an, argon);
    storeIsotopes(argon_name, argon_symbol, argon_an, argon_mass, argon_isotopes);


    Map<UInt, double> titanium_abundance;
    titanium_abundance[UInt(46)] =  0.0825;
    titanium_abundance[UInt(47)] =  0.0744;
    titanium_abundance[UInt(48)] =  0.7372;
    titanium_abundance[UInt(49)] =  0.0541;
    titanium_abundance[UInt(50)] =  0.0518;

    Map<UInt, double> titanium_mass;
    titanium_mass[UInt(46)] =  45.952631599999997;
    titanium_mass[UInt(47)] =  46.951763100000001;
    titanium_mass[UInt(48)] =  47.947946299999998;
    titanium_mass[UInt(49)] =  48.947870000000002;
    titanium_mass[UInt(50)] =  49.944791199999997;

    IsotopeDistribution titanium_isotopes = parseIsotopeDistribution_(titanium_abundance, titanium_mass);
    double titanium_avg_weight = calculateAvgWeight_(titanium_abundance, titanium_mass);
    double titanium_mono_weight = calculateMonoWeight_(titanium_mass);
    String titanium_name = "Titanium";
    String titanium_symbol = "Ti";
    UInt titanium_an = UInt(22);

    Element* titanium = new Element(titanium_name, titanium_symbol, titanium_an, titanium_avg_weight, titanium_mono_weight, titanium_isotopes);
    addElementToMaps(titanium_name, titanium_symbol, titanium_an, titanium);
    storeIsotopes(titanium_name, titanium_symbol, titanium_an, titanium_mass, titanium_isotopes);


    Map<UInt, double> sodium_abundance;
    sodium_abundance[UInt(23)] =  1.0;

    Map<UInt, double> sodium_mass;
    sodium_mass[UInt(23)] =  22.989769280899999;

    IsotopeDistribution sodium_isotopes = parseIsotopeDistribution_(sodium_abundance, sodium_mass);
    double sodium_avg_weight = calculateAvgWeight_(sodium_abundance, sodium_mass);
    double sodium_mono_weight = calculateMonoWeight_(sodium_mass);
    String sodium_name = "Sodium";
    String sodium_symbol = "Na";
    UInt sodium_an = UInt(11);

    Element* sodium = new Element(sodium_name, sodium_symbol, sodium_an, sodium_avg_weight, sodium_mono_weight, sodium_isotopes);
    addElementToMaps(sodium_name, sodium_symbol, sodium_an, sodium);
    storeIsotopes(sodium_name, sodium_symbol, sodium_an, sodium_mass, sodium_isotopes);


    Map<UInt, double> magnesium_abundance;
    magnesium_abundance[UInt(24)] =  0.7899;
    magnesium_abundance[UInt(25)] =  0.1;
    magnesium_abundance[UInt(26)] =  0.1101;

    Map<UInt, double> magnesium_mass;
    magnesium_mass[UInt(24)] =  23.985042;
    magnesium_mass[UInt(25)] =  24.985837;
    magnesium_mass[UInt(26)] =  25.982593000000001;

    IsotopeDistribution magnesium_isotopes = parseIsotopeDistribution_(magnesium_abundance, magnesium_mass);
    double magnesium_avg_weight = calculateAvgWeight_(magnesium_abundance, magnesium_mass);
    double magnesium_mono_weight = calculateMonoWeight_(magnesium_mass);
    String magnesium_name = "Magnesium";
    String magnesium_symbol = "Mg";
    UInt magnesium_an = UInt(12);

    Element* magnesium = new Element(magnesium_name, magnesium_symbol, magnesium_an, magnesium_avg_weight, magnesium_mono_weight, magnesium_isotopes);
    addElementToMaps(magnesium_name, magnesium_symbol, magnesium_an, magnesium);
    storeIsotopes(magnesium_name, magnesium_symbol, magnesium_an, magnesium_mass, magnesium_isotopes);


    Map<UInt, double> aluminium_abundance;
    aluminium_abundance[UInt(27)] =  1.0;

    Map<UInt, double> aluminium_mass;
    aluminium_mass[UInt(27)] =  26.981538629999999;

    IsotopeDistribution aluminium_isotopes = parseIsotopeDistribution_(aluminium_abundance, aluminium_mass);
    double aluminium_avg_weight = calculateAvgWeight_(aluminium_abundance, aluminium_mass);
    double aluminium_mono_weight = calculateMonoWeight_(aluminium_mass);
    String aluminium_name = "Aluminium";
    String aluminium_symbol = "Al";
    UInt aluminium_an = UInt(13);

    Element* aluminium = new Element(aluminium_name, aluminium_symbol, aluminium_an, aluminium_avg_weight, aluminium_mono_weight, aluminium_isotopes);
    addElementToMaps(aluminium_name, aluminium_symbol, aluminium_an, aluminium);
    storeIsotopes(aluminium_name, aluminium_symbol, aluminium_an, aluminium_mass, aluminium_isotopes);


    Map<UInt, double> silicon_abundance;
    silicon_abundance[UInt(28)] =  0.9221;
    silicon_abundance[UInt(29)] =  0.0467;
    silicon_abundance[UInt(30)] =  0.031;

    Map<UInt, double> silicon_mass;
    silicon_mass[UInt(28)] =  27.976926532499999;
    silicon_mass[UInt(29)] =  28.9764947;
    silicon_mass[UInt(30)] =  29.973770170000002;

    IsotopeDistribution silicon_isotopes = parseIsotopeDistribution_(silicon_abundance, silicon_mass);
    double silicon_avg_weight = calculateAvgWeight_(silicon_abundance, silicon_mass);
    double silicon_mono_weight = calculateMonoWeight_(silicon_mass);
    String silicon_name = "Silicon";
    String silicon_symbol = "Si";
    UInt silicon_an = UInt(14);

    Element* silicon = new Element(silicon_name, silicon_symbol, silicon_an, silicon_avg_weight, silicon_mono_weight, silicon_isotopes);
    addElementToMaps(silicon_name, silicon_symbol, silicon_an, silicon);
    storeIsotopes(silicon_name, silicon_symbol, silicon_an, silicon_mass, silicon_isotopes);


    Map<UInt, double> phosphorus_abundance;
    phosphorus_abundance[UInt(31)] =  1.0;

    Map<UInt, double> phosphorus_mass;
    phosphorus_mass[UInt(31)] =  30.973761490000001;

    IsotopeDistribution phosphorus_isotopes = parseIsotopeDistribution_(phosphorus_abundance, phosphorus_mass);
    double phosphorus_avg_weight = calculateAvgWeight_(phosphorus_abundance, phosphorus_mass);
    double phosphorus_mono_weight = calculateMonoWeight_(phosphorus_mass);
    String phosphorus_name = "Phosphorus";
    String phosphorus_symbol = "P";
    UInt phosphorus_an = UInt(15);

    Element* phosphorus = new Element(phosphorus_name, phosphorus_symbol, phosphorus_an, phosphorus_avg_weight, phosphorus_mono_weight, phosphorus_isotopes);
    addElementToMaps(phosphorus_name, phosphorus_symbol, phosphorus_an, phosphorus);
    storeIsotopes(phosphorus_name, phosphorus_symbol, phosphorus_an, phosphorus_mass, phosphorus_isotopes);


    Map<UInt, double> sulfur_abundance;
    sulfur_abundance[UInt(32)] =  0.9493;
    sulfur_abundance[UInt(33)] =  0.0076;
    sulfur_abundance[UInt(34)] =  0.0429;
    sulfur_abundance[UInt(36)] =  0.0002;

    Map<UInt, double> sulfur_mass;
    sulfur_mass[UInt(32)] =  31.972070729999999;
    sulfur_mass[UInt(33)] =  32.971457999999998;
    sulfur_mass[UInt(34)] =  33.967866999999998;
    sulfur_mass[UInt(36)] =  35.967081;

    IsotopeDistribution sulfur_isotopes = parseIsotopeDistribution_(sulfur_abundance, sulfur_mass);
    double sulfur_avg_weight = calculateAvgWeight_(sulfur_abundance, sulfur_mass);
    double sulfur_mono_weight = calculateMonoWeight_(sulfur_mass);
    String sulfur_name = "Sulfur";
    String sulfur_symbol = "S";
    UInt sulfur_an = UInt(16);

    Element* sulfur = new Element(sulfur_name, sulfur_symbol, sulfur_an, sulfur_avg_weight, sulfur_mono_weight, sulfur_isotopes);
    addElementToMaps(sulfur_name, sulfur_symbol, sulfur_an, sulfur);
    storeIsotopes(sulfur_name, sulfur_symbol, sulfur_an, sulfur_mass, sulfur_isotopes);


    Map<UInt, double> chlorine_abundance;
    chlorine_abundance[UInt(35)] =  0.7576;
    chlorine_abundance[UInt(37)] =  0.24239999999999998;

    Map<UInt, double> chlorine_mass;
    chlorine_mass[UInt(35)] =  34.968852679999998;
    chlorine_mass[UInt(37)] =  36.965902589999999;

    IsotopeDistribution chlorine_isotopes = parseIsotopeDistribution_(chlorine_abundance, chlorine_mass);
    double chlorine_avg_weight = calculateAvgWeight_(chlorine_abundance, chlorine_mass);
    double chlorine_mono_weight = calculateMonoWeight_(chlorine_mass);
    String chlorine_name = "Chlorine";
    String chlorine_symbol = "Cl";
    UInt chlorine_an = UInt(17);

    Element* chlorine = new Element(chlorine_name, chlorine_symbol, chlorine_an, chlorine_avg_weight, chlorine_mono_weight, chlorine_isotopes);
    addElementToMaps(chlorine_name, chlorine_symbol, chlorine_an, chlorine);
    storeIsotopes(chlorine_name, chlorine_symbol, chlorine_an, chlorine_mass, chlorine_isotopes);


    Map<UInt, double> potassium_abundance;
    potassium_abundance[UInt(39)] =  0.932581;
    potassium_abundance[UInt(40)] =  0.000117;
    potassium_abundance[UInt(41)] =  0.067302;

    Map<UInt, double> potassium_mass;
    potassium_mass[UInt(39)] =  38.963706680000001;
    potassium_mass[UInt(40)] =  39.963998480000001;
    potassium_mass[UInt(41)] =  40.961825760000004;

    IsotopeDistribution potassium_isotopes = parseIsotopeDistribution_(potassium_abundance, potassium_mass);
    double potassium_avg_weight = calculateAvgWeight_(potassium_abundance, potassium_mass);
    double potassium_mono_weight = calculateMonoWeight_(potassium_mass);
    String potassium_name = "Potassium";
    String potassium_symbol = "K";
    UInt potassium_an = UInt(19);

    Element* potassium = new Element(potassium_name, potassium_symbol, potassium_an, potassium_avg_weight, potassium_mono_weight, potassium_isotopes);
    addElementToMaps(potassium_name, potassium_symbol, potassium_an, potassium);
    storeIsotopes(potassium_name, potassium_symbol, potassium_an, potassium_mass, potassium_isotopes);


    Map<UInt, double> calcium_abundance;
    calcium_abundance[UInt(40)] =  0.96941;
    calcium_abundance[UInt(42)] =  0.00647;
    calcium_abundance[UInt(43)] =  0.00135;
    calcium_abundance[UInt(44)] =  0.02086;
    calcium_abundance[UInt(46)] =  4e-05;
    calcium_abundance[UInt(48)] =  0.00187;

    Map<UInt, double> calcium_mass;
    calcium_mass[UInt(40)] =  39.962590980000002;
    calcium_mass[UInt(42)] =  41.958618010000002;
    calcium_mass[UInt(43)] =  42.958766599999997;
    calcium_mass[UInt(44)] =  43.955481800000001;
    calcium_mass[UInt(46)] =  45.953692599999997;
    calcium_mass[UInt(48)] =  47.952534;

    IsotopeDistribution calcium_isotopes = parseIsotopeDistribution_(calcium_abundance, calcium_mass);
    double calcium_avg_weight = calculateAvgWeight_(calcium_abundance, calcium_mass);
    double calcium_mono_weight = calculateMonoWeight_(calcium_mass);
    String calcium_name = "Calcium";
    String calcium_symbol = "Ca";
    UInt calcium_an = UInt(20);

    Element* calcium = new Element(calcium_name, calcium_symbol, calcium_an, calcium_avg_weight, calcium_mono_weight, calcium_isotopes);
    addElementToMaps(calcium_name, calcium_symbol, calcium_an, calcium);
    storeIsotopes(calcium_name, calcium_symbol, calcium_an, calcium_mass, calcium_isotopes);


    Map<UInt, double> scandium_abundance;
    scandium_abundance[UInt(45)] =  1.0;

    Map<UInt, double> scandium_mass;
    scandium_mass[UInt(45)] =  44.955910000000003;

    IsotopeDistribution scandium_isotopes = parseIsotopeDistribution_(scandium_abundance, scandium_mass);
    double scandium_avg_weight = calculateAvgWeight_(scandium_abundance, scandium_mass);
    double scandium_mono_weight = calculateMonoWeight_(scandium_mass);
    String scandium_name = "Scandium";
    String scandium_symbol = "Sc";
    UInt scandium_an = UInt(21);

    Element* scandium = new Element(scandium_name, scandium_symbol, scandium_an, scandium_avg_weight, scandium_mono_weight, scandium_isotopes);
    addElementToMaps(scandium_name, scandium_symbol, scandium_an, scandium);
    storeIsotopes(scandium_name, scandium_symbol, scandium_an, scandium_mass, scandium_isotopes);


    Map<UInt, double> vanadium_abundance;
    vanadium_abundance[UInt(50)] =  0.0025;
    vanadium_abundance[UInt(51)] =  0.9975;

    Map<UInt, double> vanadium_mass;
    vanadium_mass[UInt(50)] =  49.947158500000001;
    vanadium_mass[UInt(51)] =  50.943959499999998;

    IsotopeDistribution vanadium_isotopes = parseIsotopeDistribution_(vanadium_abundance, vanadium_mass);
    double vanadium_avg_weight = calculateAvgWeight_(vanadium_abundance, vanadium_mass);
    double vanadium_mono_weight = calculateMonoWeight_(vanadium_mass);
    String vanadium_name = "Vanadium";
    String vanadium_symbol = "V";
    UInt vanadium_an = UInt(23);

    Element* vanadium = new Element(vanadium_name, vanadium_symbol, vanadium_an, vanadium_avg_weight, vanadium_mono_weight, vanadium_isotopes);
    addElementToMaps(vanadium_name, vanadium_symbol, vanadium_an, vanadium);
    storeIsotopes(vanadium_name, vanadium_symbol, vanadium_an, vanadium_mass, vanadium_isotopes);


    Map<UInt, double> chromium_abundance;
    chromium_abundance[UInt(50)] =  0.043449999999999996;
    chromium_abundance[UInt(52)] =  0.83789;
    chromium_abundance[UInt(53)] =  0.09501;
    chromium_abundance[UInt(54)] =  0.02365;

    Map<UInt, double> chromium_mass;
    chromium_mass[UInt(50)] =  49.946044200000003;
    chromium_mass[UInt(52)] =  51.940507500000003;
    chromium_mass[UInt(53)] =  52.940649399999998;
    chromium_mass[UInt(54)] =  53.938880400000002;

    IsotopeDistribution chromium_isotopes = parseIsotopeDistribution_(chromium_abundance, chromium_mass);
    double chromium_avg_weight = calculateAvgWeight_(chromium_abundance, chromium_mass);
    double chromium_mono_weight = calculateMonoWeight_(chromium_mass);
    String chromium_name = "Chromium";
    String chromium_symbol = "Cr";
    UInt chromium_an = UInt(24);

    Element* chromium = new Element(chromium_name, chromium_symbol, chromium_an, chromium_avg_weight, chromium_mono_weight, chromium_isotopes);
    addElementToMaps(chromium_name, chromium_symbol, chromium_an, chromium);
    storeIsotopes(chromium_name, chromium_symbol, chromium_an, chromium_mass, chromium_isotopes);


    Map<UInt, double> tellurium_abundance;
    tellurium_abundance[UInt(120)] =  0.0009;
    tellurium_abundance[UInt(122)] =  0.0255;
    tellurium_abundance[UInt(124)] =  0.047400000000000005;
    tellurium_abundance[UInt(125)] =  0.0707;
    tellurium_abundance[UInt(126)] =  0.1884;
    tellurium_abundance[UInt(128)] =  0.31739999999999996;
    tellurium_abundance[UInt(130)] =  0.3408;

    Map<UInt, double> tellurium_mass;
    tellurium_mass[UInt(120)] =  119.904020000000003;
    tellurium_mass[UInt(122)] =  121.9030439;
    tellurium_mass[UInt(124)] =  123.902817900000002;
    tellurium_mass[UInt(125)] =  124.904430700000006;
    tellurium_mass[UInt(126)] =  125.903311700000003;
    tellurium_mass[UInt(128)] =  127.904463100000001;
    tellurium_mass[UInt(130)] =  129.906224400000014;

    IsotopeDistribution tellurium_isotopes = parseIsotopeDistribution_(tellurium_abundance, tellurium_mass);
    double tellurium_avg_weight = calculateAvgWeight_(tellurium_abundance, tellurium_mass);
    double tellurium_mono_weight = calculateMonoWeight_(tellurium_mass);
    String tellurium_name = "Tellurium";
    String tellurium_symbol = "Te";
    UInt tellurium_an = UInt(52);

    Element* tellurium = new Element(tellurium_name, tellurium_symbol, tellurium_an, tellurium_avg_weight, tellurium_mono_weight, tellurium_isotopes);
    addElementToMaps(tellurium_name, tellurium_symbol, tellurium_an, tellurium);
    storeIsotopes(tellurium_name, tellurium_symbol, tellurium_an, tellurium_mass, tellurium_isotopes);


    Map<UInt, double> barium_abundance;
    barium_abundance[UInt(132)] =  0.00101;
    barium_abundance[UInt(134)] =  0.024169999999999997;
    barium_abundance[UInt(135)] =  0.06591999999999999;
    barium_abundance[UInt(136)] =  0.07854;
    barium_abundance[UInt(137)] =  0.11231999999999999;
    barium_abundance[UInt(138)] =  0.71698;

    Map<UInt, double> barium_mass;
    barium_mass[UInt(132)] =  131.9050613;
    barium_mass[UInt(134)] =  133.904508399999997;
    barium_mass[UInt(135)] =  134.905688599999991;
    barium_mass[UInt(136)] =  135.904575899999998;
    barium_mass[UInt(137)] =  136.905827399999993;
    barium_mass[UInt(138)] =  137.905247199999991;

    IsotopeDistribution barium_isotopes = parseIsotopeDistribution_(barium_abundance, barium_mass);
    double barium_avg_weight = calculateAvgWeight_(barium_abundance, barium_mass);
    double barium_mono_weight = calculateMonoWeight_(barium_mass);
    String barium_name = "Barium";
    String barium_symbol = "Ba";
    UInt barium_an = UInt(56);

    Element* barium = new Element(barium_name, barium_symbol, barium_an, barium_avg_weight, barium_mono_weight, barium_isotopes);
    addElementToMaps(barium_name, barium_symbol, barium_an, barium);
    storeIsotopes(barium_name, barium_symbol, barium_an, barium_mass, barium_isotopes);


    Map<UInt, double> manganese_abundance;
    manganese_abundance[UInt(55)] =  1.0;

    Map<UInt, double> manganese_mass;
    manganese_mass[UInt(55)] =  54.938049999999997;

    IsotopeDistribution manganese_isotopes = parseIsotopeDistribution_(manganese_abundance, manganese_mass);
    double manganese_avg_weight = calculateAvgWeight_(manganese_abundance, manganese_mass);
    double manganese_mono_weight = calculateMonoWeight_(manganese_mass);
    String manganese_name = "Manganese";
    String manganese_symbol = "Mn";
    UInt manganese_an = UInt(25);

    Element* manganese = new Element(manganese_name, manganese_symbol, manganese_an, manganese_avg_weight, manganese_mono_weight, manganese_isotopes);
    addElementToMaps(manganese_name, manganese_symbol, manganese_an, manganese);
    storeIsotopes(manganese_name, manganese_symbol, manganese_an, manganese_mass, manganese_isotopes);


    Map<UInt, double> ferrum_abundance;
    ferrum_abundance[UInt(54)] =  0.058449999999999995;
    ferrum_abundance[UInt(56)] =  0.91754;
    ferrum_abundance[UInt(57)] =  0.021191;
    ferrum_abundance[UInt(58)] =  0.002819;

    Map<UInt, double> ferrum_mass;
    ferrum_mass[UInt(54)] =  53.939610500000001;
    ferrum_mass[UInt(56)] =  55.934937499999997;
    ferrum_mass[UInt(57)] =  56.935394000000002;
    ferrum_mass[UInt(58)] =  57.933275600000002;

    IsotopeDistribution ferrum_isotopes = parseIsotopeDistribution_(ferrum_abundance, ferrum_mass);
    double ferrum_avg_weight = calculateAvgWeight_(ferrum_abundance, ferrum_mass);
    double ferrum_mono_weight = calculateMonoWeight_(ferrum_mass);
    String ferrum_name = "Ferrum";
    String ferrum_symbol = "Fe";
    UInt ferrum_an = UInt(26);

    Element* ferrum = new Element(ferrum_name, ferrum_symbol, ferrum_an, ferrum_avg_weight, ferrum_mono_weight, ferrum_isotopes);
    addElementToMaps(ferrum_name, ferrum_symbol, ferrum_an, ferrum);
    storeIsotopes(ferrum_name, ferrum_symbol, ferrum_an, ferrum_mass, ferrum_isotopes);


    Map<UInt, double> cobalt_abundance;
    cobalt_abundance[UInt(59)] =  1.0;

    Map<UInt, double> cobalt_mass;
    cobalt_mass[UInt(59)] =  58.933194999999998;

    IsotopeDistribution cobalt_isotopes = parseIsotopeDistribution_(cobalt_abundance, cobalt_mass);
    double cobalt_avg_weight = calculateAvgWeight_(cobalt_abundance, cobalt_mass);
    double cobalt_mono_weight = calculateMonoWeight_(cobalt_mass);
    String cobalt_name = "Cobalt";
    String cobalt_symbol = "Co";
    UInt cobalt_an = UInt(27);

    Element* cobalt = new Element(cobalt_name, cobalt_symbol, cobalt_an, cobalt_avg_weight, cobalt_mono_weight, cobalt_isotopes);
    addElementToMaps(cobalt_name, cobalt_symbol, cobalt_an, cobalt);
    storeIsotopes(cobalt_name, cobalt_symbol, cobalt_an, cobalt_mass, cobalt_isotopes);


    Map<UInt, double> nickel_abundance;
    nickel_abundance[UInt(58)] =  0.680169;
    nickel_abundance[UInt(60)] =  0.262231;
    nickel_abundance[UInt(61)] =  0.011399;
    nickel_abundance[UInt(62)] =  0.036345;
    nickel_abundance[UInt(64)] =  0.009256;

    Map<UInt, double> nickel_mass;
    nickel_mass[UInt(58)] =  57.935347999999998;
    nickel_mass[UInt(60)] =  59.930790999999999;
    nickel_mass[UInt(61)] =  60.931060000000002;
    nickel_mass[UInt(62)] =  61.928348999999997;
    nickel_mass[UInt(64)] =  63.927970000000002;

    IsotopeDistribution nickel_isotopes = parseIsotopeDistribution_(nickel_abundance, nickel_mass);
    double nickel_avg_weight = calculateAvgWeight_(nickel_abundance, nickel_mass);
    double nickel_mono_weight = calculateMonoWeight_(nickel_mass);
    String nickel_name = "Nickel";
    String nickel_symbol = "Ni";
    UInt nickel_an = UInt(28);

    Element* nickel = new Element(nickel_name, nickel_symbol, nickel_an, nickel_avg_weight, nickel_mono_weight, nickel_isotopes);
    addElementToMaps(nickel_name, nickel_symbol, nickel_an, nickel);
    storeIsotopes(nickel_name, nickel_symbol, nickel_an, nickel_mass, nickel_isotopes);


    Map<UInt, double> copper_abundance;
    copper_abundance[UInt(63)] =  0.6917;
    copper_abundance[UInt(65)] =  0.30829999999999996;

    Map<UInt, double> copper_mass;
    copper_mass[UInt(63)] =  62.929600999999998;
    copper_mass[UInt(65)] =  64.927794000000006;

    IsotopeDistribution copper_isotopes = parseIsotopeDistribution_(copper_abundance, copper_mass);
    double copper_avg_weight = calculateAvgWeight_(copper_abundance, copper_mass);
    double copper_mono_weight = calculateMonoWeight_(copper_mass);
    String copper_name = "Copper";
    String copper_symbol = "Cu";
    UInt copper_an = UInt(29);

    Element* copper = new Element(copper_name, copper_symbol, copper_an, copper_avg_weight, copper_mono_weight, copper_isotopes);
    addElementToMaps(copper_name, copper_symbol, copper_an, copper);
    storeIsotopes(copper_name, copper_symbol, copper_an, copper_mass, copper_isotopes);


    Map<UInt, double> zinc_abundance;
    zinc_abundance[UInt(64)] =  0.4863;
    zinc_abundance[UInt(66)] =  0.27899999999999997;
    zinc_abundance[UInt(67)] =  0.040999999999999995;
    zinc_abundance[UInt(68)] =  0.1875;
    zinc_abundance[UInt(70)] =  0.0062;

    Map<UInt, double> zinc_mass;
    zinc_mass[UInt(64)] =  63.929147;
    zinc_mass[UInt(66)] =  65.926036999999994;
    zinc_mass[UInt(67)] =  66.927131000000003;
    zinc_mass[UInt(68)] =  67.924847999999997;
    zinc_mass[UInt(70)] =  69.925325000000001;

    IsotopeDistribution zinc_isotopes = parseIsotopeDistribution_(zinc_abundance, zinc_mass);
    double zinc_avg_weight = calculateAvgWeight_(zinc_abundance, zinc_mass);
    double zinc_mono_weight = calculateMonoWeight_(zinc_mass);
    String zinc_name = "Zinc";
    String zinc_symbol = "Zn";
    UInt zinc_an = UInt(30);

    Element* zinc = new Element(zinc_name, zinc_symbol, zinc_an, zinc_avg_weight, zinc_mono_weight, zinc_isotopes);
    addElementToMaps(zinc_name, zinc_symbol, zinc_an, zinc);
    storeIsotopes(zinc_name, zinc_symbol, zinc_an, zinc_mass, zinc_isotopes);


    Map<UInt, double> gallium_abundance;
    gallium_abundance[UInt(69)] =  0.60108;
    gallium_abundance[UInt(71)] =  0.39892000000000005;

    Map<UInt, double> gallium_mass;
    gallium_mass[UInt(69)] =  68.925573600000007;
    gallium_mass[UInt(71)] =  70.924701299999995;

    IsotopeDistribution gallium_isotopes = parseIsotopeDistribution_(gallium_abundance, gallium_mass);
    double gallium_avg_weight = calculateAvgWeight_(gallium_abundance, gallium_mass);
    double gallium_mono_weight = calculateMonoWeight_(gallium_mass);
    String gallium_name = "Gallium";
    String gallium_symbol = "Ga";
    UInt gallium_an = UInt(31);

    Element* gallium = new Element(gallium_name, gallium_symbol, gallium_an, gallium_avg_weight, gallium_mono_weight, gallium_isotopes);
    addElementToMaps(gallium_name, gallium_symbol, gallium_an, gallium);
    storeIsotopes(gallium_name, gallium_symbol, gallium_an, gallium_mass, gallium_isotopes);


    Map<UInt, double> germanium_abundance;
    germanium_abundance[UInt(70)] =  0.20379999999999998;
    germanium_abundance[UInt(72)] =  0.2731;
    germanium_abundance[UInt(73)] =  0.0776;
    germanium_abundance[UInt(74)] =  0.36719999999999997;
    germanium_abundance[UInt(76)] =  0.0776;

    Map<UInt, double> germanium_mass;
    germanium_mass[UInt(70)] =  69.924247399999999;
    germanium_mass[UInt(72)] =  71.922075800000002;
    germanium_mass[UInt(73)] =  72.9234589;
    germanium_mass[UInt(74)] =  73.921177799999995;
    germanium_mass[UInt(76)] =  75.921401;

    IsotopeDistribution germanium_isotopes = parseIsotopeDistribution_(germanium_abundance, germanium_mass);
    double germanium_avg_weight = calculateAvgWeight_(germanium_abundance, germanium_mass);
    double germanium_mono_weight = calculateMonoWeight_(germanium_mass);
    String germanium_name = "Germanium";
    String germanium_symbol = "Ge";
    UInt germanium_an = UInt(32);

    Element* germanium = new Element(germanium_name, germanium_symbol, germanium_an, germanium_avg_weight, germanium_mono_weight, germanium_isotopes);
    addElementToMaps(germanium_name, germanium_symbol, germanium_an, germanium);
    storeIsotopes(germanium_name, germanium_symbol, germanium_an, germanium_mass, germanium_isotopes);


    Map<UInt, double> arsenic_abundance;
    arsenic_abundance[UInt(75)] =  1.0;

    Map<UInt, double> arsenic_mass;
    arsenic_mass[UInt(75)] =  74.921596500000007;

    IsotopeDistribution arsenic_isotopes = parseIsotopeDistribution_(arsenic_abundance, arsenic_mass);
    double arsenic_avg_weight = calculateAvgWeight_(arsenic_abundance, arsenic_mass);
    double arsenic_mono_weight = calculateMonoWeight_(arsenic_mass);
    String arsenic_name = "Arsenic";
    String arsenic_symbol = "As";
    UInt arsenic_an = UInt(33);

    Element* arsenic = new Element(arsenic_name, arsenic_symbol, arsenic_an, arsenic_avg_weight, arsenic_mono_weight, arsenic_isotopes);
    addElementToMaps(arsenic_name, arsenic_symbol, arsenic_an, arsenic);
    storeIsotopes(arsenic_name, arsenic_symbol, arsenic_an, arsenic_mass, arsenic_isotopes);


    Map<UInt, double> rubidium_abundance;
    rubidium_abundance[UInt(85)] =  0.7217;

    Map<UInt, double> rubidium_mass;
    rubidium_mass[UInt(85)] =  84.911789737999996;

    IsotopeDistribution rubidium_isotopes = parseIsotopeDistribution_(rubidium_abundance, rubidium_mass);
    double rubidium_avg_weight = calculateAvgWeight_(rubidium_abundance, rubidium_mass);
    double rubidium_mono_weight = calculateMonoWeight_(rubidium_mass);
    String rubidium_name = "Rubidium";
    String rubidium_symbol = "Rb";
    UInt rubidium_an = UInt(37);

    Element* rubidium = new Element(rubidium_name, rubidium_symbol, rubidium_an, rubidium_avg_weight, rubidium_mono_weight, rubidium_isotopes);
    addElementToMaps(rubidium_name, rubidium_symbol, rubidium_an, rubidium);
    storeIsotopes(rubidium_name, rubidium_symbol, rubidium_an, rubidium_mass, rubidium_isotopes);


    Map<UInt, double> strontium_abundance;
    strontium_abundance[UInt(84)] =  0.005600000000000001;
    strontium_abundance[UInt(86)] =  0.0986;
    strontium_abundance[UInt(87)] =  0.07;
    strontium_abundance[UInt(88)] =  0.8258;

    Map<UInt, double> strontium_mass;
    strontium_mass[UInt(84)] =  83.913425000000004;
    strontium_mass[UInt(86)] =  85.909260730900002;
    strontium_mass[UInt(87)] =  86.908877497000006;
    strontium_mass[UInt(88)] =  87.905612257100003;

    IsotopeDistribution strontium_isotopes = parseIsotopeDistribution_(strontium_abundance, strontium_mass);
    double strontium_avg_weight = calculateAvgWeight_(strontium_abundance, strontium_mass);
    double strontium_mono_weight = calculateMonoWeight_(strontium_mass);
    String strontium_name = "Strontium";
    String strontium_symbol = "Sr";
    UInt strontium_an = UInt(38);

    Element* strontium = new Element(strontium_name, strontium_symbol, strontium_an, strontium_avg_weight, strontium_mono_weight, strontium_isotopes);
    addElementToMaps(strontium_name, strontium_symbol, strontium_an, strontium);
    storeIsotopes(strontium_name, strontium_symbol, strontium_an, strontium_mass, strontium_isotopes);


    Map<UInt, double> yttrium_abundance;
    yttrium_abundance[UInt(89)] =  1.0;

    Map<UInt, double> yttrium_mass;
    yttrium_mass[UInt(89)] =  88.905850000000001;

    IsotopeDistribution yttrium_isotopes = parseIsotopeDistribution_(yttrium_abundance, yttrium_mass);
    double yttrium_avg_weight = calculateAvgWeight_(yttrium_abundance, yttrium_mass);
    double yttrium_mono_weight = calculateMonoWeight_(yttrium_mass);
    String yttrium_name = "Yttrium";
    String yttrium_symbol = "Y";
    UInt yttrium_an = UInt(39);

    Element* yttrium = new Element(yttrium_name, yttrium_symbol, yttrium_an, yttrium_avg_weight, yttrium_mono_weight, yttrium_isotopes);
    addElementToMaps(yttrium_name, yttrium_symbol, yttrium_an, yttrium);
    storeIsotopes(yttrium_name, yttrium_symbol, yttrium_an, yttrium_mass, yttrium_isotopes);


    Map<UInt, double> zirconium_abundance;
    zirconium_abundance[UInt(90)] =  0.5145000000000001;
    zirconium_abundance[UInt(91)] =  0.11220000000000001;
    zirconium_abundance[UInt(92)] =  0.17149999999999999;
    zirconium_abundance[UInt(94)] =  0.17379999999999998;

    Map<UInt, double> zirconium_mass;
    zirconium_mass[UInt(90)] =  89.9047044;
    zirconium_mass[UInt(91)] =  90.905645800000002;
    zirconium_mass[UInt(92)] =  91.905040799999995;
    zirconium_mass[UInt(94)] =  93.906315199999995;

    IsotopeDistribution zirconium_isotopes = parseIsotopeDistribution_(zirconium_abundance, zirconium_mass);
    double zirconium_avg_weight = calculateAvgWeight_(zirconium_abundance, zirconium_mass);
    double zirconium_mono_weight = calculateMonoWeight_(zirconium_mass);
    String zirconium_name = "Zirconium";
    String zirconium_symbol = "Zr";
    UInt zirconium_an = UInt(40);

    Element* zirconium = new Element(zirconium_name, zirconium_symbol, zirconium_an, zirconium_avg_weight, zirconium_mono_weight, zirconium_isotopes);
    addElementToMaps(zirconium_name, zirconium_symbol, zirconium_an, zirconium);
    storeIsotopes(zirconium_name, zirconium_symbol, zirconium_an, zirconium_mass, zirconium_isotopes);


    Map<UInt, double> nibium_abundance;
    nibium_abundance[UInt(93)] =  1.0;

    Map<UInt, double> nibium_mass;
    nibium_mass[UInt(93)] =  92.906378099999998;

    IsotopeDistribution nibium_isotopes = parseIsotopeDistribution_(nibium_abundance, nibium_mass);
    double nibium_avg_weight = calculateAvgWeight_(nibium_abundance, nibium_mass);
    double nibium_mono_weight = calculateMonoWeight_(nibium_mass);
    String nibium_name = "Nibium";
    String nibium_symbol = "Nb";
    UInt nibium_an = UInt(41);

    Element* nibium = new Element(nibium_name, nibium_symbol, nibium_an, nibium_avg_weight, nibium_mono_weight, nibium_isotopes);
    addElementToMaps(nibium_name, nibium_symbol, nibium_an, nibium);
    storeIsotopes(nibium_name, nibium_symbol, nibium_an, nibium_mass, nibium_isotopes);


    Map<UInt, double> ruthenium_abundance;
    ruthenium_abundance[UInt(96)] =  0.0554;
    ruthenium_abundance[UInt(98)] =  0.0187;
    ruthenium_abundance[UInt(99)] =  0.1276;
    ruthenium_abundance[UInt(100)] =  0.126;
    ruthenium_abundance[UInt(101)] =  0.17059999999999997;
    ruthenium_abundance[UInt(102)] =  0.3155;
    ruthenium_abundance[UInt(104)] =  0.1862;

    Map<UInt, double> ruthenium_mass;
    ruthenium_mass[UInt(96)] =  95.907597999999993;
    ruthenium_mass[UInt(98)] =  97.905287000000001;
    ruthenium_mass[UInt(99)] =  98.9059393;
    ruthenium_mass[UInt(100)] =  99.904219499999996;
    ruthenium_mass[UInt(101)] =  100.905582100000004;
    ruthenium_mass[UInt(102)] =  101.904349300000007;
    ruthenium_mass[UInt(104)] =  103.905433000000002;

    IsotopeDistribution ruthenium_isotopes = parseIsotopeDistribution_(ruthenium_abundance, ruthenium_mass);
    double ruthenium_avg_weight = calculateAvgWeight_(ruthenium_abundance, ruthenium_mass);
    double ruthenium_mono_weight = calculateMonoWeight_(ruthenium_mass);
    String ruthenium_name = "Ruthenium";
    String ruthenium_symbol = "Ru";
    UInt ruthenium_an = UInt(44);

    Element* ruthenium = new Element(ruthenium_name, ruthenium_symbol, ruthenium_an, ruthenium_avg_weight, ruthenium_mono_weight, ruthenium_isotopes);
    addElementToMaps(ruthenium_name, ruthenium_symbol, ruthenium_an, ruthenium);
    storeIsotopes(ruthenium_name, ruthenium_symbol, ruthenium_an, ruthenium_mass, ruthenium_isotopes);


    Map<UInt, double> tin_abundance;
    tin_abundance[UInt(112)] =  0.0097;
    tin_abundance[UInt(114)] =  0.0066;
    tin_abundance[UInt(115)] =  0.0034000000000000002;
    tin_abundance[UInt(116)] =  0.1454;
    tin_abundance[UInt(117)] =  0.0768;
    tin_abundance[UInt(118)] =  0.2422;
    tin_abundance[UInt(119)] =  0.0859;
    tin_abundance[UInt(120)] =  0.3258;
    tin_abundance[UInt(122)] =  0.0463;
    tin_abundance[UInt(124)] =  0.0579;

    Map<UInt, double> tin_mass;
    tin_mass[UInt(112)] =  111.904818000000006;
    tin_mass[UInt(114)] =  113.902777900000004;
    tin_mass[UInt(115)] =  114.903341999999995;
    tin_mass[UInt(116)] =  115.901741000000001;
    tin_mass[UInt(117)] =  116.902951999999999;
    tin_mass[UInt(118)] =  117.901602999999994;
    tin_mass[UInt(119)] =  118.903307999999996;
    tin_mass[UInt(120)] =  119.902194699999996;
    tin_mass[UInt(122)] =  121.903439000000006;
    tin_mass[UInt(124)] =  123.905273899999997;

    IsotopeDistribution tin_isotopes = parseIsotopeDistribution_(tin_abundance, tin_mass);
    double tin_avg_weight = calculateAvgWeight_(tin_abundance, tin_mass);
    double tin_mono_weight = calculateMonoWeight_(tin_mass);
    String tin_name = "Tin";
    String tin_symbol = "Sn";
    UInt tin_an = UInt(50);

    Element* tin = new Element(tin_name, tin_symbol, tin_an, tin_avg_weight, tin_mono_weight, tin_isotopes);
    addElementToMaps(tin_name, tin_symbol, tin_an, tin);
    storeIsotopes(tin_name, tin_symbol, tin_an, tin_mass, tin_isotopes);


    Map<UInt, double> antimony_abundance;
    antimony_abundance[UInt(121)] =  0.5721;
    antimony_abundance[UInt(123)] =  0.4279;

    Map<UInt, double> antimony_mass;
    antimony_mass[UInt(121)] =  120.903815699999996;
    antimony_mass[UInt(123)] =  122.904213999999996;

    IsotopeDistribution antimony_isotopes = parseIsotopeDistribution_(antimony_abundance, antimony_mass);
    double antimony_avg_weight = calculateAvgWeight_(antimony_abundance, antimony_mass);
    double antimony_mono_weight = calculateMonoWeight_(antimony_mass);
    String antimony_name = "Antimony";
    String antimony_symbol = "Sb";
    UInt antimony_an = UInt(51);

    Element* antimony = new Element(antimony_name, antimony_symbol, antimony_an, antimony_avg_weight, antimony_mono_weight, antimony_isotopes);
    addElementToMaps(antimony_name, antimony_symbol, antimony_an, antimony);
    storeIsotopes(antimony_name, antimony_symbol, antimony_an, antimony_mass, antimony_isotopes);


    Map<UInt, double> selenium_abundance;
    selenium_abundance[UInt(74)] =  0.00889;
    selenium_abundance[UInt(76)] =  0.09366;
    selenium_abundance[UInt(77)] =  0.07635;
    selenium_abundance[UInt(78)] =  0.23772;
    selenium_abundance[UInt(80)] =  0.49607;
    selenium_abundance[UInt(82)] =  0.08731;

    Map<UInt, double> selenium_mass;
    selenium_mass[UInt(74)] =  73.922476399999994;
    selenium_mass[UInt(76)] =  75.919213600000006;
    selenium_mass[UInt(77)] =  76.919914000000006;
    selenium_mass[UInt(78)] =  77.917309099999997;
    selenium_mass[UInt(80)] =  79.916521299999999;
    selenium_mass[UInt(82)] =  81.916699399999999;

    IsotopeDistribution selenium_isotopes = parseIsotopeDistribution_(selenium_abundance, selenium_mass);
    double selenium_avg_weight = calculateAvgWeight_(selenium_abundance, selenium_mass);
    double selenium_mono_weight = calculateMonoWeight_(selenium_mass);
    String selenium_name = "Selenium";
    String selenium_symbol = "Se";
    UInt selenium_an = UInt(34);

    Element* selenium = new Element(selenium_name, selenium_symbol, selenium_an, selenium_avg_weight, selenium_mono_weight, selenium_isotopes);
    addElementToMaps(selenium_name, selenium_symbol, selenium_an, selenium);
    storeIsotopes(selenium_name, selenium_symbol, selenium_an, selenium_mass, selenium_isotopes);


    Map<UInt, double> bromine_abundance;
    bromine_abundance[UInt(79)] =  0.5069;
    bromine_abundance[UInt(81)] =  0.49310000000000004;

    Map<UInt, double> bromine_mass;
    bromine_mass[UInt(79)] =  78.918337100000002;
    bromine_mass[UInt(81)] =  80.916290599999996;

    IsotopeDistribution bromine_isotopes = parseIsotopeDistribution_(bromine_abundance, bromine_mass);
    double bromine_avg_weight = calculateAvgWeight_(bromine_abundance, bromine_mass);
    double bromine_mono_weight = calculateMonoWeight_(bromine_mass);
    String bromine_name = "Bromine";
    String bromine_symbol = "Br";
    UInt bromine_an = UInt(35);

    Element* bromine = new Element(bromine_name, bromine_symbol, bromine_an, bromine_avg_weight, bromine_mono_weight, bromine_isotopes);
    addElementToMaps(bromine_name, bromine_symbol, bromine_an, bromine);
    storeIsotopes(bromine_name, bromine_symbol, bromine_an, bromine_mass, bromine_isotopes);


    Map<UInt, double> krypton_abundance;
    krypton_abundance[UInt(78)] =  0.0034999999999999996;
    krypton_abundance[UInt(80)] =  0.0225;
    krypton_abundance[UInt(82)] =  0.11599999999999999;
    krypton_abundance[UInt(83)] =  0.115;
    krypton_abundance[UInt(84)] =  0.57;
    krypton_abundance[UInt(86)] =  0.17300000000000001;

    Map<UInt, double> krypton_mass;
    krypton_mass[UInt(78)] =  77.920400000000001;
    krypton_mass[UInt(80)] =  79.916380000000004;
    krypton_mass[UInt(82)] =  81.913482000000002;
    krypton_mass[UInt(83)] =  82.914135000000002;
    krypton_mass[UInt(84)] =  83.911507;
    krypton_mass[UInt(86)] =  85.910616000000005;

    IsotopeDistribution krypton_isotopes = parseIsotopeDistribution_(krypton_abundance, krypton_mass);
    double krypton_avg_weight = calculateAvgWeight_(krypton_abundance, krypton_mass);
    double krypton_mono_weight = calculateMonoWeight_(krypton_mass);
    String krypton_name = "Krypton";
    String krypton_symbol = "Kr";
    UInt krypton_an = UInt(36);

    Element* krypton = new Element(krypton_name, krypton_symbol, krypton_an, krypton_avg_weight, krypton_mono_weight, krypton_isotopes);
    addElementToMaps(krypton_name, krypton_symbol, krypton_an, krypton);
    storeIsotopes(krypton_name, krypton_symbol, krypton_an, krypton_mass, krypton_isotopes);


    Map<UInt, double> molybdenum_abundance;
    molybdenum_abundance[UInt(92)] =  0.1484;
    molybdenum_abundance[UInt(94)] =  0.0925;
    molybdenum_abundance[UInt(95)] =  0.1592;
    molybdenum_abundance[UInt(96)] =  0.1668;
    molybdenum_abundance[UInt(97)] =  0.0955;
    molybdenum_abundance[UInt(98)] =  0.2413;
    molybdenum_abundance[UInt(100)] =  0.09630000000000001;

    Map<UInt, double> molybdenum_mass;
    molybdenum_mass[UInt(92)] =  91.906809999999993;
    molybdenum_mass[UInt(94)] =  93.905088000000006;
    molybdenum_mass[UInt(95)] =  94.905840999999995;
    molybdenum_mass[UInt(96)] =  95.904679000000002;
    molybdenum_mass[UInt(97)] =  96.906020999999996;
    molybdenum_mass[UInt(98)] =  97.905407999999994;
    molybdenum_mass[UInt(100)] =  99.907477;

    IsotopeDistribution molybdenum_isotopes = parseIsotopeDistribution_(molybdenum_abundance, molybdenum_mass);
    double molybdenum_avg_weight = calculateAvgWeight_(molybdenum_abundance, molybdenum_mass);
    double molybdenum_mono_weight = calculateMonoWeight_(molybdenum_mass);
    String molybdenum_name = "Molybdenum";
    String molybdenum_symbol = "Mo";
    UInt molybdenum_an = UInt(42);

    Element* molybdenum = new Element(molybdenum_name, molybdenum_symbol, molybdenum_an, molybdenum_avg_weight, molybdenum_mono_weight, molybdenum_isotopes);
    addElementToMaps(molybdenum_name, molybdenum_symbol, molybdenum_an, molybdenum);
    storeIsotopes(molybdenum_name, molybdenum_symbol, molybdenum_an, molybdenum_mass, molybdenum_isotopes);


    Map<UInt, double> technitium_abundance;
    technitium_abundance[UInt(97)] =  0.0;
    technitium_abundance[UInt(98)] =  0.0;
    technitium_abundance[UInt(99)] =  0.0;

    Map<UInt, double> technitium_mass;
    technitium_mass[UInt(97)] =  96.906363999999996;
    technitium_mass[UInt(98)] =  97.907214999999994;
    technitium_mass[UInt(99)] =  98.906254000000004;

    IsotopeDistribution technitium_isotopes = parseIsotopeDistribution_(technitium_abundance, technitium_mass);
    double technitium_avg_weight = calculateAvgWeight_(technitium_abundance, technitium_mass);
    double technitium_mono_weight = calculateMonoWeight_(technitium_mass);
    String technitium_name = "Technitium";
    String technitium_symbol = "Tc";
    UInt technitium_an = UInt(43);

    Element* technitium = new Element(technitium_name, technitium_symbol, technitium_an, technitium_avg_weight, technitium_mono_weight, technitium_isotopes);
    addElementToMaps(technitium_name, technitium_symbol, technitium_an, technitium);
    storeIsotopes(technitium_name, technitium_symbol, technitium_an, technitium_mass, technitium_isotopes);


    Map<UInt, double> rhodium_abundance;
    rhodium_abundance[UInt(103)] =  1.0;

    Map<UInt, double> rhodium_mass;
    rhodium_mass[UInt(103)] =  102.905500000000004;

    IsotopeDistribution rhodium_isotopes = parseIsotopeDistribution_(rhodium_abundance, rhodium_mass);
    double rhodium_avg_weight = calculateAvgWeight_(rhodium_abundance, rhodium_mass);
    double rhodium_mono_weight = calculateMonoWeight_(rhodium_mass);
    String rhodium_name = "Rhodium";
    String rhodium_symbol = "Rh";
    UInt rhodium_an = UInt(45);

    Element* rhodium = new Element(rhodium_name, rhodium_symbol, rhodium_an, rhodium_avg_weight, rhodium_mono_weight, rhodium_isotopes);
    addElementToMaps(rhodium_name, rhodium_symbol, rhodium_an, rhodium);
    storeIsotopes(rhodium_name, rhodium_symbol, rhodium_an, rhodium_mass, rhodium_isotopes);


    Map<UInt, double> palladium_abundance;
    palladium_abundance[UInt(102)] =  0.0102;
    palladium_abundance[UInt(104)] =  0.1114;
    palladium_abundance[UInt(105)] =  0.22329999999999997;
    palladium_abundance[UInt(106)] =  0.2733;
    palladium_abundance[UInt(108)] =  0.2646;
    palladium_abundance[UInt(110)] =  0.11720000000000001;

    Map<UInt, double> palladium_mass;
    palladium_mass[UInt(102)] =  101.905608999999998;
    palladium_mass[UInt(104)] =  103.904036000000005;
    palladium_mass[UInt(105)] =  104.905085;
    palladium_mass[UInt(106)] =  105.903486000000001;
    palladium_mass[UInt(108)] =  107.903891999999999;
    palladium_mass[UInt(110)] =  109.905152999999999;

    IsotopeDistribution palladium_isotopes = parseIsotopeDistribution_(palladium_abundance, palladium_mass);
    double palladium_avg_weight = calculateAvgWeight_(palladium_abundance, palladium_mass);
    double palladium_mono_weight = calculateMonoWeight_(palladium_mass);
    String palladium_name = "Palladium";
    String palladium_symbol = "Pd";
    UInt palladium_an = UInt(46);

    Element* palladium = new Element(palladium_name, palladium_symbol, palladium_an, palladium_avg_weight, palladium_mono_weight, palladium_isotopes);
    addElementToMaps(palladium_name, palladium_symbol, palladium_an, palladium);
    storeIsotopes(palladium_name, palladium_symbol, palladium_an, palladium_mass, palladium_isotopes);


    Map<UInt, double> silver_abundance;
    silver_abundance[UInt(107)] =  0.51839;
    silver_abundance[UInt(109)] =  0.48161000000000004;

    Map<UInt, double> silver_mass;
    silver_mass[UInt(107)] =  106.905092999999994;
    silver_mass[UInt(109)] =  108.904756000000006;

    IsotopeDistribution silver_isotopes = parseIsotopeDistribution_(silver_abundance, silver_mass);
    double silver_avg_weight = calculateAvgWeight_(silver_abundance, silver_mass);
    double silver_mono_weight = calculateMonoWeight_(silver_mass);
    String silver_name = "Silver";
    String silver_symbol = "Ag";
    UInt silver_an = UInt(47);

    Element* silver = new Element(silver_name, silver_symbol, silver_an, silver_avg_weight, silver_mono_weight, silver_isotopes);
    addElementToMaps(silver_name, silver_symbol, silver_an, silver);
    storeIsotopes(silver_name, silver_symbol, silver_an, silver_mass, silver_isotopes);


    Map<UInt, double> cadmium_abundance;
    cadmium_abundance[UInt(106)] =  0.0125;
    cadmium_abundance[UInt(108)] =  0.0089;
    cadmium_abundance[UInt(110)] =  0.1249;
    cadmium_abundance[UInt(111)] =  0.128;
    cadmium_abundance[UInt(112)] =  0.2413;
    cadmium_abundance[UInt(113)] =  0.1222;
    cadmium_abundance[UInt(114)] =  0.2873;
    cadmium_abundance[UInt(116)] =  0.07490000000000001;

    Map<UInt, double> cadmium_mass;
    cadmium_mass[UInt(106)] =  105.906458000000001;
    cadmium_mass[UInt(108)] =  107.904184000000001;
    cadmium_mass[UInt(110)] =  109.903002099999995;
    cadmium_mass[UInt(111)] =  110.904178099999996;
    cadmium_mass[UInt(112)] =  111.902757800000003;
    cadmium_mass[UInt(113)] =  112.904401699999994;
    cadmium_mass[UInt(114)] =  113.903358499999996;
    cadmium_mass[UInt(116)] =  115.904756000000006;

    IsotopeDistribution cadmium_isotopes = parseIsotopeDistribution_(cadmium_abundance, cadmium_mass);
    double cadmium_avg_weight = calculateAvgWeight_(cadmium_abundance, cadmium_mass);
    double cadmium_mono_weight = calculateMonoWeight_(cadmium_mass);
    String cadmium_name = "Cadmium";
    String cadmium_symbol = "Cd";
    UInt cadmium_an = UInt(48);

    Element* cadmium = new Element(cadmium_name, cadmium_symbol, cadmium_an, cadmium_avg_weight, cadmium_mono_weight, cadmium_isotopes);
    addElementToMaps(cadmium_name, cadmium_symbol, cadmium_an, cadmium);
    storeIsotopes(cadmium_name, cadmium_symbol, cadmium_an, cadmium_mass, cadmium_isotopes);


    Map<UInt, double> indium_abundance;
    indium_abundance[UInt(113)] =  0.0429;
    indium_abundance[UInt(115)] =  0.9571;

    Map<UInt, double> indium_mass;
    indium_mass[UInt(113)] =  112.904060000000001;
    indium_mass[UInt(115)] =  114.903878000000006;

    IsotopeDistribution indium_isotopes = parseIsotopeDistribution_(indium_abundance, indium_mass);
    double indium_avg_weight = calculateAvgWeight_(indium_abundance, indium_mass);
    double indium_mono_weight = calculateMonoWeight_(indium_mass);
    String indium_name = "Indium";
    String indium_symbol = "In";
    UInt indium_an = UInt(49);

    Element* indium = new Element(indium_name, indium_symbol, indium_an, indium_avg_weight, indium_mono_weight, indium_isotopes);
    addElementToMaps(indium_name, indium_symbol, indium_an, indium);
    storeIsotopes(indium_name, indium_symbol, indium_an, indium_mass, indium_isotopes);


    Map<UInt, double> iodine_abundance;
    iodine_abundance[UInt(127)] =  1.0;

    Map<UInt, double> iodine_mass;
    iodine_mass[UInt(127)] =  126.904472999999996;

    IsotopeDistribution iodine_isotopes = parseIsotopeDistribution_(iodine_abundance, iodine_mass);
    double iodine_avg_weight = calculateAvgWeight_(iodine_abundance, iodine_mass);
    double iodine_mono_weight = calculateMonoWeight_(iodine_mass);
    String iodine_name = "Iodine";
    String iodine_symbol = "I";
    UInt iodine_an = UInt(53);

    Element* iodine = new Element(iodine_name, iodine_symbol, iodine_an, iodine_avg_weight, iodine_mono_weight, iodine_isotopes);
    addElementToMaps(iodine_name, iodine_symbol, iodine_an, iodine);
    storeIsotopes(iodine_name, iodine_symbol, iodine_an, iodine_mass, iodine_isotopes);


    Map<UInt, double> xenon_abundance;
    xenon_abundance[UInt(128)] =  0.0191;
    xenon_abundance[UInt(129)] =  0.264;
    xenon_abundance[UInt(130)] =  0.040999999999999995;
    xenon_abundance[UInt(131)] =  0.212;
    xenon_abundance[UInt(132)] =  0.26899999999999996;
    xenon_abundance[UInt(134)] =  0.10400000000000001;
    xenon_abundance[UInt(136)] =  0.08900000000000001;

    Map<UInt, double> xenon_mass;
    xenon_mass[UInt(128)] =  127.903531000000001;
    xenon_mass[UInt(129)] =  128.904779999999988;
    xenon_mass[UInt(130)] =  129.903509000000014;
    xenon_mass[UInt(131)] =  130.90507199999999;
    xenon_mass[UInt(132)] =  131.904144000000002;
    xenon_mass[UInt(134)] =  133.905394999999999;
    xenon_mass[UInt(136)] =  135.90721400000001;

    IsotopeDistribution xenon_isotopes = parseIsotopeDistribution_(xenon_abundance, xenon_mass);
    double xenon_avg_weight = calculateAvgWeight_(xenon_abundance, xenon_mass);
    double xenon_mono_weight = calculateMonoWeight_(xenon_mass);
    String xenon_name = "Xenon";
    String xenon_symbol = "Xe";
    UInt xenon_an = UInt(54);

    Element* xenon = new Element(xenon_name, xenon_symbol, xenon_an, xenon_avg_weight, xenon_mono_weight, xenon_isotopes);
    addElementToMaps(xenon_name, xenon_symbol, xenon_an, xenon);
    storeIsotopes(xenon_name, xenon_symbol, xenon_an, xenon_mass, xenon_isotopes);


    Map<UInt, double> caesium_abundance;
    caesium_abundance[UInt(133)] =  1.0;

    Map<UInt, double> caesium_mass;
    caesium_mass[UInt(133)] =  132.905451932999995;

    IsotopeDistribution caesium_isotopes = parseIsotopeDistribution_(caesium_abundance, caesium_mass);
    double caesium_avg_weight = calculateAvgWeight_(caesium_abundance, caesium_mass);
    double caesium_mono_weight = calculateMonoWeight_(caesium_mass);
    String caesium_name = "Caesium";
    String caesium_symbol = "Cs";
    UInt caesium_an = UInt(55);

    Element* caesium = new Element(caesium_name, caesium_symbol, caesium_an, caesium_avg_weight, caesium_mono_weight, caesium_isotopes);
    addElementToMaps(caesium_name, caesium_symbol, caesium_an, caesium);
    storeIsotopes(caesium_name, caesium_symbol, caesium_an, caesium_mass, caesium_isotopes);


    Map<UInt, double> cerium_abundance;
    cerium_abundance[UInt(136)] =  0.00185;
    cerium_abundance[UInt(138)] =  0.00251;
    cerium_abundance[UInt(140)] =  0.8845000000000001;
    cerium_abundance[UInt(142)] =  0.11114;

    Map<UInt, double> cerium_mass;
    cerium_mass[UInt(136)] =  135.907172000000003;
    cerium_mass[UInt(138)] =  137.905991;
    cerium_mass[UInt(140)] =  139.905438699999991;
    cerium_mass[UInt(142)] =  141.909244000000001;

    IsotopeDistribution cerium_isotopes = parseIsotopeDistribution_(cerium_abundance, cerium_mass);
    double cerium_avg_weight = calculateAvgWeight_(cerium_abundance, cerium_mass);
    double cerium_mono_weight = calculateMonoWeight_(cerium_mass);
    String cerium_name = "Cerium";
    String cerium_symbol = "Ce";
    UInt cerium_an = UInt(58);

    Element* cerium = new Element(cerium_name, cerium_symbol, cerium_an, cerium_avg_weight, cerium_mono_weight, cerium_isotopes);
    addElementToMaps(cerium_name, cerium_symbol, cerium_an, cerium);
    storeIsotopes(cerium_name, cerium_symbol, cerium_an, cerium_mass, cerium_isotopes);


    Map<UInt, double> praseodymium_abundance;
    praseodymium_abundance[UInt(141)] =  1.0;

    Map<UInt, double> praseodymium_mass;
    praseodymium_mass[UInt(141)] =  140.907646999999997;

    IsotopeDistribution praseodymium_isotopes = parseIsotopeDistribution_(praseodymium_abundance, praseodymium_mass);
    double praseodymium_avg_weight = calculateAvgWeight_(praseodymium_abundance, praseodymium_mass);
    double praseodymium_mono_weight = calculateMonoWeight_(praseodymium_mass);
    String praseodymium_name = "Praseodymium";
    String praseodymium_symbol = "Pr";
    UInt praseodymium_an = UInt(59);

    Element* praseodymium = new Element(praseodymium_name, praseodymium_symbol, praseodymium_an, praseodymium_avg_weight, praseodymium_mono_weight, praseodymium_isotopes);
    addElementToMaps(praseodymium_name, praseodymium_symbol, praseodymium_an, praseodymium);
    storeIsotopes(praseodymium_name, praseodymium_symbol, praseodymium_an, praseodymium_mass, praseodymium_isotopes);


    Map<UInt, double> gadolinium_abundance;
    gadolinium_abundance[UInt(152)] =  0.002;
    gadolinium_abundance[UInt(154)] =  0.0218;
    gadolinium_abundance[UInt(155)] =  0.14800000000000002;
    gadolinium_abundance[UInt(156)] =  0.2047;
    gadolinium_abundance[UInt(157)] =  0.1565;
    gadolinium_abundance[UInt(158)] =  0.2484;
    gadolinium_abundance[UInt(160)] =  0.2186;

    Map<UInt, double> gadolinium_mass;
    gadolinium_mass[UInt(152)] =  151.919791000000004;
    gadolinium_mass[UInt(154)] =  153.920865600000013;
    gadolinium_mass[UInt(155)] =  154.92262199999999;
    gadolinium_mass[UInt(156)] =  155.922122699999989;
    gadolinium_mass[UInt(157)] =  156.923960099999988;
    gadolinium_mass[UInt(158)] =  157.924103900000006;
    gadolinium_mass[UInt(160)] =  159.927054099999992;

    IsotopeDistribution gadolinium_isotopes = parseIsotopeDistribution_(gadolinium_abundance, gadolinium_mass);
    double gadolinium_avg_weight = calculateAvgWeight_(gadolinium_abundance, gadolinium_mass);
    double gadolinium_mono_weight = calculateMonoWeight_(gadolinium_mass);
    String gadolinium_name = "Gadolinium";
    String gadolinium_symbol = "Gd";
    UInt gadolinium_an = UInt(64);

    Element* gadolinium = new Element(gadolinium_name, gadolinium_symbol, gadolinium_an, gadolinium_avg_weight, gadolinium_mono_weight, gadolinium_isotopes);
    addElementToMaps(gadolinium_name, gadolinium_symbol, gadolinium_an, gadolinium);
    storeIsotopes(gadolinium_name, gadolinium_symbol, gadolinium_an, gadolinium_mass, gadolinium_isotopes);


    Map<UInt, double> hafnium_abundance;
    hafnium_abundance[UInt(176)] =  0.0526;
    hafnium_abundance[UInt(177)] =  0.18600000000000003;
    hafnium_abundance[UInt(178)] =  0.2728;
    hafnium_abundance[UInt(179)] =  0.1362;
    hafnium_abundance[UInt(180)] =  0.3508;

    Map<UInt, double> hafnium_mass;
    hafnium_mass[UInt(176)] =  175.941408599999988;
    hafnium_mass[UInt(177)] =  176.943220700000012;
    hafnium_mass[UInt(178)] =  177.943698799999993;
    hafnium_mass[UInt(179)] =  178.945816100000002;
    hafnium_mass[UInt(180)] =  179.946550000000002;

    IsotopeDistribution hafnium_isotopes = parseIsotopeDistribution_(hafnium_abundance, hafnium_mass);
    double hafnium_avg_weight = calculateAvgWeight_(hafnium_abundance, hafnium_mass);
    double hafnium_mono_weight = calculateMonoWeight_(hafnium_mass);
    String hafnium_name = "Hafnium";
    String hafnium_symbol = "Hf";
    UInt hafnium_an = UInt(72);

    Element* hafnium = new Element(hafnium_name, hafnium_symbol, hafnium_an, hafnium_avg_weight, hafnium_mono_weight, hafnium_isotopes);
    addElementToMaps(hafnium_name, hafnium_symbol, hafnium_an, hafnium);
    storeIsotopes(hafnium_name, hafnium_symbol, hafnium_an, hafnium_mass, hafnium_isotopes);


    Map<UInt, double> tantalum_abundance;
    tantalum_abundance[UInt(181)] =  1.0;

    Map<UInt, double> tantalum_mass;
    tantalum_mass[UInt(181)] =  180.947995800000001;

    IsotopeDistribution tantalum_isotopes = parseIsotopeDistribution_(tantalum_abundance, tantalum_mass);
    double tantalum_avg_weight = calculateAvgWeight_(tantalum_abundance, tantalum_mass);
    double tantalum_mono_weight = calculateMonoWeight_(tantalum_mass);
    String tantalum_name = "Tantalum";
    String tantalum_symbol = "Ta";
    UInt tantalum_an = UInt(73);

    Element* tantalum = new Element(tantalum_name, tantalum_symbol, tantalum_an, tantalum_avg_weight, tantalum_mono_weight, tantalum_isotopes);
    addElementToMaps(tantalum_name, tantalum_symbol, tantalum_an, tantalum);
    storeIsotopes(tantalum_name, tantalum_symbol, tantalum_an, tantalum_mass, tantalum_isotopes);


    Map<UInt, double> platinum_abundance;
    platinum_abundance[UInt(192)] =  0.00782;
    platinum_abundance[UInt(194)] =  0.32966999999999996;
    platinum_abundance[UInt(195)] =  0.33832;
    platinum_abundance[UInt(196)] =  0.25242000000000003;
    platinum_abundance[UInt(198)] =  0.07163;

    Map<UInt, double> platinum_mass;
    platinum_mass[UInt(192)] =  191.961038000000002;
    platinum_mass[UInt(194)] =  193.962680299999988;
    platinum_mass[UInt(195)] =  194.964791100000014;
    platinum_mass[UInt(196)] =  195.964951500000012;
    platinum_mass[UInt(198)] =  197.967893000000004;

    IsotopeDistribution platinum_isotopes = parseIsotopeDistribution_(platinum_abundance, platinum_mass);
    double platinum_avg_weight = calculateAvgWeight_(platinum_abundance, platinum_mass);
    double platinum_mono_weight = calculateMonoWeight_(platinum_mass);
    String platinum_name = "Platinum";
    String platinum_symbol = "Pt";
    UInt platinum_an = UInt(78);

    Element* platinum = new Element(platinum_name, platinum_symbol, platinum_an, platinum_avg_weight, platinum_mono_weight, platinum_isotopes);
    addElementToMaps(platinum_name, platinum_symbol, platinum_an, platinum);
    storeIsotopes(platinum_name, platinum_symbol, platinum_an, platinum_mass, platinum_isotopes);


    Map<UInt, double> tungsten_abundance;
    tungsten_abundance[UInt(180)] =  0.0012;
    tungsten_abundance[UInt(182)] =  0.265;
    tungsten_abundance[UInt(183)] =  0.1431;
    tungsten_abundance[UInt(184)] =  0.3064;
    tungsten_abundance[UInt(186)] =  0.2843;

    Map<UInt, double> tungsten_mass;
    tungsten_mass[UInt(180)] =  179.946704000000011;
    tungsten_mass[UInt(182)] =  181.948204199999992;
    tungsten_mass[UInt(183)] =  182.950222999999994;
    tungsten_mass[UInt(184)] =  183.950930999999997;
    tungsten_mass[UInt(186)] =  185.954364099999992;

    IsotopeDistribution tungsten_isotopes = parseIsotopeDistribution_(tungsten_abundance, tungsten_mass);
    double tungsten_avg_weight = calculateAvgWeight_(tungsten_abundance, tungsten_mass);
    double tungsten_mono_weight = calculateMonoWeight_(tungsten_mass);
    String tungsten_name = "Tungsten";
    String tungsten_symbol = "W";
    UInt tungsten_an = UInt(74);

    Element* tungsten = new Element(tungsten_name, tungsten_symbol, tungsten_an, tungsten_avg_weight, tungsten_mono_weight, tungsten_isotopes);
    addElementToMaps(tungsten_name, tungsten_symbol, tungsten_an, tungsten);
    storeIsotopes(tungsten_name, tungsten_symbol, tungsten_an, tungsten_mass, tungsten_isotopes);


    Map<UInt, double> gold_abundance;
    gold_abundance[UInt(197)] =  1.0;

    Map<UInt, double> gold_mass;
    gold_mass[UInt(197)] =  196.96655100000001;

    IsotopeDistribution gold_isotopes = parseIsotopeDistribution_(gold_abundance, gold_mass);
    double gold_avg_weight = calculateAvgWeight_(gold_abundance, gold_mass);
    double gold_mono_weight = calculateMonoWeight_(gold_mass);
    String gold_name = "Gold";
    String gold_symbol = "Au";
    UInt gold_an = UInt(79);

    Element* gold = new Element(gold_name, gold_symbol, gold_an, gold_avg_weight, gold_mono_weight, gold_isotopes);
    addElementToMaps(gold_name, gold_symbol, gold_an, gold);
    storeIsotopes(gold_name, gold_symbol, gold_an, gold_mass, gold_isotopes);


    Map<UInt, double> mercury_abundance;
    mercury_abundance[UInt(196)] =  0.0015;
    mercury_abundance[UInt(198)] =  0.09970000000000001;
    mercury_abundance[UInt(199)] =  0.16870000000000002;
    mercury_abundance[UInt(200)] =  0.231;
    mercury_abundance[UInt(201)] =  0.1318;
    mercury_abundance[UInt(202)] =  0.2986;
    mercury_abundance[UInt(204)] =  0.0687;

    Map<UInt, double> mercury_mass;
    mercury_mass[UInt(196)] =  195.965833000000004;
    mercury_mass[UInt(198)] =  197.966768999999999;
    mercury_mass[UInt(199)] =  198.968279899999999;
    mercury_mass[UInt(200)] =  199.968325999999991;
    mercury_mass[UInt(201)] =  200.970302299999986;
    mercury_mass[UInt(202)] =  201.970642999999996;
    mercury_mass[UInt(204)] =  203.973493899999994;

    IsotopeDistribution mercury_isotopes = parseIsotopeDistribution_(mercury_abundance, mercury_mass);
    double mercury_avg_weight = calculateAvgWeight_(mercury_abundance, mercury_mass);
    double mercury_mono_weight = calculateMonoWeight_(mercury_mass);
    String mercury_name = "Mercury";
    String mercury_symbol = "Hg";
    UInt mercury_an = UInt(80);

    Element* mercury = new Element(mercury_name, mercury_symbol, mercury_an, mercury_avg_weight, mercury_mono_weight, mercury_isotopes);
    addElementToMaps(mercury_name, mercury_symbol, mercury_an, mercury);
    storeIsotopes(mercury_name, mercury_symbol, mercury_an, mercury_mass, mercury_isotopes);


    Map<UInt, double> thallium_abundance;
    thallium_abundance[UInt(203)] =  0.2952;
    thallium_abundance[UInt(205)] =  0.7048000000000001;

    Map<UInt, double> thallium_mass;
    thallium_mass[UInt(203)] =  202.972344200000009;
    thallium_mass[UInt(205)] =  204.97442749999999;

    IsotopeDistribution thallium_isotopes = parseIsotopeDistribution_(thallium_abundance, thallium_mass);
    double thallium_avg_weight = calculateAvgWeight_(thallium_abundance, thallium_mass);
    double thallium_mono_weight = calculateMonoWeight_(thallium_mass);
    String thallium_name = "Thallium";
    String thallium_symbol = "Tl";
    UInt thallium_an = UInt(81);

    Element* thallium = new Element(thallium_name, thallium_symbol, thallium_an, thallium_avg_weight, thallium_mono_weight, thallium_isotopes);
    addElementToMaps(thallium_name, thallium_symbol, thallium_an, thallium);
    storeIsotopes(thallium_name, thallium_symbol, thallium_an, thallium_mass, thallium_isotopes);


    Map<UInt, double> lead_abundance;
    lead_abundance[UInt(204)] =  0.013999999999999999;
    lead_abundance[UInt(206)] =  0.24100000000000002;
    lead_abundance[UInt(207)] =  0.221;
    lead_abundance[UInt(208)] =  0.524;

    Map<UInt, double> lead_mass;
    lead_mass[UInt(204)] =  203.973043600000011;
    lead_mass[UInt(206)] =  205.974465299999991;
    lead_mass[UInt(207)] =  206.975896900000009;
    lead_mass[UInt(208)] =  207.976653800000008;

    IsotopeDistribution lead_isotopes = parseIsotopeDistribution_(lead_abundance, lead_mass);
    double lead_avg_weight = calculateAvgWeight_(lead_abundance, lead_mass);
    double lead_mono_weight = calculateMonoWeight_(lead_mass);
    String lead_name = "Lead";
    String lead_symbol = "Pb";
    UInt lead_an = UInt(82);

    Element* lead = new Element(lead_name, lead_symbol, lead_an, lead_avg_weight, lead_mono_weight, lead_isotopes);
    addElementToMaps(lead_name, lead_symbol, lead_an, lead);
    storeIsotopes(lead_name, lead_symbol, lead_an, lead_mass, lead_isotopes);


    Map<UInt, double> bismuth_abundance;
    bismuth_abundance[UInt(209)] =  1.0;

    Map<UInt, double> bismuth_mass;
    bismuth_mass[UInt(209)] =  208.980398699999995;

    IsotopeDistribution bismuth_isotopes = parseIsotopeDistribution_(bismuth_abundance, bismuth_mass);
    double bismuth_avg_weight = calculateAvgWeight_(bismuth_abundance, bismuth_mass);
    double bismuth_mono_weight = calculateMonoWeight_(bismuth_mass);
    String bismuth_name = "Bismuth";
    String bismuth_symbol = "Bi";
    UInt bismuth_an = UInt(83);

    Element* bismuth = new Element(bismuth_name, bismuth_symbol, bismuth_an, bismuth_avg_weight, bismuth_mono_weight, bismuth_isotopes);
    addElementToMaps(bismuth_name, bismuth_symbol, bismuth_an, bismuth);
    storeIsotopes(bismuth_name, bismuth_symbol, bismuth_an, bismuth_mass, bismuth_isotopes);


    Map<UInt, double> rhenium_abundance;
    rhenium_abundance[UInt(185)] =  0.374;
    rhenium_abundance[UInt(187)] =  0.626;

    Map<UInt, double> rhenium_mass;
    rhenium_mass[UInt(185)] =  184.952955000000003;
    rhenium_mass[UInt(187)] =  186.95575310000001;

    IsotopeDistribution rhenium_isotopes = parseIsotopeDistribution_(rhenium_abundance, rhenium_mass);
    double rhenium_avg_weight = calculateAvgWeight_(rhenium_abundance, rhenium_mass);
    double rhenium_mono_weight = calculateMonoWeight_(rhenium_mass);
    String rhenium_name = "Rhenium";
    String rhenium_symbol = "Re";
    UInt rhenium_an = UInt(75);

    Element* rhenium = new Element(rhenium_name, rhenium_symbol, rhenium_an, rhenium_avg_weight, rhenium_mono_weight, rhenium_isotopes);
    addElementToMaps(rhenium_name, rhenium_symbol, rhenium_an, rhenium);
    storeIsotopes(rhenium_name, rhenium_symbol, rhenium_an, rhenium_mass, rhenium_isotopes);


    Map<UInt, double> neodymium_abundance;
    neodymium_abundance[UInt(142)] =  0.272;
    neodymium_abundance[UInt(143)] =  0.122;
    neodymium_abundance[UInt(144)] =  0.23800000000000002;
    neodymium_abundance[UInt(145)] =  0.083;
    neodymium_abundance[UInt(146)] =  0.172;
    neodymium_abundance[UInt(148)] =  0.057999999999999996;
    neodymium_abundance[UInt(150)] =  0.055999999999999994;

    Map<UInt, double> neodymium_mass;
    neodymium_mass[UInt(142)] =  141.907723299999987;
    neodymium_mass[UInt(143)] =  142.909814299999994;
    neodymium_mass[UInt(144)] =  143.910087299999987;
    neodymium_mass[UInt(145)] =  144.912573600000002;
    neodymium_mass[UInt(146)] =  145.913116900000006;
    neodymium_mass[UInt(148)] =  147.916892999999988;
    neodymium_mass[UInt(150)] =  149.920891000000012;

    IsotopeDistribution neodymium_isotopes = parseIsotopeDistribution_(neodymium_abundance, neodymium_mass);
    double neodymium_avg_weight = calculateAvgWeight_(neodymium_abundance, neodymium_mass);
    double neodymium_mono_weight = calculateMonoWeight_(neodymium_mass);
    String neodymium_name = "Neodymium";
    String neodymium_symbol = "Nd";
    UInt neodymium_an = UInt(60);

    Element* neodymium = new Element(neodymium_name, neodymium_symbol, neodymium_an, neodymium_avg_weight, neodymium_mono_weight, neodymium_isotopes);
    addElementToMaps(neodymium_name, neodymium_symbol, neodymium_an, neodymium);
    storeIsotopes(neodymium_name, neodymium_symbol, neodymium_an, neodymium_mass, neodymium_isotopes);


    Map<UInt, double> thorium_abundance;
    thorium_abundance[UInt(230)] =  0.0002;
    thorium_abundance[UInt(232)] =  0.9998;

    Map<UInt, double> thorium_mass;
    thorium_mass[UInt(230)] =  230.033133800000002;
    thorium_mass[UInt(232)] =  232.038055299999996;

    IsotopeDistribution thorium_isotopes = parseIsotopeDistribution_(thorium_abundance, thorium_mass);
    double thorium_avg_weight = calculateAvgWeight_(thorium_abundance, thorium_mass);
    double thorium_mono_weight = calculateMonoWeight_(thorium_mass);
    String thorium_name = "Thorium";
    String thorium_symbol = "Th";
    UInt thorium_an = UInt(90);

    Element* thorium = new Element(thorium_name, thorium_symbol, thorium_an, thorium_avg_weight, thorium_mono_weight, thorium_isotopes);
    addElementToMaps(thorium_name, thorium_symbol, thorium_an, thorium);
    storeIsotopes(thorium_name, thorium_symbol, thorium_an, thorium_mass, thorium_isotopes);


    Map<UInt, double> lanthanum_abundance;
    lanthanum_abundance[UInt(138)] =  0.00089;
    lanthanum_abundance[UInt(139)] =  0.99911;

    Map<UInt, double> lanthanum_mass;
    lanthanum_mass[UInt(138)] =  137.907112000000012;
    lanthanum_mass[UInt(139)] =  138.906353300000006;

    IsotopeDistribution lanthanum_isotopes = parseIsotopeDistribution_(lanthanum_abundance, lanthanum_mass);
    double lanthanum_avg_weight = calculateAvgWeight_(lanthanum_abundance, lanthanum_mass);
    double lanthanum_mono_weight = calculateMonoWeight_(lanthanum_mass);
    String lanthanum_name = "Lanthanum";
    String lanthanum_symbol = "La";
    UInt lanthanum_an = UInt(57);

    Element* lanthanum = new Element(lanthanum_name, lanthanum_symbol, lanthanum_an, lanthanum_avg_weight, lanthanum_mono_weight, lanthanum_isotopes);
    addElementToMaps(lanthanum_name, lanthanum_symbol, lanthanum_an, lanthanum);
    storeIsotopes(lanthanum_name, lanthanum_symbol, lanthanum_an, lanthanum_mass, lanthanum_isotopes);


    Map<UInt, double> samarium_abundance;
    samarium_abundance[UInt(144)] =  0.0308;
    samarium_abundance[UInt(147)] =  0.15;
    samarium_abundance[UInt(148)] =  0.1125;
    samarium_abundance[UInt(149)] =  0.1382;
    samarium_abundance[UInt(150)] =  0.0737;
    samarium_abundance[UInt(152)] =  0.26739999999999997;
    samarium_abundance[UInt(154)] =  0.2274;

    Map<UInt, double> samarium_mass;
    samarium_mass[UInt(144)] =  143.911999000000009;
    samarium_mass[UInt(147)] =  146.9148979;
    samarium_mass[UInt(148)] =  147.914822700000002;
    samarium_mass[UInt(149)] =  148.917184700000007;
    samarium_mass[UInt(150)] =  149.917275499999988;
    samarium_mass[UInt(152)] =  151.919732399999987;
    samarium_mass[UInt(154)] =  153.92220929999999;

    IsotopeDistribution samarium_isotopes = parseIsotopeDistribution_(samarium_abundance, samarium_mass);
    double samarium_avg_weight = calculateAvgWeight_(samarium_abundance, samarium_mass);
    double samarium_mono_weight = calculateMonoWeight_(samarium_mass);
    String samarium_name = "Samarium";
    String samarium_symbol = "Sm";
    UInt samarium_an = UInt(62);

    Element* samarium = new Element(samarium_name, samarium_symbol, samarium_an, samarium_avg_weight, samarium_mono_weight, samarium_isotopes);
    addElementToMaps(samarium_name, samarium_symbol, samarium_an, samarium);
    storeIsotopes(samarium_name, samarium_symbol, samarium_an, samarium_mass, samarium_isotopes);

}

  void ElementDB::addElementToMaps(const String& name, const String& symbol, const UInt an, const Element* e)
  {
    names_[name] = e;
    names_[symbol] = e;
    names_[an] = e;
  }

  void ElementDB::storeIsotopes(const String& name, const String& symbol, const UInt an, const Map<UInt, double>& Z_to_mass, const IsotopeDistribution& isotopes)
  {
    for (const auto& isotope : isotopes)
    {
      double atomic_mass = isotope.getMZ();
      UInt mass_number = round(atomic_mass);
      String iso_name = "(" + String(mass_number) + ")" + name;
      String iso_symbol = "(" + String(mass_number) + ")" + symbol;

      // set avg and mono to same value for isotopes (old hack...)
      double iso_avg_weight = Z_to_mass[mass_number];
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

  IsotopeDistribution ElementDB::parseIsotopeDistribution_(const Map<UInt, double>& Z_to_abundance, const Map<UInt, double>& Z_to_mass)
  {
    IsotopeDistribution::ContainerType dist;
    
    vector<UInt> keys;
    for (Map<UInt, double>::const_iterator it = Z_to_abundance.begin(); it != Z_to_abundance.end(); ++it)
    {
      keys.push_back(it->first);
    }

    // calculate weighted average
    for (vector<UInt>::iterator it = keys.begin(); it != keys.end(); ++it)
    {
      dist.push_back(Peak1D(Z_to_mass[*it] , Z_to_abundance[*it]));
    }


    IsotopeDistribution iso_dist;
    iso_dist.set(dist);
    
    return iso_dist;
  }

  void ElementDB::clear_()
  {
    // names_ has the union of all Element*, deleting this is sufficient to avoid mem leaks
    Map<String, const Element*>::Iterator it = names_.begin();
    for (; it != names_.end(); ++it)
    {
      delete it->second;
    }
    names_.clear();
    symbols_.clear();
    atomic_numbers_.clear();
  }

} // namespace OpenMS
