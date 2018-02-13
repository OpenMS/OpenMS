// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Andreas Bertsch, Timo Sachsenberg $
// --------------------------------------------------------------------------
//
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/SYSTEM/File.h>

#include <cmath>
#include <iostream>

using namespace std;

namespace OpenMS
{
  ElementDB::ElementDB()
  {
    readFromFile_("CHEMISTRY/Elements.xml");
  }

  ElementDB::~ElementDB()
  {
    clear_();
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
    // extract Zs
    vector<UInt> keys;
    for (Map<UInt, double>::const_iterator it = Z_to_abundance.begin(); it != Z_to_abundance.end(); ++it)
    {
      keys.push_back(it->first);
    }

    // calculate weighted average
    for (vector<UInt>::iterator it = keys.begin(); it != keys.end(); ++it)
    {
      avg += Z_to_mass[*it] * Z_to_abundance[*it];
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

  void ElementDB::readFromFile_(const String& file_name)
  {
    String file = File::find(file_name);

    // load elements into param object
    Param param;
    ParamXMLFile paramFile;
    paramFile.load(file, param);

    UInt an(0);
    String name, symbol;

    // determine prefix
    vector<String> split;
    param.begin().getName().split(':', split);
    String prefix("");
    for (Size i = 0; i < split.size() - 1; ++i)
    {
      prefix += split[i] + ":";
    }
    //cout << "first element prefix=" << prefix << endl;

    Map<UInt, double> Z_to_abundancy;
    Map<UInt, double> Z_to_mass;

    for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
    {
      // new element started?
      if (!it.getName().hasPrefix(prefix))
      {
        // update prefix
        it.getName().split(':', split);
        prefix = "";
        for (Size i = 0; i < split.size() - 1; ++i)
        {
          prefix += split[i] + ":";
        }
        // cout << "new element prefix=" << prefix << endl;

        // Parsing of previous element is finished. Now store data in Element object
        IsotopeDistribution isotopes = parseIsotopeDistribution_(Z_to_abundancy, Z_to_mass);
        double avg_weight = calculateAvgWeight_(Z_to_abundancy, Z_to_mass);
        double mono_weight = calculateMonoWeight_(Z_to_mass);

        /*
        // print information about elements
        cout << "Name: " << name << " AtomicNumber: " << an << " Symbol: " << symbol << " AvgWeight: " << avg_weight
             << " MonoWeight: " << mono_weight << " NIsotopes: " << isotopes.size() << endl;

        */
        Element* e = new Element(name, symbol, an, avg_weight, mono_weight, isotopes);
        names_[name] = e;
        symbols_[symbol] = e;
        atomic_numbers_[an] = e;

        // add all the individual isotopes as separate elements
        for (auto const & isotope : isotopes)
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

          /*
          // print name, symbol and atomic mass of the current isotope
          cout << "Isotope Name: " << iso_name << " Symbol: " << iso_symbol << " AtomicMass: " << iso_mono_weight << endl;
          */

          Element* iso_e = new Element(iso_name, iso_symbol, an, iso_avg_weight, iso_mono_weight, iso_isotopes);
          names_[iso_name] = iso_e;
          names_[iso_symbol] = iso_e;
        }

        Z_to_abundancy.clear();
        Z_to_mass.clear();
      }

      // top level: read the contents of the element section
      it.getName().split(':', split);
      String key = split[2];
      String value = it->value;
      value.trim();

      // cout << "Key=" << key << endl;

      if (key == "AtomicNumber")
      {
        an = (UInt)value.toInt();
      }
      else
      {
        if (key == "Isotopes")
        {
          UInt Z = UInt(split[3].toInt());
          String item = split[4];
          if (item == "RelativeAbundance")
          {
            Z_to_abundancy[Z] = double(value.toDouble() / 100.0);
          }
          else if (item == "AtomicMass")
          {
            Z_to_mass[Z] = double(value.toDouble());
          }
          else
          {
            cerr << "read unknown item in Isotopes: " << item << endl;
          }
        }
        else
        {
          if (key == "Name")
          {
            name = value;
          }
          else
          {
            if (key == "Symbol")
            {
              symbol = value;
            }
            else
            {
              cerr << "read unknown tag: " << key << endl;
            }
          }
        }
      }
    }

    // build last element
    double avg_weight(0), mono_weight(0);
    IsotopeDistribution isotopes = parseIsotopeDistribution_(Z_to_abundancy,Z_to_mass);
    Element* e = new Element(name, symbol, an, avg_weight, mono_weight, isotopes);
    names_[name] = e;
    symbols_[symbol] = e;
    atomic_numbers_[an] = e;
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
