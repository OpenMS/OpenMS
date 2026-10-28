// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SearchParameters.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <set>
#include <tuple>

using namespace std;

namespace OpenMS
{
  SearchParameters::SearchParameters() :
      db(),
      db_version(),
      taxonomy(),
      charges(),
      mass_type(PeakMassType::MONOISOTOPIC),
      fixed_modifications(),
      variable_modifications(),
      missed_cleavages(0),
      fragment_mass_tolerance(0.0),
      fragment_mass_tolerance_ppm(false),
      precursor_mass_tolerance(0.0),
      precursor_mass_tolerance_ppm(false),
      digestion_enzyme("unknown_enzyme", ""),
      enzyme_term_specificity(EnzymaticDigestion::SPEC_UNKNOWN)
  {
  }

  bool SearchParameters::operator==(const SearchParameters& rhs) const
  {
    return
        std::tie(db, db_version, taxonomy, charges, mass_type, fixed_modifications, variable_modifications,
            missed_cleavages, fragment_mass_tolerance, fragment_mass_tolerance_ppm, precursor_mass_tolerance,
            precursor_mass_tolerance_ppm, digestion_enzyme, enzyme_term_specificity) ==
        std::tie(rhs.db, rhs.db_version, rhs.taxonomy, rhs.charges, rhs.mass_type, rhs.fixed_modifications,
            rhs.variable_modifications, rhs.missed_cleavages, rhs.fragment_mass_tolerance,
            rhs.fragment_mass_tolerance_ppm, rhs.precursor_mass_tolerance,
            rhs.precursor_mass_tolerance_ppm, rhs.digestion_enzyme, rhs.enzyme_term_specificity);
  }

  bool SearchParameters::operator!=(const SearchParameters& rhs) const
  {
    return !(*this == rhs);
  }

  bool SearchParameters::mergeable(const SearchParameters& sp, const std::string& experiment_type) const
  {
    std::string spdb = sp.db;
    StringUtils::substitute(spdb, "\\","/");
    std::string pdb = this->db;
    StringUtils::substitute(pdb, "\\","/");

    if  (this->precursor_mass_tolerance != sp.precursor_mass_tolerance ||
        this->precursor_mass_tolerance_ppm != sp.precursor_mass_tolerance_ppm ||
        File::basename(pdb) != File::basename(spdb) ||
        this->db_version != sp.db_version ||
        this->fragment_mass_tolerance != sp.fragment_mass_tolerance ||
        this->fragment_mass_tolerance_ppm != sp.fragment_mass_tolerance_ppm ||
        this->charges != sp.charges ||
        this->digestion_enzyme != sp.digestion_enzyme ||
        this->taxonomy != sp.taxonomy ||
         this->enzyme_term_specificity != sp.enzyme_term_specificity)
    {
      return false;
    }

    set<std::string> fixed_mods(this->fixed_modifications.begin(), this->fixed_modifications.end());
    set<std::string> var_mods(this->variable_modifications.begin(), this->variable_modifications.end());
    set<std::string> curr_fixed_mods(sp.fixed_modifications.begin(), sp.fixed_modifications.end());
    set<std::string> curr_var_mods(sp.variable_modifications.begin(), sp.variable_modifications.end());
    if (fixed_mods != curr_fixed_mods ||
        var_mods != curr_var_mods)
    {
      if (experiment_type != "labeled_MS1")
      {
        return false;
      }
      else
      {
        //TODO actually introduce a flag for labelling modifications in the Mod datastructures?
        //OR put a unique ID for the used mod as a UserParam to the mapList entries (consensusHeaders)
        //TODO actually you would probably need an experimental design here, because
        //settings have to agree exactly in a FractionGroup but can slightly differ across runs.
        //Or just ignore labelling mods during the check
        return true;
      }
    }
    return true;
  }

  int SearchParameters::getChargeValue_(std::string& charge_str) const
  {
    // We have to do this because some people/tools put the + or - AFTER the number...
    bool neg = StringUtils::hasSubstring(charge_str, '-');
    neg ? StringUtils::remove(charge_str, '-') : StringUtils::remove(charge_str, '+');
    int val = StringUtils::toInt32(charge_str);
    return neg ? -val : val;
  }

  std::pair<int,int> SearchParameters::getChargeRange() const
  {
    std::pair<int,int> result{0,0};

    try // is there only one number (min = max)?
    {
      result.first = StringUtils::toInt32(charges);
      result.second = result.first;
    }
    catch (Exception::ConversionError&) // nope, something else
    {
      if (StringUtils::hasSubstring(charges, ',')) // it's probably a list
      {
        IntList chgs = ListUtils::create<Int>(charges);
        auto minmax = minmax_element(chgs.begin(), chgs.end());
        result.first = *minmax.first;
        result.second = *minmax.second;
      }
      else if (StringUtils::hasSubstring(charges, ':')) // it's probably a range
      {
        StringList chgs;
        StringUtils::split(charges, ':', chgs);
        if (chgs.size() > 2)
        {
          throw OpenMS::Exception::MissingInformation(
            __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Charge string in SearchParameters not parseable.");
        }
        result.first = getChargeValue_(chgs[0]);
        result.second = getChargeValue_(chgs[1]);
      }
      else
      {
        size_t pos = charges.find('-', 0);
        std::vector<size_t> minus_positions;
        while (pos != string::npos)
        {
          minus_positions.push_back(pos);
          pos = charges.find('-', pos + 1);
        }
        if (!minus_positions.empty() && minus_positions.size() <= 3) // it's probably a range with '-'
        {
          Size split_pos(0);
          if (minus_positions.size() <= 1) // split at first minus
          {
            split_pos = minus_positions[0];
          }
          else
          {
            split_pos = minus_positions[1];
          }
          std::string first = StringUtils::substr(charges, 0, split_pos);
          std::string second = StringUtils::substr(charges, split_pos + 1, string::npos);
          result.first = getChargeValue_(first);
          result.second = getChargeValue_(second);
        }
      }
    }
    return result;
  }
} // namespace OpenMS
