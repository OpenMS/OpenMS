// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AdductInfo.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

namespace OpenMS
{
  AdductInfo::AdductInfo(const String& name, const EmpiricalFormula& adduct, int charge, UInt mol_multiplier)
    :
    name_(name),
    ef_(adduct),
    charge_(charge),
    mol_multiplier_(mol_multiplier)
  {
    if (charge_ == 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Charge of 0 is not allowed for an adduct (" + ef_.toString() + ")");
    }
    if (adduct.getCharge() != 0)
    { // EmpiricalFormula adds/subtracts protons to make up the charge;
      // we just use the uncharged formula and take care of charge ourselves
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "EmpiricalFormula must not have a charge (" + ef_.toString() + "), since the internal weight computation of EF is unsuitable for adducts.");
    }
    if (mol_multiplier_ == 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Mol. multiplier of 0 is not allowed for an adduct (" + ef_.toString() + ")");
    }
    mass_ = ef_.getMonoWeight();
  }

  double AdductInfo::getNeutralMass(double observed_mz) const
  {
    // decharge and remove adduct (charge is guaranteed != 0; see C'tor)
    double mass = observed_mz * abs(charge_) - mass_;

    // correct for electron masses
    // (positive charge means there are electrons missing!)
    // (negative charge requires increasing the mass by X electrons)
    // --> looking at observed m/z, we thus need to decharge to get equal protons and electrons
    mass += charge_ * 1 * Constants::ELECTRON_MASS_U;

    // the Mol multiplier determines if we assume to be looking at dimers or higher
    // Currently, we just want the monomer, to compare its mass to a DB entry
    mass /= mol_multiplier_;

    return mass;
  }

  double AdductInfo::getMZ(double neutral_mass) const
  {
    // this is the inverse of getNeutralMass()
    double neutral_nmer_mass_with_adduct = (neutral_mass * mol_multiplier_ + mass_);  // [nM+adduct]

    // correct for electron masses
    // (positive charge means there are electrons missing!)
    // (negative charge requires increasing the mass by X electrons)
    neutral_nmer_mass_with_adduct -= charge_ * Constants::ELECTRON_MASS_U;

    return neutral_nmer_mass_with_adduct / abs(charge_);
  }

  double AdductInfo::getMassShift(bool use_avg_mass) const
  {
    double mass = use_avg_mass ? ef_.getAverageWeight() : mass_;
    // intrinsic adduct charge comes from additional/missing electrons, but for
    // mass shift must be compensated by adding/removing hydrogens:
    return mass - charge_ * (Constants::PROTON_MASS_U + Constants::ELECTRON_MASS_U);
  }

  /// checks if an adduct (e.g.a 'M+2K-H;1+') is valid, i.e. if the losses (==negative amounts) can actually be lost by the compound given in @p db_entry.
  /// If the negative parts are present in @p db_entry, true is returned.
  bool AdductInfo::isCompatible(const EmpiricalFormula& db_entry) const
  {
    return db_entry.contains(ef_ * -1);
  }

  int AdductInfo::getCharge() const
  {
    return charge_;
  }

  const String& AdductInfo::getName() const
  {
    return name_;
  }

  const EmpiricalFormula& AdductInfo::getEmpiricalFormula() const
  {
    return ef_;
  }

  UInt AdductInfo::getMolMultiplier() const
  {
    return mol_multiplier_;
  }

  bool AdductInfo::operator==(const AdductInfo& other) const
  {
    return (name_ == other.name_) && (ef_ == other.ef_) &&
      (charge_ == other.charge_) && (mol_multiplier_ == other.mol_multiplier_);
  }

  AdductInfo AdductInfo::parseAdductString(const String& adduct)
  {
    // adduct string looks like this:
    // M+2K-H;1+   or
    // 2M+CH3CN+Na;1+  (i.e. multimers are supported)

    // do some sanity checks on the string

    // retrieve adduct and charge
    String cp_str(adduct);
    cp_str.removeWhitespaces();
    StringList list;
    cp_str.split(";", list);
    // split term into formula and charge, e.g. "M-H" and "1-"
    String mol_formula, charge_str;
    if (list.size() == 2)
    {
      mol_formula = list[0];
      charge_str = list[1];
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not detect molecular ion; charge in '" + cp_str + "'. Got semicolon right?", cp_str);
    }

    // check if charge string is formatted correctly
    if ((!charge_str.hasSuffix("+")) && (!charge_str.hasSuffix("-")))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Charge sign +/- in the end of the string is missing! ", charge_str);
    }

    // get charge and sign (throws ConversionError if not an integer)
    int charge = charge_str.substr(0, charge_str.size() - 1).toInt();

    if (charge_str.suffix(1) == "+")
    {
      if (charge < 0)
      {
        charge *= -1;
      }
    }
    else
    {
      if (charge > 0)
      {
        charge *= -1;
      }
    }

    // not allowing double ++ or -- or +- or -+
    String op_str(mol_formula);
    op_str.substitute('-', '+');
    if (op_str.hasSubstring("++") || op_str.hasSuffix("+") || op_str.hasPrefix("+"))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "+/- operator must be surrounded by a chemical formula. Offending string: ", mol_formula);
    }

    // split by + and -
    op_str = mol_formula;
    if (op_str.has('%'))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Character '%' not allowed within chemical formula. Offending string: ", mol_formula);
    }
    // ... we want to keep the - and +, so we add extra chars around, which we use as splitter later
    op_str.substitute("-", "%-%");
    op_str.substitute("+", "%+%");
    // split while keeping + and - as separate entries
    op_str.split("%", list);

    // some further sanity check if adduct formula is correct
    String m_part(list[0]);
    // std::cout << m_part.at(m_part.size() - 1) << std::endl;

    if (!m_part.hasSuffix("M"))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "First term of adduct string must contain the molecular entity 'M', optionally prefixed by a multiplier (e.g. '2M'); not found in ", m_part);
    }

    int mol_multiplier(1);
    // check if M has a multiplier in front
    if (m_part.length() > 1)
    { // will throw conversion error of not a number
      mol_multiplier = m_part.prefix(m_part.length()-1).toDouble();
    }

    // evaluate the adduct string ...
    // ... add/subtract each adduct compound
    bool op_plus(false);
    EmpiricalFormula ef; // will remain empty if there are no explicit adducts (e.g. 'M;+1')
    for (Size part_idx = 1 /* omit 0 index, since its 'M' */; part_idx < list.size(); ++part_idx)
    {
      if (list[part_idx] == "+")
      {
        op_plus = true;
        continue;
      }
      else if (list[part_idx] == "-")
      {
        op_plus = false;
        continue;
      }

      // std::cout << "putting " << tmpvec2[part_idx] << " into a formula with mass ";

      // check if formula has got a stoichiometry factor in front
      String formula_str(list[part_idx]);
      int stoichio_factor(1);
      int idx(0);
      while (isdigit(formula_str[idx])) ++idx;
      if (idx > 0)
      {
        stoichio_factor = formula_str.substr(0, idx).toInt();
        formula_str = formula_str.substr(idx, formula_str.size());
      }

      EmpiricalFormula ef_part(formula_str);
      OPENMS_LOG_DEBUG << "Adducts: " << stoichio_factor << "*" << formula_str << " == " << stoichio_factor * ef_part.getMonoWeight() << std::endl;

      if (op_plus)
      {
        ef += ef_part * stoichio_factor;
      }
      else // "-" operator
      {
        ef -= ef_part * stoichio_factor;
      }
    }

    return AdductInfo(cp_str, ef, charge, mol_multiplier);
  }
} // closing namespace OpenMS
