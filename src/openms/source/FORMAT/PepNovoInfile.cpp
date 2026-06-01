// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PepNovoInfile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <fstream>

using namespace std;

namespace OpenMS
{

  PepNovoInfile::PepNovoInfile() = default;

  PepNovoInfile::PepNovoInfile(const PepNovoInfile& pepnovo_infile)
  {
    mods_ = pepnovo_infile.mods_;
    mods_and_keys_ = pepnovo_infile.mods_and_keys_;
    ptm_file_ = pepnovo_infile.ptm_file_;
  }

  PepNovoInfile::~PepNovoInfile() = default;

  PepNovoInfile& PepNovoInfile::operator=(const PepNovoInfile& pepnovo_infile)
  {
    if (this != &pepnovo_infile)
    {
      mods_ = pepnovo_infile.mods_;
      mods_and_keys_ = pepnovo_infile.mods_and_keys_;
      ptm_file_ = pepnovo_infile.ptm_file_;
    }
    return *this;
  }

  bool PepNovoInfile::operator==(const PepNovoInfile& pepnovo_infile) const
  {
    if (this != &pepnovo_infile)
    {
      //return ( PTMname_residues_mass_type_ == pepnovo_infile.getModifications() );
      return mods_ == pepnovo_infile.mods_;
    }
    return true;
  }

  String PepNovoInfile::handlePTMs_(const String& modification, const bool variable)
  {
    String locations, key, type;

    ResidueModification::TermSpecificity ts = ModificationsDB::getInstance()->getModification(modification)->getTermSpecificity();
    String origin = ModificationsDB::getInstance()->getModification(modification)->getOrigin();
    double mass = ModificationsDB::getInstance()->getModification(modification)->getDiffMonoMass();
    String full_name = ModificationsDB::getInstance()->getModification(modification)->getFullName();
    String full_id = ModificationsDB::getInstance()->getModification(modification)->getFullId();


    if (variable)
    {
      type = "OPTIONAL";
    }
    else
    {
      type = "FIXED";
    }

    switch (ts)
    {
    case ResidueModification::C_TERM: locations = "C_TERM";
      break;

    case ResidueModification::N_TERM: locations = "N_TERM";
      break;

    case ResidueModification::ANYWHERE: locations = "ALL";
      break;

    default: throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid term specificity", String(ts));
    }

    if (ts == ResidueModification::C_TERM)
    {
      key = "$";
    }
    else if (ts == ResidueModification::N_TERM)
    {
      key = "^";
    }

    //cout<<"origin: "<<origin<<"    loc: "<<locations<<endl;
    if ((ts == ResidueModification::C_TERM) && (origin == "X"))
    {
      origin = "C_TERM";
    }
    else if ((ts == ResidueModification::N_TERM) && (origin == "X"))
    {
      origin = "N_TERM";
    }
    else
    {
      key = origin;
    }

    if (mass >= 0)
    {
      key += "+";
    }
    key += String(int(Math::round(mass)));


    String line = "";
    line += origin.toUpper();
    line += "\t";
    line += mass;
    line += "\t";
    line += type;
    line += "\t";
    line += locations;
    line += "\t";
    line += key;
    line += "\t";
    line += full_name;

    mods_and_keys_[key] = full_id;

    return line;
  }

  void PepNovoInfile::store(const String& filename)
  {
    ptm_file_.store(filename);
  }

  void PepNovoInfile::setModifications(const StringList& fixed_mods, const StringList& variable_mods)
  {
    mods_.setModifications(fixed_mods, variable_mods);
    mods_and_keys_.clear();

    //TextFile ptm_file_;
    ptm_file_.addLine("#AA\toffset\ttype\tlocations\tsymbol\tPTM\tname");

    // fixed modifications
    std::set<String> fixed_modifications = mods_.getFixedModificationNames();
    for (std::set<String>::const_iterator it = fixed_modifications.begin(); it != fixed_modifications.end(); ++it)
    {
      ptm_file_.addLine(handlePTMs_(*it, false));
    }
    // variable modifications
    std::set<String> variable_modifications = mods_.getVariableModificationNames();
    for (std::set<String>::const_iterator it = variable_modifications.begin(); it != variable_modifications.end(); ++it)
    {
      ptm_file_.addLine(handlePTMs_(*it, true));
    }
  }

  void PepNovoInfile::getModifications(std::map<String, String>& modification_key_map) const
  {
    modification_key_map = mods_and_keys_;
  }

} //namespace
