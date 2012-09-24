// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Sandro Andreotti $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PepNovoInfile.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <algorithm>
#include <set>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

namespace OpenMS
{

  PepNovoInfile::PepNovoInfile()
  {
  }

  PepNovoInfile::PepNovoInfile(const PepNovoInfile & pepnovo_infile)
  {
    mods_ = pepnovo_infile.mods_;
    mods_and_keys_ = pepnovo_infile.mods_and_keys_;
    ptm_file_ = pepnovo_infile.ptm_file_;
  }

  PepNovoInfile::~PepNovoInfile()
  {
  }

  PepNovoInfile & PepNovoInfile::operator=(const PepNovoInfile & pepnovo_infile)
  {
    if (this != &pepnovo_infile)
    {
      mods_ = pepnovo_infile.mods_;
      mods_and_keys_ = pepnovo_infile.mods_and_keys_;
      ptm_file_ = pepnovo_infile.ptm_file_;
    }
    return *this;
  }

  bool PepNovoInfile::operator==(const PepNovoInfile & pepnovo_infile) const
  {
    if (this != &pepnovo_infile)
    {
      //return ( PTMname_residues_mass_type_ == pepnovo_infile.getModifications() );
      return mods_ == pepnovo_infile.mods_;
    }
    return true;
  }

  String PepNovoInfile::handlePTMs_(const String & modification, const bool variable)
  {
    String locations, key, type;

    ResidueModification::Term_Specificity ts = ModificationsDB::getInstance()->getModification(modification).getTermSpecificity();
    String origin = ModificationsDB::getInstance()->getModification(modification).getOrigin();
    DoubleReal mass = ModificationsDB::getInstance()->getModification(modification).getDiffMonoMass();
    String full_name = ModificationsDB::getInstance()->getModification(modification).getFullName();
    String full_id = ModificationsDB::getInstance()->getModification(modification).getFullId();


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

    default: throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid Term_Specificity", String(ts));
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
    if (origin == "C-term")
    {
      origin = "C_TERM";
    }
    else if (origin == "N-term")
    {
      origin = "N_TERM";
    }
    else
    {
      key = origin;
    }

    if (mass >= 0)
    {
      key += "+" + String(Math::round(mass));
    }
    else
    {
      key += String(Math::round(mass));
    }

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

  void PepNovoInfile::store(const String & filename)
  {
    ptm_file_.store(filename);
  }

  void PepNovoInfile::setModifications(const StringList & fixed_mods, const StringList & variable_mods)
  {
    mods_.setModifications(fixed_mods, variable_mods);
    mods_and_keys_.clear();

    //TextFile ptm_file_;
    ptm_file_.reserve(mods_.getNumberOfModifications() + 1);
    ptm_file_.push_back("#AA\toffset\ttype\tlocations\tsymbol\tPTM\tname");

    // fixed modifications
    std::set<String> fixed_modifications = mods_.getFixedModificationNames();
    for (std::set<String>::const_iterator it = fixed_modifications.begin(); it != fixed_modifications.end(); ++it)
    {
      ptm_file_.push_back(handlePTMs_(*it, false));
    }
    // variable modifications
    std::set<String> variable_modifications = mods_.getVariableModificationNames();
    for (std::set<String>::const_iterator it = variable_modifications.begin(); it != variable_modifications.end(); ++it)
    {
      ptm_file_.push_back(handlePTMs_(*it, true));
    }
  }

  void PepNovoInfile::getModifications(std::map<String, String> & modification_key_map) const
  {
    modification_key_map = mods_and_keys_;
  }

} //namespace
