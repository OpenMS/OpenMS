// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <fstream>
#include <limits>
#include <utility>

using namespace std;

namespace OpenMS
{

  bool ModificationsDB::residuesMatch_(const char residue, const ResidueModification* curr_mod) const
  {
    const char origin = curr_mod->getOrigin();

    if (origin != 'X')
    {
      // residues match if they are equal or they match everything (X/.)
      return (origin == residue || residue == 'X' || residue == '.' || residue == '?');
    }
    else
    {
      // origin is X, this usually means that the modification can be at any amino acid

      // residues do NOT match if the modification is user-defined and has origin
      // X (which here means an actual input AA X and it does *not* mean "match
      // all AA") while the current residue is not X. Make sure we don't match things like
      // PEPN[400] and PEPX[400] since these have very different masses.
      bool non_matching_user_defined = (
           curr_mod->isUserDefined() &&
           residue != '?' &&
           origin != residue );

      return !non_matching_user_defined;
    }
  }

  bool ModificationsDB::is_instantiated_ = false;

  ModificationsDB* ModificationsDB::getInstance()
  {
    static ModificationsDB* db_ = ModificationsDB::initializeModificationsDB();
    return db_;
  }

  ModificationsDB* ModificationsDB::initializeModificationsDB(OpenMS::String unimod_file, OpenMS::String custommod_file, OpenMS::String psimod_file, OpenMS::String xlmod_file)
  {
    // Currently its not possible to check for double initialization since getInstance() also calls this function.
    // if (is_instantiated_)
    // {
    //   throw Exception::FailedAPICall(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot initialize ModificationsDB twice");
    // }

    static ModificationsDB* db_ = new ModificationsDB(std::move(unimod_file), std::move(custommod_file), std::move(psimod_file), std::move(xlmod_file));
    return db_;
  }

  ModificationsDB::ModificationsDB(const OpenMS::String& unimod_file, const OpenMS::String& custommod_file, const OpenMS::String& psimod_file, const OpenMS::String& xlmod_file)
  {
    if (!unimod_file.empty())
    {
      readFromUnimodXMLFile(unimod_file);
    }

    if(!custommod_file.empty())
    {
      readFromUnimodXMLFile(custommod_file); 
    }

    if (!psimod_file.empty())
    {
      readFromOBOFile(psimod_file);
    }

    if (!xlmod_file.empty())
    {
      readFromOBOFile(xlmod_file);
    }
    is_instantiated_ = true;
  }

  ModificationsDB::~ModificationsDB()
  {
    modification_names_.clear();
    for (auto it = mods_.begin(); it != mods_.end(); ++it)
    {
      delete *it;
    }
  }

  bool ModificationsDB::isInstantiated()
  {
    return is_instantiated_;
  }

  Size ModificationsDB::getNumberOfModifications() const
  {
    Size s;
    #pragma omp critical (OpenMS_ModificationsDB)
    {
      s = mods_.size();
    }
    return s;
  }

  const ResidueModification* ModificationsDB::searchModificationsFast(const String& mod_name_,
                                                                      bool& multiple_matches,
                                                                      const String& residue,
                                                                      ResidueModification::TermSpecificity term_spec
                                                                      ) const
  {
    const ResidueModification* mod(nullptr);

    String mod_name = mod_name_;
    multiple_matches = false;

    char res = '?'; // empty
    if (!residue.empty()) res = residue[0];

    #pragma omp critical(OpenMS_ModificationsDB)
    {
      bool found = true;
      auto modifications = modification_names_.find(mod_name);
      if (modifications == modification_names_.end())
      {
        // Try to fix things, Skyline for example uses unimod:10 and not UniMod:10 syntax
        if (mod_name.size() > 6 && mod_name.prefix(6).toLower() == "unimod")
        {
          mod_name = "UniMod" + mod_name.substr(6, mod_name.size() - 6);
        }

        modifications = modification_names_.find(mod_name);
        if (modifications == modification_names_.end())
        {
          OPENMS_LOG_WARN << OPENMS_PRETTY_FUNCTION << "Modification not found: " << mod_name << endl;
          found = false; 
        }
      }

      int nr_mods = 0;
      if (found)
      {
        for (const auto& it : modifications->second)
        {
          if ( residuesMatch_(res, it) &&
               (term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY ||
               (term_spec == it->getTermSpecificity())))
          {
            mod = it;
            nr_mods++;
          }
        }
      }
      if (nr_mods > 1) multiple_matches = true;
    }
    return mod;
  }

  const ResidueModification* ModificationsDB::searchModification(const ResidueModification& mod_in) const
  {
    const ResidueModification* mod(nullptr);

    const String& mod_name = mod_in.getFullId();

    #pragma omp critical(OpenMS_ModificationsDB)
    {
      bool found = true;
      auto modifications = modification_names_.find(mod_name);

      if (modifications == modification_names_.end())
      {
        OPENMS_LOG_WARN << OPENMS_PRETTY_FUNCTION << "Modification not found: " << mod_name << endl;
        found = false; 
      }

      if (found)
      {
        for (const auto& mod_indb : modifications->second)
        {
          if (mod_in == *mod_indb)
          {
            mod = mod_indb;
            break;
          }
        }
      }
    }
    return mod;
  }

  const ResidueModification* ModificationsDB::getModification(Size index) const
  {
    OPENMS_PRECONDITION(index < mods_.size(), "Index out of bounds in ModificationsDB::getModification(Size index)." );
    return mods_[index];
  }

  void ModificationsDB::searchModifications(set<const ResidueModification*>& mods,
                                            const String& mod_name_,
                                            const String& residue,
                                            ResidueModification::TermSpecificity term_spec) const
  {
    mods.clear();

    String mod_name = mod_name_;

    char res = '?'; // empty
    if (!residue.empty()) res = residue[0];

    #pragma omp critical(OpenMS_ModificationsDB)
    {
      bool found = true;
      auto modifications = modification_names_.find(mod_name);
      if (modifications == modification_names_.end())
      {
        // Try to fix things, Skyline for example uses unimod:10 and not UniMod:10 syntax
        if (mod_name.size() > 6 && mod_name.prefix(6).toLower() == "unimod")
        {
          mod_name = "UniMod" + mod_name.substr(6, mod_name.size() - 6);
        }

        modifications = modification_names_.find(mod_name);
        if (modifications == modification_names_.end())
        {
          OPENMS_LOG_WARN << OPENMS_PRETTY_FUNCTION << "Modification not found: " << mod_name << endl;
          found = false; 
        }
      }

      if (found)
      {
        for (const auto& it : modifications->second)
        {
          if ( residuesMatch_(res, it) &&
               (term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY ||
               (term_spec == it->getTermSpecificity())))
          {
            mods.insert(it);
          }
        }
      }
    } 
  }

  const ResidueModification* ModificationsDB::getModification(const String& mod_name, const String& residue, ResidueModification::TermSpecificity term_spec) const
  {
    const ResidueModification* mod(nullptr);
    // if residue is specified, try residue-specific search first to avoid
    // ambiguities (e.g. "Carbamidomethyl (N-term)"/"Carbamidomethyl (C)"):
    bool multiple_matches = false;
    if (!residue.empty() &&
        (term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY))
    {
      mod = searchModificationsFast(mod_name, multiple_matches, residue,
                          ResidueModification::ANYWHERE);
    }
    if (mod == nullptr) mod = searchModificationsFast(mod_name, multiple_matches, residue, term_spec);
    if (mod == nullptr)
    {
      String message = String("Retrieving the modification failed. It is not available for the residue '") + residue 
        + "' and term specificity '" + ResidueModification().getTermSpecificityName(term_spec) + "'. ";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, mod_name);
    }
    if (multiple_matches)
    {
      OPENMS_LOG_WARN << "Warning (ModificationsDB::getModification): more than one modification with name '" + mod_name + "', residue '" + residue + "', specificity '" + String(Int(term_spec)) << "' found, picking the first one only.";
      // for (auto it = mods.begin(); it != mods.end(); ++it)
      // {
      //   OPENMS_LOG_WARN << " " << (*it)->getFullId();
      // }
      OPENMS_LOG_WARN << "\n";
    }
    return mod;
  }


  bool ModificationsDB::has(const String & modification) const
  {
    bool has_mod;
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      has_mod = (modification_names_.find(modification) != modification_names_.end());
    }
    return has_mod;
  }

  Size ModificationsDB::findModificationIndex(const String & mod_name) const
  {
    if (!has(mod_name))
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Modification not found: " + mod_name);
    }

    bool one_mod(true);
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      if (modification_names_.find(mod_name)->second.size() > 1)
      {
        one_mod = false;
      }
    }
    if (!one_mod) 
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "More than one modification with name: " + mod_name);
    }

    Size index(numeric_limits<Size>::max());
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      const ResidueModification* mod = *(modification_names_.find(mod_name)->second.begin());
      for (Size i = 0; i != mods_.size(); ++i)
      {
        if (mods_[i] == mod)
        {
          index = i;
          break;
        }
      }
    }

    if (index == numeric_limits<Size>::max())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Modification name found but modification not found: " + mod_name);
    }
    return index;
  }

  void ModificationsDB::searchModificationsByDiffMonoMass(vector<String>& mods, double mass, double max_error, const String& residue, ResidueModification::TermSpecificity term_spec)
  {
    mods.clear();
    char res = '?'; // empty
    if (!residue.empty()) res = residue[0];
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      for (auto const & m : mods_)
      {
        if ((fabs(m->getDiffMonoMass() - mass) <= max_error) &&
            residuesMatch_(res, m) &&
            ((term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY) ||
             (term_spec == m->getTermSpecificity())))
        {
          mods.push_back(m->getFullId());
        }
      }
    }
  }

  void ModificationsDB::searchModificationsByDiffMonoMass(vector<const ResidueModification*>& mods, double mass, double max_error, const String& residue, ResidueModification::TermSpecificity term_spec)
  {
    mods.clear();
    char res = '?'; // empty
    if (!residue.empty()) res = residue[0];
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      for (auto const & m : mods_)
      {
        if ((fabs(m->getDiffMonoMass() - mass) <= max_error) &&
            residuesMatch_(res, m) &&
            ((term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY) ||
             (term_spec == m->getTermSpecificity())))
        {
          mods.push_back(m);
        }
      }
    }
  }

  void ModificationsDB::searchModificationsByDiffMonoMassSorted(vector<String>& mods, double mass, double max_error, const String& residue, ResidueModification::TermSpecificity term_spec)
  {
    mods.clear();
    std::map<std::pair<double,Size>, const String&> diff_idx2mods;
    char res = '?'; // empty
    if (!residue.empty()) res = residue[0];
    double diff = 0;
    Size cnt = 0;
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      for (auto const & m : mods_)
      {
        diff = fabs(m->getDiffMonoMass() - mass);
        if ((diff <= max_error) &&
            residuesMatch_(res, m) &&
            ((term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY) ||
             (term_spec == m->getTermSpecificity())))
        {
          diff_idx2mods.emplace(make_pair(diff, cnt++), m->getFullId());
        }
      }
    }
    for (const auto& foo_mod : diff_idx2mods)
    {
      mods.push_back(foo_mod.second);
    }
  }

  void ModificationsDB::searchModificationsByDiffMonoMassSorted(vector<const ResidueModification*>& mods, double mass, double max_error, const String& residue, ResidueModification::TermSpecificity term_spec)
  {
    mods.clear();
    std::map<std::pair<double,Size>, const ResidueModification*> diff_idx2mods;
    char res = '?'; // empty
    if (!residue.empty()) res = residue[0];
    double diff = 0;
    Size cnt = 0;
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      for (auto const & m : mods_)
      {
        diff = fabs(m->getDiffMonoMass() - mass);
        if ((diff <= max_error) &&
            residuesMatch_(res, m) &&
            ((term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY) ||
             (term_spec == m->getTermSpecificity())))
        {
          diff_idx2mods.emplace(make_pair(diff, cnt++), m);
        }
      }
    }
    for (const auto& foo_mod : diff_idx2mods)
    {
      mods.push_back(foo_mod.second);
    }
  }


  const ResidueModification* ModificationsDB::getBestModificationByDiffMonoMass(double mass, double max_error, const String& residue, ResidueModification::TermSpecificity term_spec)
  {
    double min_error = max_error;
    const ResidueModification* mod = nullptr;
    char res = '?'; // empty
    if (!residue.empty())
    {
      res = residue[0];
    }
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      for (auto const & m : mods_)
      {
        // using less instead of less-or-equal will pick the first matching
        // modification of equally heavy modifications (in our case this is the
        // first matching UniMod entry)
        double mass_error = fabs(m->getDiffMonoMass() - mass);
        if ((mass_error < min_error) &&
            residuesMatch_(res, m) &&
            ((term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY) ||
             (term_spec == m->getTermSpecificity())))
        {
          min_error = mass_error;
          mod = m;
        }
      }
    }
    return mod;
  }

  void ModificationsDB::readFromUnimodXMLFile(const String& filename)
  {
    vector<ResidueModification*> new_mods;
    UnimodXMLFile().load(filename, new_mods);

    for (auto & m : new_mods)
    {
      // create full ID based on other information:
      m->setFullId();

      #pragma omp critical(OpenMS_ModificationsDB)
      {
        // e.g. Oxidation (M)
        modification_names_[m->getFullId()].insert(m);
        // e.g. Oxidation
        modification_names_[m->getId()].insert(m);
        // e.g. Oxidized
        modification_names_[m->getFullName()].insert(m);
        // e.g. UniMod:312
        modification_names_[m->getUniModAccession()].insert(m);
        mods_.push_back(m);
      }
    }
  }

  const ResidueModification* ModificationsDB::addModification(std::unique_ptr<ResidueModification> new_mod)
  {
    const ResidueModification* ret;
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      auto it = modification_names_.find(new_mod->getFullId());
      if (it != modification_names_.end())
      {
        OPENMS_LOG_WARN << "Modification already exists in ModificationsDB. Skipping." << new_mod->getFullId() << endl;
        ret = *(it->second.begin()); // returning from omp critical is not allowed
      }
      else
      {
        modification_names_[new_mod->getFullId()].insert(new_mod.get());
        modification_names_[new_mod->getId()].insert(new_mod.get());
        modification_names_[new_mod->getFullName()].insert(new_mod.get());
        modification_names_[new_mod->getUniModAccession()].insert(new_mod.get());
        mods_.push_back(new_mod.get());
        new_mod.release(); // do not delete the object;
        ret = mods_.back();
      }
    }
    return ret;
  }

  const ResidueModification* ModificationsDB::addModification(const ResidueModification& new_mod)
  {
    const ResidueModification* ret = new ResidueModification(new_mod);
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      auto it = modification_names_.find(new_mod.getFullId());
      if (it != modification_names_.end())
      {
        OPENMS_LOG_WARN << "Modification already exists in ModificationsDB. Skipping." << new_mod.getFullId() << endl;
        ret = *(it->second.begin()); // returning from omp critical is not allowed
      }
      else
      {
        modification_names_[ret->getFullId()].insert(ret);
        modification_names_[ret->getId()].insert(ret);
        modification_names_[ret->getFullName()].insert(ret);
        modification_names_[ret->getUniModAccession()].insert(ret);
        mods_.push_back(const_cast<ResidueModification*>(ret));
        ret = mods_.back();
      }
    }
    return ret;
  }

  const ResidueModification* ModificationsDB::addNewModification_(const ResidueModification& new_mod)
  {
    const ResidueModification* ret = new ResidueModification(new_mod);
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      modification_names_[ret->getFullId()].insert(ret);
      modification_names_[ret->getId()].insert(ret);
      modification_names_[ret->getFullName()].insert(ret);
      modification_names_[ret->getUniModAccession()].insert(ret);
      mods_.push_back(const_cast<ResidueModification*>(ret));
      ret = mods_.back();
    }
    return ret;
  }

  void ModificationsDB::readFromOBOFile(const String& filename)
  {
    ResidueModification mod;
    // add multiple mods for multiple specificities
    //Map<String, ResidueModification> all_mods;
    multimap<String, ResidueModification> all_mods;

    ifstream is(File::find(filename).c_str());
    String line, line_wo_spaces, id;
    String origin = "";

    bool reading_cross_link = false;

    //parse file
    while (getline(is, line, '\n'))
    {
      line.trim();
      line_wo_spaces = line;
      line_wo_spaces.removeWhitespaces();

      if (line.empty() || line[0] == '!') //skip empty lines and comments
      {
        continue;
      }

      if (line_wo_spaces == "[Term]")       //new term
      {
        // if the last [Term] was a moon-link, then it does not belong in CrossLinksDB
        if (!id.empty() && !reading_cross_link) //store last term
        {
          // split into single residues and make unique (for XL-MS, where equal specificities for both sides are possible)
          vector<String> origins;
          origin.split(",", origins);

          std::sort(origins.begin(), origins.end());
          vector<String>::iterator unique_end = unique(origins.begin(), origins.end());
          origins.resize(distance(origins.begin(), unique_end));

          for (vector<String>::iterator orig_it = origins.begin(); orig_it != origins.end(); ++orig_it)
          {
            // we don't allow modifications with ambiguity codes as origin (except "X"):
            if ((orig_it->size() == 1) && (*orig_it != "B") && (*orig_it != "J") && (*orig_it != "Z"))
            {
              mod.setOrigin((*orig_it)[0]);
              all_mods.insert(make_pair(id, mod));
            }
          }

          // for mono-links from XLMOD.obo:
          if (origin.hasSubstring("ProteinN-term"))
          {
            mod.setTermSpecificity(ResidueModification::PROTEIN_N_TERM);
            mod.setOrigin('X');
            all_mods.insert(make_pair(id, mod));
          }
          if (origin.hasSubstring("ProteinC-term"))
          {
            mod.setTermSpecificity(ResidueModification::PROTEIN_C_TERM);
            mod.setOrigin('X');
            all_mods.insert(make_pair(id, mod));
          }

          id = "";
          origin = "";
          mod = ResidueModification();
        }
        else if (reading_cross_link) // re-initialize before reading next [Term]
        {
          id = "";
          origin = "";
          mod = ResidueModification();
          reading_cross_link = false;
        }
      }

      //new id line
      else if (line_wo_spaces.hasPrefix("id:"))
      {
        id = line.substr(line.find(':') + 1).trim();
        mod.setId(id);
        mod.setPSIMODAccession(id);
      }
      else if (line_wo_spaces.hasPrefix("name:"))
      {
        String name = line.substr(line.find(':') + 1).trim();
        mod.setFullName(name);
        if (mod.getId().hasSubstring("XLMOD"))
        {
          mod.setName(name);
          mod.setId(name);
          mod.setFullName(name);
        }
      }
      else if (line_wo_spaces.hasPrefix("is_a:"))
      {
        // TODO
      }
      else if (line_wo_spaces.hasPrefix("def:"))
      {
        line.remove('[');
        line.remove(']');
        line.remove(',');
        vector<String> split;
        line.split(' ', split);
        for (Size i = 0; i != split.size(); ++i)
        {
          if (split[i].hasPrefix("UniMod:"))
          {
            // Parse UniMod identifier to int
            String identifier = split[i].substr(7, split[i].size());
            mod.setUniModRecordId(identifier.toInt());
          }
        }
      }
      else if (line_wo_spaces.hasPrefix("comment:"))
      {
        // TODO
      }
      else if (line_wo_spaces.hasPrefix("synonym:"))
      {
        vector<String> val_split;
        line.split('"', val_split);
        if (val_split.size() < 3)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "missing \" characters to enclose argument!");
        }
        mod.addSynonym(val_split[1]);

        if (line_wo_spaces.hasSubstring("PSI-MOD-label"))
        {
          mod.setName(val_split[1]);
        }
      }
      else if (line_wo_spaces.hasPrefix("property_value:"))
      {
        String val = line_wo_spaces.substr(15, line_wo_spaces.size() - 15);
        val.trim();

        if (val.hasSubstring("\"none\""))
        {
          continue;
        }

        vector<String> val_split;
        val.split('"', val_split);
        if (val_split.size() != 3)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "missing \" characters to enclose argument!");
        }
        if (val.hasPrefix("DiffAvg:"))
        {
          mod.setDiffAverageMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("DiffFormula:"))
        {
          vector<String> tmp_split;
          line.split('"', tmp_split);
          tmp_split[1].removeWhitespaces();
          mod.setDiffFormula(EmpiricalFormula(tmp_split[1]));
        }
        else if (val.hasPrefix("DiffMono:"))
        {
          mod.setDiffMonoMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("Formula:"))
        {
          mod.setFormula(val_split[1]);
        }
        else if (val.hasPrefix("MassAvg:"))
        {
          mod.setAverageMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("MassMono:"))
        {
          mod.setMonoMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("Origin:"))
        {
          //mod.setOrigin(val_split[1]);
          origin = val_split[1];
        }
        else if (val.hasPrefix("Source:"))
        {
          mod.setSourceClassification(val_split[1]);
        }
        else if (val.hasPrefix("TermSpec:"))
        {
          mod.setTermSpecificity(val_split[1]);
        }
        // XLMOD specific fields
        else if (val.hasPrefix("reactionSites:"))
        {
          if (val_split[1] == "2")
          {
            reading_cross_link = true;
          }
        }
        else if (val.hasPrefix("monoisotopicMass:"))
        {
          mod.setDiffMonoMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("specificities:"))
        {
          // TODO cross-linker specificities can be different for both chain sides, right now the union of both sides is used
          // Input parameters of the cross-link search tool make sure, that the chemistry is not violated
          origin = val_split[1];

          // remove brackets
          origin.remove('(');
          origin.remove(')');
          origin.substitute("&", ",");
        }
      }
    }

    if (!id.empty()) //store last term
    {
      // split into single residues and make unique (for XL-MS, where equal specificities for both sides are possible)
      vector<String> origins;
      origin.split(",", origins);

      std::sort(origins.begin(), origins.end());
      vector<String>::iterator unique_end = unique(origins.begin(), origins.end());
      origins.resize(distance(origins.begin(), unique_end));

      for (vector<String>::iterator orig_it = origins.begin(); orig_it != origins.end(); ++orig_it)
      {
        // we don't allow modifications with ambiguity codes as origin (except "X"):
        if ((orig_it->size() == 1) && (*orig_it != "B") && (*orig_it != "J") && (*orig_it != "Z"))
        {
          mod.setOrigin((*orig_it)[0]);
          all_mods.insert(make_pair(id, mod));
        }
      }

      // for mono-links from XLMOD.obo:
      if (origin.hasSubstring("ProteinN-term"))
      {
        mod.setTermSpecificity(ResidueModification::N_TERM);
        mod.setOrigin('X');
        all_mods.insert(make_pair(id, mod));
      }
      if (origin.hasSubstring("ProteinC-term"))
      {
        mod.setTermSpecificity(ResidueModification::C_TERM);
        mod.setOrigin('X');
        all_mods.insert(make_pair(id, mod));
      }

      id = "";
      origin = "";
      mod = ResidueModification();
    }

    // now use the term and all synonyms to build the database
    #pragma omp critical(OpenMS_ModificationsDB)
    {
      for (multimap<String, ResidueModification>::const_iterator it = all_mods.begin(); it != all_mods.end(); ++it)
      {
        // check whether a unimod definition already exists, then simply add synonyms to it
        if (it->second.getUniModRecordId() > 0)
        {
          //cerr << "Found UniMod PSI-MOD mapping: " << it->second.getPSIMODAccession() << " " << it->second.getUniModAccession() << endl;
          set<const ResidueModification*> mods = modification_names_[it->second.getUniModAccession()];
          for (set<const ResidueModification*>::const_iterator mit = mods.begin(); mit != mods.end(); ++mit)
          {
            //cerr << "Adding PSIMOD accession: " << it->second.getPSIMODAccession() << " " << it->second.getUniModAccession() << endl;
            modification_names_[it->second.getPSIMODAccession()].insert(*mit);
          }
        }
        else
        {
          // the mod has so far not been mapped to a unimod mod
          // first check whether the mod is specific
          if ((it->second.getOrigin() != 'X') ||
             ((it->second.getTermSpecificity() != ResidueModification::ANYWHERE) &&
             (it->second.getDiffMonoMass() != 0)))
          {
            mods_.push_back(new ResidueModification(it->second));

            set<String> synonyms = it->second.getSynonyms();
            synonyms.insert(it->first);
            synonyms.insert(it->second.getFullName());
            //synonyms.insert(it->second.getUniModAccession());
            synonyms.insert(it->second.getPSIMODAccession());
            // full ID is auto-generated based on (short) ID, but we want the name instead:
            mods_.back()->setId(it->second.getFullName());
            mods_.back()->setFullId();
            mods_.back()->setId(it->second.getId());
            synonyms.insert(mods_.back()->getFullId());

            // now check each of the names and link it to the residue modification
            for (set<String>::const_iterator nit = synonyms.begin(); nit != synonyms.end(); ++nit)
            {
              modification_names_[*nit].insert(mods_.back());
            }
          }
        }
      }
    }
  }

  void ModificationsDB::getAllSearchModifications(vector<String>& modifications) const
  {
    modifications.clear();

    #pragma omp critical(OpenMS_ModificationsDB)
    {
      for (auto const & m : mods_)
      {
        if (m->getUniModRecordId() > 0)
        {
          modifications.push_back(m->getFullId());
        }
      }
    }

    // sort by name (case INsensitive)
    sort(modifications.begin(), modifications.end(), [&](const String& a, const String& b) {
      size_t i(0);
      while (i < a.size() && i < b.size())
      {
        if (tolower(a[i]) == tolower(b[i]))
        {
          ++i;
        }
        else
        {
          return tolower(a[i]) < tolower(b[i]);
        }
      }
      return a.size() < b.size();
    });
  }

  void ModificationsDB::writeTSV(String const& filename)
  {
    std::ofstream ofs(filename, std::ofstream::out);
    ofs << "FullId\tFullName\tUnimodAccession\tOrigin/AA\tTerminusSpecificity\tDiffMonoMass\n";
    ResidueModification tmp;
    for (const auto& mod : mods_)
    {
      ofs << mod->getFullId() << "\t" << mod->getFullName() << "\t" << mod->getUniModAccession() << "\t" << mod->getOrigin() << "\t"
      << tmp.getTermSpecificityName(mod->getTermSpecificity()) << "\t"
      << mod->getDiffMonoMass() << "\n";
    }
  }
} // namespace OpenMS
