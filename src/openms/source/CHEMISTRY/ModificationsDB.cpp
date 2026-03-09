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

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CHEMISTRY/ModificationDataProvider.h>
#include <OpenMS/CHEMISTRY/UnimodXMLDataProvider.h>
#include <OpenMS/CHEMISTRY/OBODataProvider.h>

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
    std::vector<std::unique_ptr<ModificationDataProvider>> providers;
    if (!unimod_file.empty())
    {
      providers.push_back(std::make_unique<UnimodXMLDataProvider>(std::move(unimod_file)));
    }
    if (!custommod_file.empty())
    {
      providers.push_back(std::make_unique<UnimodXMLDataProvider>(std::move(custommod_file)));
    }
    if (!psimod_file.empty())
    {
      providers.push_back(std::make_unique<OBODataProvider>(std::move(psimod_file)));
    }
    if (!xlmod_file.empty())
    {
      providers.push_back(std::make_unique<OBODataProvider>(std::move(xlmod_file)));
    }

    static ModificationsDB* db_ = new ModificationsDB(std::move(providers));
    is_instantiated_ = true;
    return db_;
  }

  ModificationsDB::ModificationsDB(std::vector<std::unique_ptr<ModificationDataProvider>> providers)
  {
    loadFromProviders_(providers);
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

  void ModificationsDB::loadFromProviders_(std::vector<std::unique_ptr<ModificationDataProvider>>& providers)
  {
    for (auto& provider : providers)
    {
      auto mods = provider->loadModifications();
      for (auto& m : mods)
      {
        // OBO mods with UniMod record ID: alias resolution
        // Map PSI-MOD accession to existing UniMod entries
        if (m->getUniModRecordId() > 0 && !m->getPSIMODAccession().empty())
        {
          #pragma omp critical(OpenMS_ModificationsDB)
          {
            auto existing = modification_names_.find(m->getUniModAccession());
            if (existing != modification_names_.end())
            {
              for (const auto* existing_mod : existing->second)
              {
                modification_names_[m->getPSIMODAccession()].insert(existing_mod);
              }
            }
          }
          // Whether alias resolved or not, don't add as a new modification.
          // Matches original readFromOBOFile behavior: mods with UniMod record IDs
          // are only used for alias mapping, never added as separate entries.
          continue;
        }

        // Regular modification: add to database and register lookup keys
        #pragma omp critical(OpenMS_ModificationsDB)
        {
          ResidueModification* raw = m.release();
          mods_.push_back(raw);

          // Register standard lookup keys
          modification_names_[raw->getFullId()].insert(raw);
          modification_names_[raw->getId()].insert(raw);
          modification_names_[raw->getFullName()].insert(raw);
          if (!raw->getUniModAccession().empty())
          {
            modification_names_[raw->getUniModAccession()].insert(raw);
          }
          if (!raw->getPSIMODAccession().empty())
          {
            modification_names_[raw->getPSIMODAccession()].insert(raw);
          }

          // Register additional synonyms (OBO providers attach these)
          for (const auto& syn : raw->getSynonyms())
          {
            modification_names_[syn].insert(raw);
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
