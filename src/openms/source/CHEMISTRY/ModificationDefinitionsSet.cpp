// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

using namespace std;

namespace OpenMS
{
  ModificationDefinitionsSet::ModificationDefinitionsSet() :
    max_mods_per_peptide_(0)
  {
  }

  ModificationDefinitionsSet::ModificationDefinitionsSet(const ModificationDefinitionsSet& rhs) = default;

  ModificationDefinitionsSet::ModificationDefinitionsSet(const StringList& fixed_modifications, const StringList& variable_modifications) :
    max_mods_per_peptide_(0)
  {
    setModifications(fixed_modifications, variable_modifications);
  }

  ModificationDefinitionsSet::~ModificationDefinitionsSet() = default;

  void ModificationDefinitionsSet::setMaxModifications(Size max_mod)
  {
    max_mods_per_peptide_ = max_mod;
  }

  Size ModificationDefinitionsSet::getMaxModifications() const
  {
    return max_mods_per_peptide_;
  }

  Size ModificationDefinitionsSet::getNumberOfModifications() const
  {
    return variable_mods_.size() + fixed_mods_.size();
  }

  Size ModificationDefinitionsSet::getNumberOfFixedModifications() const
  {
    return fixed_mods_.size();
  }

  Size ModificationDefinitionsSet::getNumberOfVariableModifications() const
  {
    return variable_mods_.size();
  }

  void ModificationDefinitionsSet::addModification(const ModificationDefinition& mod_def)
  {
    if (mod_def.isFixedModification())
    {
      fixed_mods_.insert(mod_def);
    }
    else
    {
      variable_mods_.insert(mod_def);
    }
    return;
  }

  void ModificationDefinitionsSet::setModifications(const set<ModificationDefinition>& mods)
  {
    fixed_mods_.clear();
    variable_mods_.clear();

    for (set<ModificationDefinition>::const_iterator it = mods.begin(); it != mods.end(); ++it)
    {
      if (it->isFixedModification())
      {
        fixed_mods_.insert(*it);
      }
      else
      {
        variable_mods_.insert(*it);
      }
    }
    return;
  }

  void ModificationDefinitionsSet::setModifications(const String& fixed_modifications, const String& variable_modifications)
  {
    setModifications(ListUtils::create<String>(fixed_modifications), ListUtils::create<String>(variable_modifications));
  }

  void ModificationDefinitionsSet::setModifications(const StringList& fixed_modifications, const StringList& variable_modifications)
  {
    fixed_mods_.clear();
    variable_mods_.clear();

    for (StringList::const_iterator it = fixed_modifications.begin(); it != fixed_modifications.end(); ++it)
    {
      ModificationDefinition def(*it, true);
      fixed_mods_.insert(def);
    }

    for (StringList::const_iterator it = variable_modifications.begin(); it != variable_modifications.end(); ++it)
    {
      ModificationDefinition def(*it, false);
      variable_mods_.insert(def);
    }
  }

  set<ModificationDefinition> ModificationDefinitionsSet::getModifications() const
  {
    set<ModificationDefinition> mods = fixed_mods_;
    for (set<ModificationDefinition>::const_iterator it = variable_mods_.begin(); it != variable_mods_.end(); ++it)
    {
      mods.insert(*it);
    }

    return mods;
  }

  set<String> ModificationDefinitionsSet::getModificationNames() const
  {
    set<String> mod_names;
    for (set<ModificationDefinition>::const_iterator it = variable_mods_.begin(); it != variable_mods_.end(); ++it)
    {
      mod_names.insert(it->getModificationName());
    }
    for (set<ModificationDefinition>::const_iterator it = fixed_mods_.begin(); it != fixed_mods_.end(); ++it)
    {
      mod_names.insert(it->getModificationName());
    }
    return mod_names;
  }

  void ModificationDefinitionsSet::getModificationNames(StringList& fixed_modifications, StringList& variable_modifications) const
  {
    fixed_modifications.clear();
    fixed_modifications.reserve(fixed_mods_.size());
    for (set<ModificationDefinition>::const_iterator it = fixed_mods_.begin(); it != fixed_mods_.end(); ++it)
    {
      fixed_modifications.push_back(it->getModificationName());
    }
    variable_modifications.clear();
    variable_modifications.reserve(variable_mods_.size());
    for (set<ModificationDefinition>::const_iterator it = variable_mods_.begin(); it != variable_mods_.end(); ++it)
    {
      variable_modifications.push_back(it->getModificationName());
    }
  }

  const set<ModificationDefinition>& ModificationDefinitionsSet::getFixedModifications() const
  {
    return fixed_mods_;
  }

  const set<ModificationDefinition>& ModificationDefinitionsSet::getVariableModifications() const
  {
    return variable_mods_;
  }

  set<String> ModificationDefinitionsSet::getFixedModificationNames() const
  {
    set<String> mod_names;
    for (set<ModificationDefinition>::const_iterator it = fixed_mods_.begin(); it != fixed_mods_.end(); ++it)
    {
      mod_names.insert(it->getModificationName());
    }
    return mod_names;
  }

  set<String> ModificationDefinitionsSet::getVariableModificationNames() const
  {
    set<String> mod_names;
    for (set<ModificationDefinition>::const_iterator it = variable_mods_.begin(); it != variable_mods_.end(); ++it)
    {
      mod_names.insert(it->getModificationName());
    }
    return mod_names;
  }

  ModificationDefinitionsSet& ModificationDefinitionsSet::operator=(const ModificationDefinitionsSet& rhs)
  {
    if (this != &rhs)
    {
      variable_mods_ = rhs.variable_mods_;
      fixed_mods_ = rhs.fixed_mods_;
      max_mods_per_peptide_ = rhs.max_mods_per_peptide_;
    }
    return *this;
  }

  bool ModificationDefinitionsSet::isCompatible(const AASequence& peptide) const
  {
    set<String> var_names(getVariableModificationNames()), fixed_names(getFixedModificationNames());
    // no modifications present and needed
    if (fixed_names.empty() && !peptide.isModified())
    {
      return true;
    }

    // check whether the fixed modifications are fulfilled
    for (set<String>::const_iterator it1 = fixed_names.begin(); it1 != fixed_names.end(); ++it1)
    {
      String origin = ModificationsDB::getInstance()->getModification(*it1)->getOrigin();
      // only single 1lc amino acids are allowed
      if (origin.size() != 1) continue;
     
      for (AASequence::ConstIterator it2 = peptide.begin(); it2 != peptide.end(); ++it2)
      {
        if (origin == it2->getOneLetterCode())
        {
          // check whether the residue is modified (has to be)
          if (!it2->isModified())
          {
            return false;
          }
          // check whether the modification is the same
          if (ModificationsDB::getInstance()->getModification(*it1)->getId() != it2->getModificationName())
          {
            return false;
          }
        }
      }
    }

    // check whether other modifications than the variable are present
    for (AASequence::ConstIterator it = peptide.begin(); it != peptide.end(); ++it)
    {
      if (it->isModified())
      {
        String mod = it->getModification()->getFullId();
        if (var_names.find(mod) == var_names.end() &&
            fixed_names.find(mod) == fixed_names.end())
        {
          return false;
        }
      }
    }

    if (peptide.hasNTerminalModification())
    {
      String mod = peptide.getNTerminalModification()->getFullId();
      if (var_names.find(mod) == var_names.end() &&
          fixed_names.find(mod) == fixed_names.end())
      {
        return false;
      }
    }

    if (peptide.hasCTerminalModification())
    {
      String mod = peptide.getCTerminalModification()->getFullId();
      if (var_names.find(mod) == var_names.end() &&
          fixed_names.find(mod) == fixed_names.end())
      {
        return false;
      }
    }

    return true;
  }

  bool ModificationDefinitionsSet::operator==(const ModificationDefinitionsSet& rhs) const
  {
    return variable_mods_ == rhs.variable_mods_ &&
           fixed_mods_ == rhs.fixed_mods_ &&
           max_mods_per_peptide_ == rhs.max_mods_per_peptide_;
  }

  bool ModificationDefinitionsSet::operator!=(const ModificationDefinitionsSet& rhs) const
  {
    return !(*this == rhs);
  }

  void ModificationDefinitionsSet::addMatches_(multimap<double, ModificationDefinition>& matches, double mass, const String& residue, ResidueModification::TermSpecificity term_spec, const set<ModificationDefinition>& source, bool is_delta, double tolerance)
  {
    for (set<ModificationDefinition>::const_iterator it = source.begin();
         it != source.end(); ++it)
    {
      const ResidueModification& mod = it->getModification();
      // do the residues match?
      char origin = mod.getOrigin();
      if (!(residue.empty() || (origin == 'X') || (residue[0] == origin) ||
            (residue == ".") || (residue == "X"))) continue;
      // do the term specificities match?
      if (!((term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY) ||
            (term_spec == mod.getTermSpecificity()))) continue;
      // do the masses match?
      double mass_error = tolerance;
      if (is_delta)
      {
        mass_error = fabs(mod.getDiffMonoMass() - mass);
        if (mass_error > tolerance) continue;
      }
      else
      {
        double mod_mass = mod.getMonoMass();
        if ((mod_mass <= 0) && !residue.empty())
        {
          // no absolute mass stored? - calculate it based on the residue
          // (see 'ModificationsDB::getBestModificationByMonoMass'):
          const Residue* res = ResidueDB::getInstance()->getResidue(residue);
          if (res == nullptr) continue;
          double weight = (res->getMonoWeight() - 
                           res->getInternalToFull().getMonoWeight());
          mod_mass = mod.getDiffMonoMass() + weight;
        }
        mass_error = fabs(mod_mass - mass);
        if (mass_error > tolerance) continue;
      }
      matches.insert(make_pair(mass_error, *it));
    }
  }

  void ModificationDefinitionsSet::findMatches(multimap<double, ModificationDefinition>& matches, double mass, const String& residue, ResidueModification::TermSpecificity term_spec, bool consider_fixed, bool consider_variable, bool is_delta, double tolerance) const
  {
    if (!consider_variable && !consider_fixed)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No modifications to consider - set 'consider_variable' and/or 'consider_fixed' to true.");
    }

    matches.clear();
    if (consider_fixed)
    {
      addMatches_(matches, mass, residue, term_spec, fixed_mods_, is_delta, tolerance);
    }
    if (consider_variable)
    {
      addMatches_(matches, mass, residue, term_spec, variable_mods_, is_delta, tolerance);
    }
  }

  // @TODO: should this function handle "max_mods_per_peptide_" as well?
  void ModificationDefinitionsSet::inferFromPeptides(const vector<PeptideIdentification>& peptides)
  {
    // amino acid (or terminus) -> set of modifications (incl. no mod. = 0):
    map<String, set<const ResidueModification*> > mod_map;

    for (const PeptideIdentification& pep : peptides)
    {
      for (const PeptideHit& hit : pep.getHits())
      {
        const AASequence& seq = hit.getSequence();
        mod_map["N-term"].insert(seq.getNTerminalModification());
        mod_map["C-term"].insert(seq.getCTerminalModification());
        for (AASequence::ConstIterator seq_it = seq.begin();
             seq_it != seq.end(); ++seq_it)
        {
          mod_map[seq_it->getOneLetterCode()].insert(seq_it->getModification());
        }
      }
    }

    fixed_mods_.clear();
    variable_mods_.clear();
    for (map<String, set<const ResidueModification*> >::const_iterator map_it =
           mod_map.begin(); map_it != mod_map.end(); ++map_it)
    {
      set<const ResidueModification*>::const_iterator set_it =
        map_it->second.begin();
      // if there's only one mod, it's probably a fixed one:
      if ((map_it->second.size() == 1) && (*set_it != 0))
      {
        ModificationDefinition mod_def(**set_it, true);
        fixed_mods_.insert(mod_def);
      }
      else // variable mod(s)
      {
        for (; set_it != map_it->second.end(); ++set_it)
        {
          if (*set_it != 0)
          {
            ModificationDefinition mod_def(**set_it, false);
            variable_mods_.insert(mod_def);
          }
        }
      }
    }
  }


} // namespace OpenMS
