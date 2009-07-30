// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

using namespace std;

namespace OpenMS
{
	ModificationDefinitionsSet::ModificationDefinitionsSet()
		: max_mods_per_peptide_(0)
	{
	}

	ModificationDefinitionsSet::ModificationDefinitionsSet(const ModificationDefinitionsSet& rhs)
		: variable_mods_(rhs.variable_mods_),
			fixed_mods_(rhs.fixed_mods_),
			max_mods_per_peptide_(rhs.max_mods_per_peptide_)
	{
	}

	ModificationDefinitionsSet::ModificationDefinitionsSet(const String& fixed_modifications, const String& variable_modifications)
		: max_mods_per_peptide_(0)
	{
		setModifications(fixed_modifications, variable_modifications);
	}

	ModificationDefinitionsSet::ModificationDefinitionsSet(const StringList& fixed_modifications, const StringList& variable_modifications)
		: max_mods_per_peptide_(0)
	{
		setModifications(fixed_modifications, variable_modifications);
	}

	ModificationDefinitionsSet::~ModificationDefinitionsSet()
	{
	}

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
		setModifications(StringList::create(fixed_modifications), StringList::create(variable_modifications));
	}


	void ModificationDefinitionsSet::setModifications(const StringList& fixed_modifications, const StringList& variable_modifications)
	{
		fixed_mods_.clear();
		variable_mods_.clear();

		if (fixed_modifications.size() != 0)
		{
			for (StringList::const_iterator it = fixed_modifications.begin(); it != fixed_modifications.end(); ++it)
			{
				ModificationDefinition def;
				def.setModification(*it);
				def.setFixedModification(true);
				fixed_mods_.insert(def);
			}
		}

		if (variable_modifications.size() != 0)
		{
			for (StringList::const_iterator it = variable_modifications.begin(); it != variable_modifications.end(); ++it)
			{
				ModificationDefinition def;
				def.setModification(*it);
				def.setFixedModification(false);
				variable_mods_.insert(def);
			}
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
			mod_names.insert(it->getModification());
		}
		for (set<ModificationDefinition>::const_iterator it = fixed_mods_.begin(); it != fixed_mods_.end(); ++it)
		{
			mod_names.insert(it->getModification());
		}
		return mod_names;
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
			mod_names.insert(it->getModification());
		}
		return mod_names;
	}

	set<String> ModificationDefinitionsSet::getVariableModificationNames() const
  {
    set<String> mod_names;
    for (set<ModificationDefinition>::const_iterator it = variable_mods_.begin(); it != variable_mods_.end(); ++it)
    {
      mod_names.insert(it->getModification());
    }
    return mod_names;
  }

	ModificationDefinitionsSet& ModificationDefinitionsSet::operator = (const ModificationDefinitionsSet& rhs)
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
		if (fixed_names.size() == 0 && !peptide.isModified())
		{
			return true;
		}

		// check whether the fixed modifications are fulfilled
		if (fixed_names.size() != 0)
		{
			for (set<String>::const_iterator it1 = fixed_names.begin(); it1 != fixed_names.end(); ++it1)
			{
				String origin = ModificationsDB::getInstance()->getModification(*it1).getOrigin();
				// only single 1lc amino acids are allowed
				if (origin.size() != 1)
				{
					continue;
				}
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
						if (*it1 != it2->getModification())
						{
							return false;
						}
					}
				}
			}
		}
		
		// check wether other modifications than the variable are present
		for (AASequence::ConstIterator it = peptide.begin(); it != peptide.end(); ++it)
		{
			if (it->isModified())
			{
				String mod = it->getModification();
				if (var_names.find(mod) == var_names.end() &&
						fixed_names.find(mod) == fixed_names.end())
				{
					return false;
				}
			}
		}

		if (peptide.hasNTerminalModification())
		{
			String mod = peptide.getNTerminalModification();
			if (var_names.find(mod) == var_names.end() &&
					fixed_names.find(mod) == fixed_names.end())
			{
				return false;
			}
		}

		if (peptide.hasCTerminalModification())
		{
			String mod = peptide.getCTerminalModification();
			if (var_names.find(mod) == var_names.end() &&
					fixed_names.find(mod) == fixed_names.end())
			{
				return false;
			}
		}

		return true;
	}

	bool ModificationDefinitionsSet::operator == (const ModificationDefinitionsSet& rhs) const
	{
		return variable_mods_ == rhs.variable_mods_ &&
					 fixed_mods_ == rhs.fixed_mods_ &&
					 max_mods_per_peptide_ == rhs.max_mods_per_peptide_;
	}

	bool ModificationDefinitionsSet::operator != (const ModificationDefinitionsSet& rhs) const
	{
		return !(*this == rhs);
	}

} // namespace OpenMS

