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
		fixed_mods_.clear();
		variable_mods_.clear();

		if (fixed_modifications != "")
		{
			vector<String> split;
			fixed_modifications.split(',', split);
			if (split.size() == 0)
			{
				split.push_back(fixed_modifications);
			}
			for (vector<String>::const_iterator it = split.begin(); it != split.end(); ++it)
			{
				ModificationDefinition def;
				def.setModification(*it);
				def.setFixedModification(true);
				fixed_mods_.insert(def);
			}
		}
		if (variable_modifications != "")
		{
			vector<String> split;
			variable_modifications.split(',', split);
			if (split.size() == 0)
			{
				split.push_back(variable_modifications);
			}
			for (vector<String>::const_iterator it = split.begin(); it != split.end(); ++it)
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

