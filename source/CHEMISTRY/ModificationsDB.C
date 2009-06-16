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
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

namespace OpenMS
{
	ModificationsDB::ModificationsDB()
	{
		readFromOBOFile("CHEMISTRY/PSI-MOD.obo");
	}

	ModificationsDB::~ModificationsDB()
	{
		modification_names_.clear();
		for (vector<ResidueModification*>::iterator it = mods_.begin(); it != mods_.end(); ++it)
		{
			delete *it;
		}
	}

	Size ModificationsDB::getNumberOfModifications() const
	{
		return mods_.size();
	}

	const ResidueModification& ModificationsDB::getModification(Size index) const
	{
		if (index >= mods_.size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, mods_.size());
		}
		return *mods_[index];
	}

	
	set<String> ModificationsDB::searchModifications(const String& name) const
	{
		if (!modification_names_.has(name))
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, name);
		}
		
		set<const ResidueModification*> mods = modification_names_[name];
		set<String> new_mods;
		for (set<const ResidueModification*>::const_iterator it = mods.begin(); it != mods.end(); ++it)
		{
			new_mods.insert((*it)->getId());
		}
		return new_mods;
	}

	const ResidueModification& ModificationsDB::getModification(const String& mod_name) const
	{
		set<const ResidueModification*> mods;
		if (modification_names_.has(mod_name))
		{
			mods = modification_names_[mod_name];
		}
		if (mods.size() != 1)
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, mod_name);
		}
		return **mods.begin();
	}
	
	const ResidueModification& ModificationsDB::getModification(const String& residue_name, const String& mod_name) const
	{
		if (ResidueDB::getInstance()->getResidue(residue_name) == 0)
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, residue_name);
		}
		String res = ResidueDB::getInstance()->getResidue(residue_name)->getOneLetterCode();
		set<String> mods = searchModifications(mod_name);
		const ResidueModification* mod_res = 0;
		const ResidueModification* mod_x = 0;
		for (set<String>::const_iterator it = mods.begin(); it != mods.end(); ++it)
		{
			if (res == ModificationsDB::getInstance()->getModification(*it).getOrigin())
			{
				if (mod_res == 0)
				{
					mod_res = &ModificationsDB::getInstance()->getModification(*it);
				}
				else
				{
					cerr << "ModificationsDB::getModification: Warning more than one modification found for modification '" << mod_name << "' at residue '" << residue_name << "', picking first" << endl;
				}
			}
			if ("X" == ModificationsDB::getInstance()->getModification(*it).getOrigin() && mod_x == 0)
			{
				mod_x = &ModificationsDB::getInstance()->getModification(*it);
			}
		}

		if (mod_res == 0 && mod_x == 0)
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, mod_name + "@" + residue_name);
		}
		if (mod_res != 0)
		{
			return *mod_res;
		}
		return *mod_x;
	}

	Size ModificationsDB::findModificationIndex(const String& mod_name) const
	{
		Size idx(0);
		if (modification_names_.has(mod_name))
		{
			if (modification_names_[mod_name].size() > 1)
			{
				throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, "more than one element of name '" + mod_name + "' found!");
			}
		}
		else
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, mod_name);
		}

		const ResidueModification* mod = *modification_names_[mod_name].begin();
		for (Size i = 0; i != mods_.size(); ++i)
		{
			if (mods_[i] == mod)
			{
				idx = i;
				break;
			}
		}
		return idx;
	}
	
  void ModificationsDB::getModificationsByDiffMonoMass(vector<String>& mods, double mass, double error)
	{
		for (vector<ResidueModification*>::const_iterator it = mods_.begin(); it != mods_.end(); ++it)
		{
			if (fabs((*it)->getDiffMonoMass() - mass) <= error)
			{
				mods.push_back((*it)->getId());
			}
		}
	}

	void ModificationsDB::getModificationsByDiffMonoMass(vector<String>& mods, const String& residue, double mass, double error)
	{
		for (vector<ResidueModification*>::const_iterator it = mods_.begin(); it != mods_.end(); ++it)
		{
			if (fabs((*it)->getDiffMonoMass() - mass) <= error)
			{
				String origin = (*it)->getOrigin();
				if (ResidueDB::getInstance()->getResidue(origin) == ResidueDB::getInstance()->getResidue(residue))
				{
					mods.push_back((*it)->getId());
				}
			}
		}

		// no specific mod found? Then use 'X' as origin
		if (mods.size() == 0)
		{
			for (vector<ResidueModification*>::const_iterator it = mods_.begin(); it != mods_.end(); ++it)
			{
				if (fabs((*it)->getDiffMonoMass() - mass) <= error)
				{
					if ((*it)->getOrigin() == "X")
					{
						mods.push_back((*it)->getId());
					}
				}
			}
		}
	}
	
	void ModificationsDB::readFromUnimodXMLFile(const String& filename)
	{
		String file = File::find(filename);

		UnimodXMLFile().load(file, mods_);

    for (vector<ResidueModification*>::const_iterator it = mods_.begin(); it !=mods_.end(); ++it)
    {
      modification_names_[(*it)->getFullName()].insert(*it);
    }

		return;
	}

	void ModificationsDB::readFromOBOFile(const String& filename)
	{
		ResidueModification mod;
		Map<String, ResidueModification> all_mods;
	
 		ifstream is(File::find(filename).c_str());
    String line, line_wo_spaces, id;
      
		//parse file
    while(getline(is, line, '\n'))
    {
    	line.trim();
      line_wo_spaces = line;
      line_wo_spaces.removeWhitespaces();

      if (line == "" || line[0] == '!') //skip empty lines and comments
			{
				continue;
			}
      
			if (line_wo_spaces == "[Term]") //new term
      {
        if (id != "") //store last term
        {
          all_mods[id] = mod;
          id="";
					mod = ResidueModification();
        }
      }
      //new id line
      else if (line_wo_spaces.hasPrefix("id:"))
      {
        id = line.substr(line.find(':')+1).trim();
				mod.setId(id);
      }
      else if (line_wo_spaces.hasPrefix("name:"))
      {
        mod.setFullName(line.substr(line.find(':')+1).trim());
      }
      else if (line_wo_spaces.hasPrefix("is_a:"))
      {
				// TODO
      }
      else if (line_wo_spaces.hasPrefix("def:"))
      {
        // TODO
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
					 Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, line, "missing \" characters to enclose argument!");
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
					Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, line, "missing \" characters to enclose argument!");
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
					mod.setDiffFormula(tmp_split[1]);
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
					mod.setOrigin(val_split[1]);
        }
        else if (val.hasPrefix("Source:"))
        {
					mod.setSourceClassification(val_split[1]);
        }
        else if (val.hasPrefix("TermSpec:"))
        {
					mod.setTermSpecificity(val_split[1]);
        }
			}
    }

		if (id!="") //store last term
		{
			all_mods[id] = mod;
		}

		// now use the term and all synonyms to build the database
		for (Map<String, ResidueModification>::ConstIterator it = all_mods.begin(); it != all_mods.end(); ++it)
		{
			mods_.push_back(new ResidueModification(it->second));
			
			set<String> synonyms = it->second.getSynonyms();
			synonyms.insert(it->first);
			synonyms.insert(it->second.getFullName());

			// now check each of the names and link it to the residue modification
			for (set<String>::const_iterator nit = synonyms.begin(); nit != synonyms.end(); ++nit)
			{
				modification_names_[*nit].insert(mods_.back());
			}
		}

		return;
	}
} // namespace OpenMS

