// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification2.h>
#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <vector>
#include <fstream>

using namespace std;

namespace OpenMS
{
	ModificationsDB::ModificationsDB()
	{
		// read the modifications from unimod.xml
		//UnimodXMLFile().load("CHEMISTRY/unimod.xml", mods_);
		readFromOBOFile("CHEMISTRY/PSI-MOD.obo");

		//for (vector<ResidueModification2>::const_iterator it = mods_.begin(); it !=mods_.end(); ++it)
		//{
		//	modification_names_[it->getFullName()] = &*it;
		//}
	}

	ModificationsDB::ModificationsDB(const ModificationsDB& /*res_db*/)
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}
	

	ModificationsDB::~ModificationsDB()
	{
		modification_names_.clear();
		for (vector<ResidueModification2*>::iterator it = mods_.begin(); it != mods_.end(); ++it)
		{
			delete *it;
		}
	}

	ModificationsDB& ModificationsDB::operator = (const ModificationsDB& /*res_db*/)
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		return *this;
	}

	UInt ModificationsDB::getNumberOfModifications() const
	{
		return mods_.size();
	}

	const ResidueModification2& ModificationsDB::getModification(UInt index) const
	{
		if (index >= mods_.size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, mods_.size());
		}
		return *mods_[index];
	}

	const ResidueModification2& ModificationsDB::getModification(const String& name) const
	{
		if (!modification_names_.has(name))
		{
			throw Exception::ElementNotFound<String>(__FILE__, __LINE__, __PRETTY_FUNCTION__, name);
		}
		return *modification_names_[name];
		
	}

	void ModificationsDB::readFromUnimodXMLFile(const String& filename)
	{
		UnimodXMLFile().load(filename, mods_);

    for (vector<ResidueModification2*>::const_iterator it = mods_.begin(); it !=mods_.end(); ++it)
    {
      modification_names_[(*it)->getFullName()] = *it;
    }

		return;
	}

	void ModificationsDB::readFromOBOFile(const String& filename)
	{
		ResidueModification2 mod;
		Map<String, ResidueModification2> all_mods;
	
		String file = File::find(filename);
 		ifstream is(file.c_str());
    if (!is)
    {
    	throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

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
        
			if (line_wo_spaces.toLower()=="[term]") //new term
      {
        if (id != "") //store last term
        {
          all_mods[id] = mod;
          id="";
        }
      }
      //new id line
      else if (line_wo_spaces.hasPrefix("id:"))
      {
        id = line.substr(line.find(':')+1).trim();
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
				line_wo_spaces.split('"', val_split);
				if (val_split.size() < 3)
				{
					 Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, line, "missing \" characters to enclose argument!");
				}
				mod.addSynonym(val_split[1]);
			}
			else if (line_wo_spaces.hasPrefix("property_value:"))
			{
				String val = line_wo_spaces.substr(14, line_wo_spaces.size() - 14);
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
		for (Map<String, ResidueModification2>::ConstIterator it = all_mods.begin(); it != all_mods.end(); ++it)
		{
			mods_.push_back(new ResidueModification2(it->second));
			
			set<String> synonyms = it->second.getSynonyms();
			synonyms.insert(it->first);

			// now check each of the names and link it to the residue modification
			for (set<String>::const_iterator nit = synonyms.begin(); nit != synonyms.end(); ++nit)
			{
				if (!modification_names_.has(*nit))
				{
					modification_names_[*nit] = mods_.back();
				}
				else
				{
					cerr << "ModificationsDB: Warning: synonym of Modifications '" << *nit << "'already in use!" << endl;
				}
			}
		}

		return;
	}
} // namespace OpenMS

