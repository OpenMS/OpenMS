// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/SYSTEM/File.h>

using namespace std;

namespace OpenMS
{
	ResidueDB::ResidueDB()
	{
		readResiduesFromFile_("CHEMISTRY/Residues.xml" );
		buildResidueNames_();
	}

	ResidueDB::~ResidueDB()
	{
		clear_();
	}


	const Residue* ResidueDB::getResidue(const String& name) const
	{
		if (residue_names_.has(name))
		{
			return residue_names_[name];
		}
		return 0;
	}

	Size ResidueDB::getNumberOfResidues() const
	{
		return residues_.size();
	}

	Size ResidueDB::getNumberOfModifiedResidues() const
	{
		return modified_residues_.size();
	}
		
	const set<const Residue*> ResidueDB::getResidues(const String& residue_set) const
	{
		if (!residues_by_set_.has(residue_set))
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Residue set cannot be found: '" + residue_set + "'");
		}
		
		return residues_by_set_[residue_set];
	}

	void ResidueDB::setResidues(const String& file_name)
	{
		clearResidues_();
		readResiduesFromFile_(file_name);
		buildResidueNames_();
	}

	void ResidueDB::addResidue(const Residue& residue)
	{
		Residue* r = new Residue(residue);
		addResidue_(r);
	}

	void ResidueDB::addResidue_(Residue* r)
	{
		vector<String> names;
		if (r->getName() != "")
		{
			names.push_back(r->getName());
		}
		if (r->getShortName() != "")
		{
			names.push_back(r->getShortName());
		}
		set<String> synonyms = r->getSynonyms();
    for (set<String>::iterator it = synonyms.begin(); it != synonyms.end(); ++it)
    {
    	names.push_back(*it);
    }

		
		if (!r->isModified())
		{
			for (vector<String>::const_iterator it = names.begin(); it != names.end(); ++it)
			{
				residue_names_[*it] = r;
			}
			residues_.insert(r);
			const_residues_.insert(r);
		}
		else
		{
			modified_residues_.insert(r);
			const_modified_residues_.insert(r);

			// get all modification names
			vector<String> mod_names;
			const ResidueModification& mod = ModificationsDB::getInstance()->getModification(r->getOneLetterCode(), r->getModification(), ResidueModification::ANYWHERE);

			mod_names.push_back(mod.getId());
			mod_names.push_back(mod.getFullName());
			set<String> mod_synonyms = mod.getSynonyms();
			for (set<String>::iterator it = mod_synonyms.begin(); it != mod_synonyms.end(); ++it)
			{
				mod_names.push_back(*it);
			}

			for (vector<String>::const_iterator it = names.begin(); it != names.end(); ++it)
			{
				for (vector<String>::const_iterator mod_it = mod_names.begin(); mod_it != mod_names.end(); ++mod_it)
				{
					residue_mod_names_[*it][*mod_it] = r;
				}
			}
		}
		buildResidueNames_();
		return;
	}

	bool ResidueDB::hasResidue(const String& res_name) const
	{
		if (residue_names_.has(res_name))
		{
			return true;
		}
		return false;
	}

	bool ResidueDB::hasResidue(const Residue* residue) const
	{
		if (const_residues_.find(residue) != const_residues_.end() ||
				const_modified_residues_.find(residue) != const_modified_residues_.end())
		{
			return true;
		}
		return false;
	}
	
	void ResidueDB::readResiduesFromFile_(const String& file_name)
	{
		String file = File::find(file_name);
		
		Param param;
		param.load(file);
		
		if (!param.begin().getName().hasPrefix("Residues"))
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", "");
		}

		try
		{
			vector<String> split;
			param.begin().getName().split(':',split);
			String prefix = split[0] + split[1];
			Residue * res_ptr = 0;
		
			Map<String, String> values;
			
			for (Param::ParamIterator it=param.begin(); it!=param.end(); ++it)
			{
				it.getName().split(':',split);
				if (prefix != split[0] + split[1])
				{
					// add residue
					res_ptr = parseResidue_(values);
					values.clear();
					residues_.insert(res_ptr);
					const_residues_.insert(res_ptr);
					prefix = split[0] + split[1];
				}
				
				String value = it->value;
				String key = it.getName();
				values[key] = value;

			}

			// add last residue
			res_ptr = parseResidue_(values);
			residues_.insert(res_ptr);
			const_residues_.insert(res_ptr);
		}
		catch (Exception::BaseException& e)
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, e.what(), "");
		}
	}

	void ResidueDB::clear_()
	{
		clearResidues_();
		//clearResidueModifications_();
	}

	void ResidueDB::clearResidues_()
	{
		set<Residue*>::iterator it;
		for (it=residues_.begin();it!=residues_.end();++it)
		{
			delete *it;
		}
		residues_.clear();
		residue_names_.clear();
		const_residues_.clear();
	}

	Residue* ResidueDB::parseResidue_(Map<String, String>& values) 
	{
		vector<EmpiricalFormula> low_mass_ions;
		Residue * res_ptr = new Residue();
		
		for (Map<String, String>::iterator it=values.begin();it!=values.end();++it)
		{
			String key(it->first);
			String value(it->second);

			if (key.hasSuffix(":Name"))
			{
				res_ptr->setName(value);
				continue;
			}
			if (key.hasSuffix(":ShortName"))
			{
				res_ptr->setShortName(value);
				continue;
			}
			if (key.hasSuffix(":ThreeLetterCode"))
			{
				res_ptr->setThreeLetterCode(value);
				continue;
			}
			if (key.hasSuffix(":OneLetterCode"))
			{
				res_ptr->setOneLetterCode(value);
				continue;
			}
			if (key.hasSuffix(":Formula"))
			{
				EmpiricalFormula formula(value);
				res_ptr->setFormula(value);
				res_ptr->setAverageWeight(formula.getAverageWeight());
				res_ptr->setMonoWeight(formula.getMonoWeight());
				continue;
			}

			if (key.hasSubstring(":Losses:LossName"))
			{
				res_ptr->addLossName(value);
				continue;
			}
			if (key.hasSubstring(":Losses:LossFormula"))
			{
				EmpiricalFormula loss(value);
				res_ptr->addLossFormula(loss);
				continue;
			}

			if (key.hasSubstring("NTermLosses:LossName"))
			{
				res_ptr->addNTermLossName(value);
				continue;
			}

			if (key.hasSubstring("NTermLosses:LossFormula"))
			{
				EmpiricalFormula loss(value);
				res_ptr->addNTermLossFormula(loss);
				continue;
			}
			
			if (key.hasSubstring("LowMassIons"))
			{
				// no markers defined?
				if (!key.hasSuffix(":"))
				{
					low_mass_ions.push_back(value);
				}
				continue;
			}
			if (key.hasSubstring("Synonyms"))
			{
				// no synonyms defined?
				if (!key.hasSuffix(":"))
				{
					res_ptr->addSynonym(value);
				}
				continue;
			}
			if (key.hasSubstring("pka"))
			{
				// no pka defined?
				if (!key.hasSuffix(":"))
				{
					res_ptr->setPka(value.toDouble());
				}
				continue;
			}
			if (key.hasSubstring("pkb"))
			{
				// no pkb defined?
				if (!key.hasSuffix(":"))
				{
					res_ptr->setPkb(value.toDouble());
				}
				continue;
			}
			if (key.hasSubstring("pkc"))
			{
				// no pkc defined?
				if (!key.hasSuffix(":"))
				{
					res_ptr->setPkc(value.toDouble());
				}
				continue;
			}
			if (key.hasSubstring("GB_SC"))
			{
				res_ptr->setSideChainBasicity(value.toDouble());
				continue;
			}
			if (key.hasSubstring("GB_BB_L"))
			{
				res_ptr->setBackboneBasicityLeft(value.toDouble());
				continue;
			}
			if (key.hasSubstring("GB_BB_R"))
			{
				res_ptr->setBackboneBasicityRight(value.toDouble());
				continue;
			}
			if (key.hasSubstring("ResidueSets"))
			{
				StringList residue_sets = StringList::create(value);
				for (StringList::const_iterator it = residue_sets.begin(); it != residue_sets.end(); ++it)
				{
					res_ptr->addResidueSet(*it);
					residue_sets_.insert(*it);
				}
				continue;
			}
			cerr << "unknown key: " << key << ", with value: " << value << endl;
		}
		
		if (low_mass_ions.size() != 0)
		{
			res_ptr->setLowMassIons(low_mass_ions);
		}

		for (set<String>::const_iterator it = res_ptr->getResidueSets().begin(); it != res_ptr->getResidueSets().end(); ++it)
		{
			residues_by_set_[*it].insert(res_ptr);
		}
		
		return res_ptr;
	}

	const set<String>& ResidueDB::getResidueSets() const
	{
		return residue_sets_;
	}
	
	void ResidueDB::buildResidueNames_()
	{
		set<Residue*>::iterator it;
		for (it = residues_.begin(); it!=residues_.end(); ++it)
		{
			residue_names_[(*it)->getName()] = *it;
			if ((*it)->getThreeLetterCode() != "")
			{
				residue_names_[(*it)->getThreeLetterCode()] = *it;
			}
			if ((*it)->getOneLetterCode() != "")
			{
				residue_names_[(*it)->getOneLetterCode()] = *it;
			}
			if ((*it)->getShortName() != "")
			{
				residue_names_[(*it)->getOneLetterCode()] = *it;
			}
			set<String>::iterator sit;
			set<String> syn = (*it)->getSynonyms();
			for (sit = syn.begin(); sit!=syn.end(); ++sit)
			{
				if (*sit != "")
				{
					residue_names_[*sit] = *it;
				}
			}
		}
	}


	const Residue* ResidueDB::getModifiedResidue(const String& modification)
	{
		const ResidueModification& mod = ModificationsDB::getInstance()->getModification(modification);
		return getModifiedResidue(getResidue(mod.getOrigin()), mod.getFullId());
	}
	
	
	const Residue* ResidueDB::getModifiedResidue(const Residue* residue, const String& modification)
	{
		// search if the mod already exists
		String res_name(residue->getName());

		if (!residue_names_.has(res_name))
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Residue with name " 
											+ res_name + " was not registered in residue DB, register first!").c_str());
		}

		String id = ModificationsDB::getInstance()->getModification(res_name, modification, ResidueModification::ANYWHERE).getId();
		
		if (residue_mod_names_.has(res_name) && residue_mod_names_[res_name].has(id))
		{
			return residue_mod_names_[res_name][id];
		}

		
		Residue* res = new Residue(*residue_names_[res_name]);
		res->setModification(id);
		//res->setLossFormulas(vector<EmpiricalFormula>());
		//res->setLossNames(vector<String>());

		// now register this modified residue 
		addResidue_(res);
		return res;
	}
}

