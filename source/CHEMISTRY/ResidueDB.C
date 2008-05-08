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

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification2.h>
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
		//readResidueModificationsFromFile_("CHEMISTRY/Modifications.xml" );
		//buildResidueModificationNames_();
		//buildModifiedResidues_();
	}

	/*
	ResidueDB::ResidueDB(const String& res_filename, const String& mod_filename) 
		throw (Exception::FileNotFound, Exception::ParseError)
	{
		readResiduesFromFile_(res_filename);
		buildResidueNames_();
		//readResidueModificationsFromFile_(mod_filename);
		//buildResidueModificationNames_();
		//buildModifiedResidues_();
	}
	*/
	

	//ResidueDB::ResidueDB(const ResidueDB& /*res_db*/)
	//{
	//	throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	//}
	

	ResidueDB::~ResidueDB()
	{
		clear_();
	}


	//ResidueDB& ResidueDB::operator = (const ResidueDB& /*res_db*/)
	//{
	//	throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	//	return *this;
	//}

	const Residue* ResidueDB::getResidue(const String& name) const
	{
		if (residue_names_.has(name))
		{
			return residue_names_[name];
		}
		return 0;
	}

	UInt ResidueDB::getNumberOfResidues() const
	{
		return residues_.size();
	}
		
	/*
	const ResidueModification* ResidueDB::getModification(const String& name) const
	{
		if (modification_names_.has(name))
		{
			return modification_names_[name];
		}
		return 0;
	}

	set<const ResidueModification*> ResidueDB::getModifications(const Residue* residue) const
	{
		return getModifications(residue->getName());
	}

	set<const ResidueModification*> ResidueDB::getModifications(const String& res_name) const
	{
		set <const ResidueModification*> mods;
		/// @improvement speed up computation! (Andreas)
		set<ResidueModification*>::iterator mit = modifications_.begin();
		for (;mit!=modifications_.end();++mit)
		{
			const set<Residue*> res = (*mit)->getValidResidues();
			set<Residue*>::const_iterator rit = res.begin();
			for (;rit!=res.end();++rit)
			{
				if (residue_names_.has((*rit)->getUnmodifiedName()) &&
						residue_names_.has(res_name))
				{
					if (residue_names_[(*rit)->getUnmodifiedName()] == residue_names_[res_name])
					{
						mods.insert(*mit);
					}
				}
			}
		}
		return mods;
	}
	*/
	
	/*
	const set<const ResidueModification*>& ResidueDB::getModifications() const
	{
		return const_modifications_;
	}

	set<const Residue*> ResidueDB::getResidues(const ResidueModification* modification) const
	{
		set<const Residue*> res;
		set<const Residue*>::const_iterator it;
		for (it=const_residues_.begin();it!=const_residues_.end();++it)
		{
			if ((*it)->isModified() && (*it)->getModification() == modification)
			{
				res.insert(*it);
			}
		}
		return res;
	}

	set<const Residue*> ResidueDB::getResidues(const String& mod_name) const
	{
		set<const Residue*> res;
		if (modification_names_.has(mod_name))
		{
			const ResidueModification* mod = modification_names_[mod_name];
			set<const Residue*>::const_iterator it;
			for (it=const_residues_.begin();it!=const_residues_.end();++it)
			{
				if ((*it)->isModified() && (*it)->getModification() == mod)
				{
					res.insert(*it);
				}
			}
		}
		return res;
	}
	*/

	const set<const Residue*>& ResidueDB::getResidues() const
	{
		return const_residues_;
	}

	/*
	UInt ResidueDB::getNumberOfResidueModifications() const
	{
		return modifications_.size();
	}*/

	/*
	void ResidueDB::setModifications(const String& file_name) 
		throw(Exception::FileNotFound, Exception::ParseError)
	{
		clearResidueModifications_();
		readResidueModificationsFromFile_(file_name);
		buildResidueModificationNames_();
		buildModifiedResidues_();
	}*/
	
	/*
	void ResidueDB::addResidueModification(ResidueModification modification )
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);		
	}
*/
	
	void ResidueDB::setResidues(const String& file_name)
		throw(Exception::FileNotFound, Exception::ParseError)
	{
		clearResidues_();
		readResiduesFromFile_(file_name);
		buildResidueNames_();
		//buildModifiedResidues_();
	}

	void ResidueDB::addResidue(const Residue& residue)
	{
		Residue * r = new Residue(residue);
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
			const ResidueModification2& mod = ModificationsDB::getInstance()->getModification(r->getModification());

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
		return;
	}

	/*
	bool ResidueDB::hasResidueModification(const String& mod_name) const
	{
		if (modification_names_.has(mod_name))
		{
			return true;
		}
		return false;
	}
*/
	bool ResidueDB::hasResidue(const String& res_name) const
	{
		if (residue_names_.has(res_name))
		{
			return true;
		}
		return false;
	}
	
	void ResidueDB::readResiduesFromFile_(const String& file_name)
		throw(Exception::FileNotFound, Exception::ParseError)
	{
		String file = File::find(file_name);

		if (!File::exists(file))
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, file_name);
		}
		
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
		catch (...)
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", "");
		}
	}

	/*
	void ResidueDB::readResidueModificationsFromFile_(const String& filename)
		throw(Exception::FileNotFound, Exception::ParseError)
	{
		String file = File::find(filename);
		// try filename
		if (!File::exists(file))
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}

		// seems to be ok, open it
		Param param;
		param.load(file);

		if (!param.begin().getName().hasPrefix("Modifications"))
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", "");
		}
		
		try
		{
			vector<String> split;
			param.begin().getName().split(':',split);
			String prefix = split[0] + split[1];
			ResidueModification* mod_ptr = new ResidueModification();

			modifications_.insert(mod_ptr);
			const_modifications_.insert(mod_ptr);

			Map<String, Map<String, String> > valid_res;
			
			for (Param::ParamIterator it=param.begin(); it!=param.end(); ++it)
			{
				it.getName().split(':',split);
				if (prefix != split[0] + split[1])
				{
					prefix = split[0] + split[1];
					Map<String, Map<String, String> >::iterator hit;
					for (hit = valid_res.begin(); hit != valid_res.end(); ++hit)
					{
						mod_ptr->addValidResidue(parseResidue_(hit->second));
					}
					valid_res.clear();
					
					mod_ptr = new ResidueModification();
					modifications_.insert(mod_ptr);
					const_modifications_.insert(mod_ptr);
				}

				String value = it->value;
				String key(split[2]);

				if (key == "Name")
				{
					mod_ptr->setName(value);
					continue;
				}
				if (key == "ShortName")
				{
					mod_ptr->setShortName(value);
					continue;
				}
				if (key == "NamePrefix")
				{
					mod_ptr->setNamePrefix(value);
					continue;
				}
				if (key == "DelFormula")
				{
					EmpiricalFormula del(value);
					mod_ptr->setDelFormula(del);
					mod_ptr->setDelAverageWeight(del.getAverageWeight());
					mod_ptr->setDelMonoWeight(del.getAverageWeight());
					continue;
				}
				if (key == "AddFormula")
				{
					EmpiricalFormula add(value);
					mod_ptr->setAddFormula(value);
					mod_ptr->setAddAverageWeight(add.getAverageWeight());
					mod_ptr->setAddMonoWeight(add.getMonoWeight());
					continue;
				}
				if (key == "ValidResidues")
				{
					valid_res[String(split[3])][it.getName()] = value;
					continue;
				}
				if (key == "Synonyms")
				{
					// no synonyms defined?
					if (!key.hasSuffix(":"))
					{
						mod_ptr->addSynonym(value);
					}
					continue;
				}
				cerr << "unknown key: " << key << ", with value: " << value << endl;
			}
			
			// add valid residues of the last modification
			Map<String, Map<String, String> >::iterator hit;
			for (hit = valid_res.begin(); hit != valid_res.end(); ++hit)
			{
				mod_ptr->addValidResidue(parseResidue_(hit->second));
			}
		}
		catch (...)
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "" , "");
		}
	}
*/

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

	/*
	void ResidueDB::clearResidueModifications_()
	{
		set<ResidueModification*>::iterator it;
		for (it=modifications_.begin();it!=modifications_.end();++it)
		{
			delete *it;
		}
		modifications_.clear();
		modification_names_.clear();
		const_modifications_.clear();

		// delete modified residues as well
		set<Residue*>::iterator rit;
		set<Residue*> to_delete;
		for (rit=residues_.begin();rit!=residues_.end();++rit)
		{
			if ((*rit)->isModified())
			{
				to_delete.insert(*rit);
			}
		}

		for (rit=to_delete.begin();rit!=to_delete.end();++rit)
		{
			delete *rit;
			residues_.erase(*rit);
			const_residues_.erase(*rit);
		}
		residue_names_.clear();
		buildResidueNames_();
	}
*/
	
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
			if (key.hasSuffix(":LossName"))
			{
				res_ptr->setLossName(value);
				continue;
			}
			if (key.hasSuffix(":LossFormula"))
			{
				EmpiricalFormula loss(value);
				res_ptr->setLossFormula(loss);
				res_ptr->setLossAverageWeight(loss.getAverageWeight());
				res_ptr->setLossMonoWeight(loss.getMonoWeight());
				continue;
			}
			/*
			if (key.hasSuffix(":UnmodifiedName"))
			{
				res_ptr->setUnmodifiedName(value);
				continue;
			}
			*/
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
				res_ptr->setSideChainBasicity(value.toFloat());
				continue;
			}
			if (key.hasSubstring("GB_BB_L"))
			{
				res_ptr->setBackboneBasicityLeft(value.toFloat());
				continue;
			}
			if (key.hasSubstring("GB_BB_R"))
			{
				res_ptr->setBackboneBasicityRight(value.toFloat());
				continue;
			}
			cerr << "unknown key: " << key << ", with value: " << value << endl;
		}
		
		if (low_mass_ions.size() != 0)
		{
			res_ptr->setLowMassIons(low_mass_ions);
		}
		
		return res_ptr;
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

	/*
	void ResidueDB::buildResidueModificationNames_()
	{
		set<ResidueModification*>::iterator it;
		for (it=modifications_.begin(); it!=modifications_.end(); ++it)
		{
			if ((*it)->getName() != "")
			{
				modification_names_[(*it)->getName()] = *it;
			}
			if ((*it)->getNamePrefix() != "")
			{
				modification_names_[(*it)->getNamePrefix()] = *it;
			}
			if ((*it)->getShortName() != "")
			{
				modification_names_[(*it)->getShortName()] = *it;
			}
			set<String>::iterator sit;
			set<String> syn = (*it)->getSynonyms();
			for (sit = syn.begin(); sit!=syn.end(); ++sit)
			{
				if (*sit != "")
				{
					modification_names_[*sit] = *it;
				}
			}
		}
	}

	void ResidueDB::buildModifiedResidues_()
	{
		set<ResidueModification*>::iterator mit = modifications_.begin();
		for (;mit!=modifications_.end();++mit)
		{
			set<Residue*> residues = (*mit)->getValidResidues();
			set<Residue*>::iterator rit = residues.begin();
			for (;rit!=residues.end();++rit)
			{
				const Residue * orig = getResidue((*rit)->getUnmodifiedName());
				if (orig == 0)
				{
					cerr << "cannot find residue for modification: " << (*rit)->getUnmodifiedName();
				}
				else
				{
					Residue * r = new Residue(**rit);
					ResidueModification * m = *mit;
					r->setModification(m);
					r->setAverageWeight(orig->getAverageWeight()+m->getAddAverageWeight()-m->getDelAverageWeight());
					r->setMonoWeight(orig->getMonoWeight()+m->getAddMonoWeight()-m->getDelMonoWeight());
					r->setFormula(orig->getFormula()+m->getAddFormula()-m->getDelFormula());
					residues_.insert(r);
					const_residues_.insert(r);
				}
			}
		}
	}
*/
	bool ResidueDB::operator == (const ResidueDB& /*rhs*/) const
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		return false;
	}

	bool ResidueDB::operator != (const ResidueDB& /*rhs*/) const
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		return false;
	}

	const Residue* ResidueDB::getModifiedResidue(const String& modification)
	{
		const ResidueModification2& mod = ModificationsDB::getInstance()->getModification(modification);
		
		String origin = mod.getOrigin();
		if (origin.size() == 1)
		{
			if (!residue_names_.has(origin))
			{
				throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Residue with name "
							+ origin + " was not registered in residue DB, register first; needed for modification '" + modification + "'!").c_str());
			}
		}
		else
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Origin of a modification is only allowed to" 
															" have exactly one amino acid given, instead of '" + origin  + "' of modification '" + modification + "'!").c_str());
		}

		return getModifiedResidue(getResidue(origin), modification);
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
		
		if (residue_mod_names_.has(res_name) && residue_mod_names_[res_name].has(modification))
		{
			return residue_mod_names_[res_name][modification];
		}

		
		Residue* res = new Residue(*residue_names_[res_name]);
		res->setModification(modification);
		
		// now register this modified residue 
		addResidue_(res);
		return res;
	}
}

