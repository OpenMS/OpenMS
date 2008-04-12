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

#include <OpenMS/CHEMISTRY/ResidueModification2.h>

using namespace std;

namespace OpenMS
{

	ResidueModification2::ResidueModification2()
		: allowed_position_(ResidueModification2::ANYWHERE),
			average_mass_(0.0),
			mono_mass_(0.0)
	{
	}

	ResidueModification2::ResidueModification2(const ResidueModification2& rhs)
		: title_(rhs.title_),
			full_name_(rhs.full_name_),
			allowed_position_(rhs.allowed_position_),
			site_(rhs.site_),
			classification_(rhs.classification_),
			average_mass_(rhs.average_mass_),
			mono_mass_(rhs.mono_mass_),
			composition_(rhs.composition_)
	{
	}
	
	ResidueModification2& ResidueModification2::operator = (const ResidueModification2& rhs)
  {
    title_ = rhs.title_;
		full_name_ = rhs.full_name_;
		allowed_position_ = rhs.allowed_position_;
		site_ = rhs.site_;
		classification_ = rhs.classification_;
		average_mass_ = rhs.average_mass_;
		mono_mass_ = rhs.mono_mass_;
		composition_ = rhs.composition_;

		return *this;
  }
	
	/*
	bool ResidueModification2::operator == (const ResidueModification2& rhs)
	{
		return  title_ == rhs.title_ &&
						full_name_ == rhs.full_name_ &&
						allowed_position_ == rhs.allowed_position_ &&
						site_ == rhs.site_ &&
						classification_ == rhs.classification_ &&
						average_mass_ == rhs.average_mass_ &&
						mono_mass_ == rhs.mono_mass_ && 
						composition_ == rhs.composition_;
																											
	}*/
	
	ResidueModification2::~ResidueModification2()
	{

	}

	void ResidueModification2::setTitle(const String& title)
	{
		title_ = title;
	}

	const String& ResidueModification2::getTitle() const
	{
		return title_;
	}

	void ResidueModification2::setFullName(const String& full_name)
	{
		full_name_ = full_name;
	}

	const String& ResidueModification2::getFullName() const
	{
		return full_name_;
	}

	void ResidueModification2::setAllowedPosition(ResidueModification2::AllowedPosition position)
	{
		allowed_position_ = position;
	}
	
	ResidueModification2::AllowedPosition ResidueModification2::getAllowedPosition() const
	{
		return allowed_position_;
	}

	String ResidueModification2::getAllowedPositionName() const
	{
		switch(allowed_position_)
		{
			case ANY_C_TERM: return "Any C-term";
			case ANY_N_TERM: return "Any N-term";
			case PROTEIN_C_TERM: return "Protein C-term";
			case PROTEIN_N_TERM: return "Protein N-term";
			default: // ANYWHERE
				return "Anywhere";
		}
	}
	
	void ResidueModification2::setSite(const String& site)
	{
		site_ = site;
	}

	const String& ResidueModification2::getSite() const
	{
		return site_;
	}

	void ResidueModification2::setClassification(const String& classification)
	{
		classification_ = classification;
	}

	const String& ResidueModification2::getClassification() const
	{
		return classification_;
	}

	void ResidueModification2::setAverageMass(double mass)
	{
		average_mass_ = mass;
	}

	double ResidueModification2::getAverageMass() const
	{
		return average_mass_;
	}

	void ResidueModification2::setMonoMass(double mass)
	{
		mono_mass_ = mass;
	}

	double ResidueModification2::getMonoMass() const
	{
		return mono_mass_;
	}

	void ResidueModification2::setComposition(const String& composition)
	{
		composition_ = composition;
	}

	const String& ResidueModification2::getComposition() const
	{
		return composition_;
	}
/*
	// modification
	ResidueModification2::ResidueModification2()
		: add_average_weight_(0.0f),
			add_mono_weight_(0.0f),
			del_average_weight_(0.0f),
			del_mono_weight_(0.0f)
	{
	}

	ResidueModification2::ResidueModification2(const ResidueModification2& modification)
		:	name_(modification.name_),
			short_name_(modification.short_name_),
			name_prefix_(modification.name_prefix_),
			synonyms_(modification.synonyms_),
			add_formula_(modification.add_formula_),
			add_average_weight_(modification.add_average_weight_),
			add_mono_weight_(modification.add_mono_weight_),
			del_formula_(modification.del_formula_),
			del_average_weight_(modification.del_average_weight_),
			del_mono_weight_(modification.del_mono_weight_),
			valid_residues_(modification.valid_residues_)
	{
	}

	ResidueModification2::~ResidueModification2()
	{
	}

	ResidueModification2& ResidueModification2::operator = (const ResidueModification2& modification)
	{
		if (this != &modification)
		{
			name_ = modification.name_;
			name_prefix_ = modification.name_prefix_;
			short_name_ = modification.short_name_;
			synonyms_ = modification.synonyms_;
			add_formula_ = modification.add_formula_;
			add_average_weight_ = modification.add_average_weight_;
			add_mono_weight_ = modification.add_mono_weight_;
			del_formula_ = modification.del_formula_;
			del_average_weight_ = modification.del_average_weight_;
			del_mono_weight_ = modification.del_mono_weight_;
			valid_residues_ = modification.valid_residues_;
		}
		return *this;
	}
	
	void ResidueModification2::setName(const String& name)
	{
		name_ = name;
	}

	const String& ResidueModification2::getName() const
	{
		return name_;
	}

	void ResidueModification2::setShortName(const String& short_name)
	{
		short_name_ = short_name;
	}

	const String& ResidueModification2::getShortName() const
	{
		return short_name_;
	}

	void ResidueModification2::setNamePrefix(const String& prefix)
	{
		name_prefix_ = prefix;
	}

	const String& ResidueModification2::getNamePrefix() const
	{
		return name_prefix_;
	}

	void ResidueModification2::setSynonyms(const std::set<String>& synonyms)
	{
		synonyms_ = synonyms;
	}

	void ResidueModification2::addSynonym(const String& synonym)
	{
		synonyms_.insert(synonym);
	}

	const std::set<String>& ResidueModification2::getSynonyms() const
	{
		return synonyms_;
	}
	
	void ResidueModification2::setAddFormula(const EmpiricalFormula& formula)
	{
		add_formula_ = formula;
	}

	const EmpiricalFormula& ResidueModification2::getAddFormula() const
	{
		return add_formula_;
	}

	void ResidueModification2::setAddAverageWeight(DoubleReal weight)
	{
		add_average_weight_ = weight;
	}

	DoubleReal ResidueModification2::getAddAverageWeight() const
	{
		return add_average_weight_;
	}

	void ResidueModification2::setAddMonoWeight(DoubleReal weight)
	{
		add_mono_weight_ = weight;
	}

	DoubleReal ResidueModification2::getAddMonoWeight() const
	{
		return add_mono_weight_;
	}

	void ResidueModification2::setDelFormula(const EmpiricalFormula& formula)
	{
		del_formula_ = formula;
	}

	const EmpiricalFormula& ResidueModification2::getDelFormula() const
	{
		return del_formula_;
	}

	void ResidueModification2::setDelAverageWeight(DoubleReal weight)
	{
		del_average_weight_ = weight;
	}

	DoubleReal ResidueModification2::getDelAverageWeight() const
	{
		return del_average_weight_;
	}
	
	void ResidueModification2::setDelMonoWeight(DoubleReal weight)
	{
		del_mono_weight_ = weight;
	}

	DoubleReal ResidueModification2::getDelMonoWeight() const
	{
		return del_mono_weight_;
	}

	void ResidueModification2::setValidResidues(const set<Residue*>& valid_residues)
	{
		valid_residues_ = valid_residues;
	}

	void ResidueModification2::addValidResidue(Residue* valid_residue)
	{
		valid_residues_.insert(valid_residue);
	}
	
	const set<Residue*>& ResidueModification2::getValidResidues() const
	{
		return valid_residues_;
	}

	bool ResidueModification2::operator == (const ResidueModification2& modification) const
	{
		return 	name_ == modification.name_ &&
						name_prefix_ == modification.name_prefix_ &&
						synonyms_ == modification.synonyms_ &&
						add_formula_ == modification.add_formula_ &&
						add_average_weight_ == modification.add_average_weight_ &&
						add_mono_weight_ == modification.add_mono_weight_ &&
						del_formula_ == modification.del_formula_ &&
						del_average_weight_ == modification.del_average_weight_ &&
						del_mono_weight_ == modification.del_mono_weight_ &&
						valid_residues_ == modification.valid_residues_;
	}
	
	bool ResidueModification2::operator != (const ResidueModification2& modification) const
	{
		return !(*this == modification);
	}
*/
}

