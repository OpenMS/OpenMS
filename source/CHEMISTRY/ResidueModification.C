// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CHEMISTRY/ResidueModification.h>

using namespace std;

namespace OpenMS
{
	// modification
	ResidueModification::ResidueModification()
		: add_average_weight_(0.0f),
			add_mono_weight_(0.0f),
			del_average_weight_(0.0f),
			del_mono_weight_(0.0f)
	{
	}

	ResidueModification::ResidueModification(const ResidueModification& modification)
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

	ResidueModification::~ResidueModification()
	{
	}

	ResidueModification& ResidueModification::operator = (const ResidueModification& modification)
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
	
	void ResidueModification::setName(const String& name)
	{
		name_ = name;
	}

	const String& ResidueModification::getName() const
	{
		return name_;
	}

	void ResidueModification::setShortName(const String& short_name)
	{
		short_name_ = short_name;
	}

	const String& ResidueModification::getShortName() const
	{
		return short_name_;
	}

	void ResidueModification::setNamePrefix(const String& prefix)
	{
		name_prefix_ = prefix;
	}

	const String& ResidueModification::getNamePrefix() const
	{
		return name_prefix_;
	}

	void ResidueModification::setSynonyms(const std::set<String>& synonyms)
	{
		synonyms_ = synonyms;
	}

	void ResidueModification::addSynonym(const String& synonym)
	{
		synonyms_.insert(synonym);
	}

	const std::set<String>& ResidueModification::getSynonyms() const
	{
		return synonyms_;
	}
	
	void ResidueModification::setAddFormula(const EmpiricalFormula& formula)
	{
		add_formula_ = formula;
	}

	const EmpiricalFormula& ResidueModification::getAddFormula() const
	{
		return add_formula_;
	}

	void ResidueModification::setAddAverageWeight(DoubleReal weight)
	{
		add_average_weight_ = weight;
	}

	DoubleReal ResidueModification::getAddAverageWeight() const
	{
		return add_average_weight_;
	}

	void ResidueModification::setAddMonoWeight(DoubleReal weight)
	{
		add_mono_weight_ = weight;
	}

	DoubleReal ResidueModification::getAddMonoWeight() const
	{
		return add_mono_weight_;
	}

	void ResidueModification::setDelFormula(const EmpiricalFormula& formula)
	{
		del_formula_ = formula;
	}

	const EmpiricalFormula& ResidueModification::getDelFormula() const
	{
		return del_formula_;
	}

	void ResidueModification::setDelAverageWeight(DoubleReal weight)
	{
		del_average_weight_ = weight;
	}

	DoubleReal ResidueModification::getDelAverageWeight() const
	{
		return del_average_weight_;
	}
	
	void ResidueModification::setDelMonoWeight(DoubleReal weight)
	{
		del_mono_weight_ = weight;
	}

	DoubleReal ResidueModification::getDelMonoWeight() const
	{
		return del_mono_weight_;
	}

	void ResidueModification::setValidResidues(const set<Residue*>& valid_residues)
	{
		valid_residues_ = valid_residues;
	}

	void ResidueModification::addValidResidue(Residue* valid_residue)
	{
		valid_residues_.insert(valid_residue);
	}
	
	const set<Residue*>& ResidueModification::getValidResidues() const
	{
		return valid_residues_;
	}

	bool ResidueModification::operator == (const ResidueModification& modification) const
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
	
	bool ResidueModification::operator != (const ResidueModification& modification) const
	{
		return !(*this == modification);
	}

}

