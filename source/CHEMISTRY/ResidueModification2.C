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
		: term_spec_(ResidueModification2::ANYWHERE),
			average_mass_(0.0),
			mono_mass_(0.0),
			diff_average_mass_(0.0),
			diff_mono_mass_(0.0)
	{
	}

	ResidueModification2::ResidueModification2(const ResidueModification2& rhs)
		: title_(rhs.title_),
			full_name_(rhs.full_name_),
			term_spec_(rhs.term_spec_),
			origin_(rhs.origin_),
			classification_(rhs.classification_),
			average_mass_(rhs.average_mass_),
			mono_mass_(rhs.mono_mass_),
			diff_average_mass_(rhs.diff_average_mass_),
			diff_mono_mass_(rhs.diff_mono_mass_),
			formula_(rhs.formula_),
			valid_residues_(rhs.valid_residues_)
	{
	}
	
	ResidueModification2& ResidueModification2::operator = (const ResidueModification2& rhs)
  {
    title_ = rhs.title_;
		full_name_ = rhs.full_name_;
		term_spec_ = rhs.term_spec_;
		origin_ = rhs.origin_;
		classification_ = rhs.classification_;
		average_mass_ = rhs.average_mass_;
		mono_mass_ = rhs.mono_mass_;
		diff_average_mass_ = rhs.diff_average_mass_;
		diff_mono_mass_ = rhs.diff_mono_mass_;
		formula_ = rhs.formula_;
		valid_residues_ = rhs.valid_residues_;
		
		return *this;
  }
	
	bool ResidueModification2::operator == (const ResidueModification2& rhs) const
	{
		return  title_ == rhs.title_ &&
						full_name_ == rhs.full_name_ &&
						term_spec_ == rhs.term_spec_ &&
						origin_ == rhs.origin_ &&
						classification_ == rhs.classification_ &&
						average_mass_ == rhs.average_mass_ &&
						mono_mass_ == rhs.mono_mass_ && 
						diff_average_mass_ == rhs.diff_average_mass_ &&
						diff_mono_mass_ == rhs.diff_mono_mass_ &&
						formula_ == rhs.formula_ &&
						valid_residues_ == rhs.valid_residues_;
																											
	}
	
	bool ResidueModification2::operator != (const ResidueModification2& rhs) const
	{
		return !(*this == rhs);
	}
	
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

	void ResidueModification2::setTermSpecificity(Term_Specificity term_spec)
	{
		term_spec_ = term_spec;
	}

	void ResidueModification2::setTermSpecificity(const String& term_spec)
	{
		// TODO
		return;
	}
	
	ResidueModification2::Term_Specificity ResidueModification2::getTermSpecificity() const
	{
		return term_spec_;
	}

	String ResidueModification2::getTermSpecitificityName(Term_Specificity term_spec) const
	{
		if (term_spec == NUMBER_OF_TERM_SPECIFICITY)
		{
			term_spec = term_spec_;	
		}
		switch(term_spec)
		{
			case C_TERM: return "Any C-term";
			case N_TERM: return "Any N-term";
			default: // ANYWHERE
				return "Anywhere";
		}
	}
	
	void ResidueModification2::setOrigin(const String& origin)
	{
		origin_ = origin;
	}

	const String& ResidueModification2::getOrigin() const
	{
		return origin_;
	}

	void ResidueModification2::setSourceClassification(Source_Classification classification)
	{
		classification_ = classification;
	}

	void ResidueModification2::setSourceClassification(const String& classification)
	{
		// TODO
	}
	
	ResidueModification2::Source_Classification ResidueModification2::getSourceClassification() const
	{
		return classification_;
	}

	String ResidueModification2::getSourceClassificationName(Source_Classification classification) const
	{
		// TODO
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


	void ResidueModification2::setDiffAverageMass(double mass)
	{
		diff_average_mass_ = mass;
	}

	double ResidueModification2::getDiffAverageMass() const
	{
		return diff_average_mass_;
	}
	
	void ResidueModification2::setDiffMonoMass(double mass)
	{
		diff_mono_mass_ = mass;
	}

	double ResidueModification2::getDiffMonoMass() const
	{
		return diff_mono_mass_;
	}
	
	void ResidueModification2::setFormula(const String& composition)
	{
		formula_ = composition;
	}

	const String& ResidueModification2::getFormula() const
	{
		return formula_;
	}

	/*void ResidueModification2::setValidResidues(const vector<String>& valid_residues)
	{
		valid_residues_ = valid_residues;
	}

	const vector<String>& ResidueModification2::getValidResidues() const
	{
		return valid_residues_;
	}*/

	void ResidueModification2::addSynonym(const String& synonym)
	{
		synonyms_.insert(synonym);
	}

	void ResidueModification2::setSynonyms(const set<String>& synonyms)
	{
		synonyms_ = synonyms;
	}

	const set<String>& ResidueModification2::getSynonyms() const
	{
		return synonyms_;
	}
}

