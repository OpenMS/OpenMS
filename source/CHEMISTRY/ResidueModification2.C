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
			composition_(rhs.composition_),
			valid_residues_(rhs.valid_residues_)
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
		valid_residues_ = rhs.valid_residues_;
		
		return *this;
  }
	
	bool ResidueModification2::operator == (const ResidueModification2& rhs) const
	{
		return  title_ == rhs.title_ &&
						full_name_ == rhs.full_name_ &&
						allowed_position_ == rhs.allowed_position_ &&
						site_ == rhs.site_ &&
						classification_ == rhs.classification_ &&
						average_mass_ == rhs.average_mass_ &&
						mono_mass_ == rhs.mono_mass_ && 
						composition_ == rhs.composition_ &&
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

	void ResidueModification2::setValidResidues(const vector<String>& valid_residues)
	{
		valid_residues_ = valid_residues;
	}

	const vector<String>& ResidueModification2::getValidResidues() const
	{
		return valid_residues_;
	}
}

