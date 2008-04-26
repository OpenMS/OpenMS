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

#include <OpenMS/CHEMISTRY/ModificationDefinition.h>

using namespace std;

namespace OpenMS
{
	ModificationDefinition::ModificationDefinition()
		: allowed_position_(ResidueModification2::ANYWHERE),
			mod_(0),
			fixed_modification_(true),
			max_occurences_(0)
	{
	}

	ModificationDefinition::ModificationDefinition(const ModificationDefinition& rhs)
		: allowed_position_(rhs.allowed_position_),
			mod_(rhs.mod_),
			fixed_modification_(rhs.fixed_modification_),
			max_occurences_(rhs.max_occurences_)
	{
	}

	ModificationDefinition::ModificationDefinition(const String& mod)
	{
		
	}
	
	ModificationDefinition::~ModificationDefinition()
	{
	}

	void ModificationDefinition::setAllowedPosition(ResidueModification2::AllowedPosition pos)
	{
		allowed_position_ = pos;
	}

	ResidueModification2::AllowedPosition ModificationDefinition::getAllowedPosition() const
	{
		return allowed_position_;
	}
	
	void ModificationDefinition::setFixedModification(bool fixed_mod)
	{
		fixed_modification_ = fixed_mod;
	}

	bool ModificationDefinition::isFixedModification() const
	{
		return fixed_modification_;
	}

	void ModificationDefinition::setMaxOccurences(UInt max_occurences)
	{
		max_occurences_ = max_occurences;
	}

	UInt ModificationDefinition::getMaxOccurences() const
	{
		return max_occurences_;
	}
	
} // namespace OpenMS

