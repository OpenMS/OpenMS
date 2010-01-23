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

#include <OpenMS/DATASTRUCTURES/CVReference.h>

using namespace std;

namespace OpenMS
{
	// CV reference implementation
	CVReference::CVReference()
	{
	}

	CVReference::~CVReference()
	{
	}
	
	CVReference::CVReference(const CVReference& rhs)
		: name_(rhs.name_),
			identifier_(rhs.identifier_)
	{
	}

	CVReference& CVReference::operator = (const CVReference& rhs)
	{
		if (this != &rhs)
		{
			name_ = rhs.name_;
			identifier_ = rhs.identifier_;
		}
		return *this;
	}

	bool CVReference::operator == (const CVReference& rhs) const
	{
		return name_ == rhs.name_ && identifier_ == rhs.identifier_;
	}

	bool CVReference::operator != (const CVReference& rhs) const
	{
		return !(*this == rhs);
	}

	void CVReference::setName(const String& name)
	{
		name_ = name;
	}

	const String& CVReference::getName() const
	{
		return name_;
	}

	void CVReference::setIdentifier(const String& identifier)
	{
		identifier_ = identifier;
	}

	const String& CVReference::getIdentifier() const
	{
		return identifier_;
	}

} // namespace OpenMS



