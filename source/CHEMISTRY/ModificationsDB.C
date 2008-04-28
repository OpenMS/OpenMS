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
#include <vector>

using namespace std;

namespace OpenMS
{
	ModificationsDB::ModificationsDB()
	{
		// read the modifications from unimod.xml
		UnimodXMLFile().load("CHEMISTRY/unimod.xml", mods_);

		for (vector<ResidueModification2>::const_iterator it = mods_.begin(); it !=mods_.end(); ++it)
		{
			modification_names_[it->getFullName()] = &*it;
		}
	}

	ModificationsDB::ModificationsDB(const ModificationsDB& /*res_db*/)
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}
	

	ModificationsDB::~ModificationsDB()
	{
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
		return mods_[index];
	}
} // namespace OpenMS

