// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cstdlib>

namespace OpenMS
{

	String VersionInfo::getVersionAndTime()
	{
		return PACKAGE_VERSION " ("__DATE__", " __TIME__ ")";
	}

	String VersionInfo::getVersion()
	{
		return PACKAGE_VERSION;
	}

	Int VersionInfo::getMajorRevision()
	{
		static std::string release(PACKAGE_VERSION);
		static Int major = -1;
		if (major < 0)
		{
			size_t first_dot = release.find(".");
			std::string major_str(release, first_dot);
			major = atoi(major_str.c_str());
		}
		return major;
	}
	
	Int VersionInfo::getMinorRevision()
	{
		static std::string release(PACKAGE_VERSION);
		static Int minor = -1;
		if (minor < 0)
		{
			size_t first_dot = release.find(".");
			size_t second_dot = release.find(".", first_dot + 1);
			std::string minor_str(release, first_dot + 1, second_dot - first_dot - 1);
			minor = atoi(minor_str.c_str());
		}
		return minor;
	}
}
