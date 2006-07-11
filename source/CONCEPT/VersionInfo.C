// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

namespace OpenMS
{

	const char* VersionInfo::getVersion() throw()
	{
		return OPENMS_RELEASE_STRING " ("__DATE__", " __TIME__ ")";
	}

	int VersionInfo::getMinorRevision()
	{
		static std::string release(OPENMS_RELEASE_STRING);
		static int minor = -1;
		if (minor < 0)
		{
			size_t first_dot = release.find(".");
			std::string minor_str(release, first_dot);
			minor = atoi(minor_str.c_str());
		}
		return minor;
	}
	
	int VersionInfo::getMajorRevision()
	{
		static std::string release(OPENMS_RELEASE_STRING);
		static int major = -1;
		if (major < 0)
		{
			size_t first_dot = release.find(".");
			size_t second_dot = release.find(".", first_dot + 1);
			std::string major_str(release, first_dot + 1, second_dot - first_dot - 1);
			major = atoi(major_str.c_str());
		}
		return major;
	}
	
	const char* VersionInfo::RELEASE_DATE_AND_TIME = "@(#)OpenMS kernel compiled on: " __DATE__ ", " __TIME__;

	const char* VersionInfo::RELEASE_WHAT_STRING = "@(#)OpenMS kernel version: V." OPENMS_RELEASE_STRING;					  

}
