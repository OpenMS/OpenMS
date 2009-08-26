// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm, Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cstdlib>
#include <fstream>

using namespace std;
//#ifndef  OPENMS_REVISION
// #warning is not a standard preprocessor directive, but supported by many compilers, including GCC.
//#ifdef   OPENMS_COMPILER_GXX
//#warning "Note: OPENMS_REVISION is undefined.  OpenMS uses the subversion software for revision control, but now revision info is unavailable.  This is normally the case when you are compiling a released version (without .svn subdirectories)."
//#endif
//#define OPENMS_REVISION ""
//#endif

namespace OpenMS
{
	

	String VersionInfo::getTime()
	{
		static bool is_initialized = false;
		static String result;
		if ( !is_initialized )
		{
			result = String(__DATE__) + ", " + __TIME__;
			is_initialized = true;
		}
		return result;
	}

	String VersionInfo::getVersion()
	{
		static bool is_initialized = false;
		static String result;
		if ( !is_initialized )
		{
			result = OPENMS_PACKAGE_VERSION;
			result.trim();
			is_initialized = true;
		}
		return result;
	}

	Int VersionInfo::getMajorVersion()
	{
		static Int major = -1;
		if (major < 0)
		{
			const String version = getVersion();
			size_t first_dot = version.find('.');
			major = atoi(version.substr(0,first_dot).c_str());
		}
		return major;
	}
	
	Int VersionInfo::getMinorVersion()
	{
		static Int minor = -1;
		if (minor < 0)
		{
			const String version = getVersion();
			size_t first_dot = version.find('.');
			size_t second_dot = version.find('.', first_dot + 1);
			if ( second_dot > 1000 ) second_dot = 1000; // avoid arithmetic overflow if there is no second dot
			minor = atoi(version.substr(first_dot+1,second_dot-first_dot+1).c_str());
		}
		return minor;
	}

	String VersionInfo::getRevision()
	{
		static bool is_initialized = false;
		static String result;
		if ( !is_initialized )
		{
#ifdef OPENMS_HAS_SVNVERSION			
			ifstream in(PACKAGE_REVISION_FILE);
			getline(in, result, '\n');
			result.trim();
#else
			result = "<unknown>";
#endif
			is_initialized = true;
		}
		return result;
	}

}
