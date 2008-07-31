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

#ifndef OPENMS_CONCEPT_VERSIONINFO_H
#define OPENMS_CONCEPT_VERSIONINFO_H

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
	class String;
	
	/**	
		@brief Version information class.
		
		The OpenMS release and version data can be retrieved as a string
		or as integers.
			
		The VersionInfo class contains only static methods.

		@ingroup Concept
	*/
	class VersionInfo
	{
		public:

		///Return the version number and the build time of OpenMS
		static String getVersionAndTime();

		///Return the version number of OpenMS
		static String getVersion();

		/// Return the major revision number. The part of the release number before the dot.
		static Int getMajorRevision();

		///Return the minor revision number. The part of the release number after the dot.
		static Int getMinorRevision();
	};
	
}

#endif // OPENMS_CONCEPT_VERSIONINFO_H
