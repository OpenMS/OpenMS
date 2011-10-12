// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors:  Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_VERSIONINFO_H
#define OPENMS_CONCEPT_VERSIONINFO_H

#include<OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
	class String;

	/**
		@brief Version information class.

		The OpenMS release version and revision data can be retrieved as a string
		or as integers.

		Note that the term <i>"version"</i> refers to releases (such as 1.0, 1.1, 1.1.1,
		1.2, ...),  whereas the term <i>"revision"</i> refers to a revision control system
		such as subversion and is mainly of interest for developers.

		The VersionInfo class contains only static methods.

		@ingroup Concept
	*/
	class OPENMS_DLLAPI VersionInfo
	{
		public:

    struct OPENMS_DLLAPI VersionDetails
    {
      Int version_major;
      Int version_minor;
      Int version_patch;

      VersionDetails()
        : version_major(0), version_minor(0), version_patch(0) 
      {
      }

      /** @brief parse String and return as proper struct 
      
          @returns VersionInfo::empty on failure

      */
      static VersionDetails create(const String& version);

      bool operator<(const VersionDetails& rhs) const;
      bool operator==(const VersionDetails& rhs) const;
      bool operator>(const VersionDetails& rhs) const;

      static const VersionDetails EMPTY; // 0.0.0 version for comparison
    };

		/// Return the build time of OpenMS
		static String getTime();

		/// Return the version number of OpenMS
		static String getVersion();

    /// Return the version number of OpenMS
    static VersionDetails getVersionStruct();

    /**
      @brief Return the revision number from revision control system, e.g. Subversion.

      On released versions of OpenMS (not from SVN), the result is "exported".
      The result can be possibly be "" on some platforms, which means that
      revision info is unavailable.  You should check for both cases in your
      code.

      @internal The current svn version is queried by the build system regularly and
      the result is written as a header file which is
      included by VersionInfo.C.
		*/
		static String getRevision();

	};

}

#endif // OPENMS_CONCEPT_VERSIONINFO_H
