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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace OpenMS;
using namespace std;

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

int main( int /*argc*/, const char** /*argv*/ )
{
	cout << "OpenMS Version:" << endl;
	cout << "==================" << endl;
	cout << "Version      : " << VersionInfo::getVersion() << endl;
	cout << "Build time   : " << VersionInfo::getTime() << endl;
	cout << "SVN revision : " << VersionInfo::getRevision() << endl;
	cout << endl;
	cout << "Installation information:" << endl;
	cout << "==================" << endl;
	cout << "Data path    : " << OPENMS_DATA_PATH << endl;
	cout << endl;
	cout << "Build information:" << endl;
	cout << "==================" << endl;
	cout << "Source path  : " << OPENMS_SOURCE_PATH << endl;
	cout << "Binary path  : " << OPENMS_BINARY_PATH << endl;
	cout << endl;
	
	return 0;
}

/// @endcond



