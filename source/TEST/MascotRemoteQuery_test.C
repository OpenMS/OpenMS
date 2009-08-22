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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/MascotRemoteQuery.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MascotRemoteQuery, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MascotRemoteQuery* ptr = 0;
START_SECTION(MascotRemoteQuery())
{
	ptr = new MascotRemoteQuery();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(virtual ~MascotRemoteQuery())
{
	delete ptr;
}
END_SECTION

START_SECTION((MascotRemoteQuery(QObject *parent=0)))
{
  // TODO
}
END_SECTION

START_SECTION((MascotRemoteQuery(const MascotRemoteQuery &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void setQuerySpectra(const String &exp)))
{
  // TODO
}
END_SECTION

START_SECTION((const QByteArray& getMascotXMLResponse() const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool hasError() const ))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getErrorMessage() const ))
{
  // TODO
}
END_SECTION

START_SECTION((MascotRemoteQuery& operator=(const MascotRemoteQuery &rhs)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



