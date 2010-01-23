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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/DeNovoIdentification.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DeNovoIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(DeNovoIdentification())
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual ~DeNovoIdentification()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DeNovoIdentification(const DeNovoIdentification &rhs)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DeNovoIdentification& operator=(const DeNovoIdentification &rhs)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void getIdentifications(std::vector< PeptideIdentification > &ids, const RichPeakMap &exp)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void getIdentification(PeptideIdentification &id, const RichPeakSpectrum &spectrum)=0))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



