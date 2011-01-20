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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/DeNovoAlgorithm.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DeNovoAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(DeNovoAlgorithm())
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual ~DeNovoAlgorithm()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DeNovoAlgorithm(const DeNovoAlgorithm &rhs)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DeNovoAlgorithm& operator=(const DeNovoAlgorithm &rhs)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void generateCandidates(std::vector< PeptideIdentification > &candidates, const std::vector< std::vector< DeNovoIonScoring::IonScore > > &ion_scores, const RichPeakMap &exp)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void generateCandidates(PeptideIdentification &candidates, std::vector< DeNovoIonScoring::IonScore > &ion_scores, const RichPeakSpectrum &spec)=0))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



