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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringBase.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIonScoringBase, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(CompNovoIonScoringBase())
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((CompNovoIonScoringBase(const CompNovoIonScoringBase &source)))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual ~CompNovoIonScoringBase()))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((DoubleReal scoreIsotopes(const PeakSpectrum &CID_spec, PeakSpectrum::ConstIterator it, Size charge)))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((CompNovoIonScoringBase& operator=(const CompNovoIonScoringBase &source)))
{
	NOT_TESTABLE
}
END_SECTION

CompNovoIonScoringBase::IonScore * ptr;
START_SECTION([CompNovoIonScoringBase::IonScore] IonScore())
	ptr=new CompNovoIonScoringBase::IonScore();
	TEST_NOT_EQUAL(ptr,0)
END_SECTION

START_SECTION([CompNovoIonScoringBase::IonScore] IonScore(const IonScore &rhs))
	CompNovoIonScoringBase::IonScore ion_score;
	ion_score.s_bion=5.0;
	TEST_EQUAL(CompNovoIonScoringBase::IonScore(ion_score).s_bion, 5.0)
END_SECTION

START_SECTION([CompNovoIonScoringBase::IonScore] IonScore& operator=(const IonScore &rhs))
	CompNovoIonScoringBase::IonScore ion_score, copy;
	ion_score.s_bion=5.0;
	copy = ion_score;
	TEST_EQUAL(copy.s_bion, ion_score.s_bion)
END_SECTION

START_SECTION([CompNovoIonScoringBase::IonScore] virtual ~IonScore())
	delete ptr;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



