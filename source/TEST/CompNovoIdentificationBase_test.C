// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationBase.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIdentificationBase, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(CompNovoIdentificationBase())
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((CompNovoIdentificationBase(const CompNovoIdentificationBase &source)))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual ~CompNovoIdentificationBase()))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void getIdentifications(std::vector< PeptideIdentification > &ids, const PeakMap &exp)=0))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((CompNovoIdentificationBase& operator=(const CompNovoIdentificationBase &source)))
{
	NOT_TESTABLE
}
END_SECTION


std::set<String>str_set;
str_set.insert("TESTSTRING");
std::set<String>::const_iterator it = str_set.begin();


START_SECTION([CompNovoIdentificationBase::Permut] Permut(const std::set< String >::const_iterator &permut, DoubleReal s))
	CompNovoIdentificationBase::Permut perm(it, 50.0);
	TEST_EQUAL(perm.getScore(), 50.0)
	TEST_EQUAL(*perm.getPermut(), "TESTSTRING")
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] Permut(const Permut &rhs))
	CompNovoIdentificationBase::Permut perm(it, 50.0);
	CompNovoIdentificationBase::Permut copy(perm);
	TEST_EQUAL(perm.getScore(), copy.getScore())
	TEST_EQUAL(*perm.getPermut(), *copy.getPermut())
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] Permut& operator=(const Permut &rhs))
	CompNovoIdentificationBase::Permut perm(it, 50.0);
	CompNovoIdentificationBase::Permut copy(it, 0.0);
	copy=perm;
	TEST_EQUAL(perm.getScore(), copy.getScore())
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] virtual ~Permut())
	CompNovoIdentificationBase::Permut * ptr = new CompNovoIdentificationBase::Permut(it, 50.0);
	delete ptr;
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] void setPermut(const std::set< String >::const_iterator &it))
  std::set<String>str_set;
  str_set.insert("zero");
  std::set<String>::const_iterator it_zero = str_set.begin();
	CompNovoIdentificationBase::Permut perm(it_zero, 50.0);
	perm.setPermut(it);
	TEST_EQUAL(*perm.getPermut(), "TESTSTRING");
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] void setScore(DoubleReal score))
	CompNovoIdentificationBase::Permut perm(it, 50.0);
	perm.setScore(0.0);
	TEST_EQUAL(perm.getScore(), 0.0)
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] DoubleReal getScore() const)
	NOT_TESTABLE //already tested above
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] const std::set<String>::const_iterator& getPermut() const)
	NOT_TESTABLE //already tested above
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



