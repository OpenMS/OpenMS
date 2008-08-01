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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <iostream>
///////////////////////////
#include <OpenMS/FORMAT/FastaIterator.h>
#include <OpenMS/CHEMISTRY/PepIterator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FastaIterator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
typedef std::pair <String, String> FASTAEntry;

FastaIterator* ptr = 0;
CHECK(FastaIterator())
	ptr = new FastaIterator();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~FastaIterator())
	delete ptr;
RESULT

CHECK(virtual void setFastaFile(const String &f))
	ptr = new FastaIterator();
	TEST_EXCEPTION (Exception::FileNotFound,ptr->setFastaFile("FileThatNotExists"));
	TEST_EXCEPTION (Exception::FileNotFound,ptr->setFastaFile(""));
	ptr->setFastaFile("data/FastaIterator_test.fasta");
RESULT

CHECK(String getFastaFile ())
	ptr = new FastaIterator();
	TEST_EQUAL (ptr->getFastaFile(),"");
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	TEST_EQUAL (ptr->getFastaFile(),"data/FastaIterator_test.fasta");
RESULT


CHECK(virtual FASTAEntry operator *())
	ptr = new FastaIterator();
	TEST_EXCEPTION (Exception::InvalidIterator,**ptr);
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	ptr->begin();
	FASTAEntry fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 1");
	TEST_EQUAL(fe.second,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
RESULT

CHECK(virtual PepIterator& operator++())
	ptr = new FastaIterator();
	TEST_EXCEPTION (Exception::InvalidIterator, ++(*ptr));
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	ptr->begin();
	PepIterator & pepIt = ++(*ptr);
	TEST_EQUAL ((*pepIt).first,(**ptr).first);
	TEST_EQUAL ((*pepIt).second,(**ptr).second);
	pepIt = ++(*ptr);
	TEST_EQUAL ((*pepIt).first,(**ptr).first);
	TEST_EQUAL ((*pepIt).second,(**ptr).second);
RESULT

CHECK(virtual PepIterator* operator++(int i))
	ptr = new FastaIterator();
	TEST_EXCEPTION (Exception::InvalidIterator, (*ptr)++);
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	ptr->begin();
	FASTAEntry fe = **ptr;
	PepIterator * pepIt = (*ptr)++;
	TEST_EQUAL ((**pepIt).first,fe.first);
	TEST_EQUAL ((**pepIt).second,fe.second);
RESULT

CHECK(virtual bool begin())
	ptr = new FastaIterator();
	TEST_EXCEPTION (Exception::InvalidIterator, (*ptr).begin());
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	ptr->begin();
	FASTAEntry fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 1");
	TEST_EQUAL(fe.second,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
RESULT

CHECK(virtual bool isAtEnd())
	ptr = new FastaIterator();
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	ptr->begin();
	for (int i = 0; i < 5;i++)
	{
		TEST_EQUAL(ptr->isAtEnd(),0);
		++(*ptr);
	}
	TEST_EQUAL(ptr->isAtEnd(),1);
RESULT

CHECK(FastaIterator(const FastaIterator &))
	ptr = new FastaIterator();
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	ptr->begin();
	++*ptr;
	++*ptr;
	FastaIterator copy (*ptr);
	TEST_EQUAL((**ptr).first,(*copy).first);
	TEST_EQUAL((**ptr).second,(*copy).second);
	TEST_EQUAL((*ptr).getFastaFile(),(copy).getFastaFile());
RESULT



CHECK(virtual void setSpectrum(const std::vector< float > &))
	const std::vector<float> spec;
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).setSpectrum(spec));
RESULT

CHECK(virtual const std::vector<float>& getSpectrum())
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).getSpectrum());
RESULT

CHECK(virtual void setTolerance(float))
	float t = 0.5;
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).setTolerance(t));
RESULT

CHECK(virtual float getTolerance())
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).getTolerance());
RESULT

CHECK(static const String getProductName())
	ptr = new FastaIterator();
	TEST_EQUAL(ptr->getProductName(),"FastaIterator");
RESULT

CHECK(static PepIterator* create())
	ptr = new FastaIterator();
	TEST_NOT_EQUAL(ptr->create(),0);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
