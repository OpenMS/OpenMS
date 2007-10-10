// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bauer$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <iostream>
///////////////////////////
#include <OpenMS/FORMAT/FastaIteratorIntern.h>
#include <OpenMS/CHEMISTRY/PepIterator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FastaIteratorIntern, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
typedef std::pair <String, String> FASTAEntry;

FastaIteratorIntern* ptr = 0;

CHECK(FastaIteratorIntern())
	ptr = new FastaIteratorIntern();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~FastaIteratorIntern())
	delete ptr;
RESULT

CHECK(virtual void setFastaFile(const String &f) throw (Exception::FileNotFound, Exception::ParseError))
	ptr = new FastaIteratorIntern();
	TEST_EXCEPTION (Exception::FileNotFound,ptr->setFastaFile("FileThatNotExists"));
	TEST_EXCEPTION (Exception::FileNotFound,ptr->setFastaFile(""));
	ptr->setFastaFile("data/FastaIterator_test.fasta");
RESULT

CHECK(String getFastaFile ())
	ptr = new FastaIteratorIntern();
	TEST_EQUAL (ptr->getFastaFile(),"");
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	TEST_EQUAL (ptr->getFastaFile(),"data/FastaIterator_test.fasta");
RESULT


CHECK(virtual FASTAEntry operator *() throw (Exception::InvalidIterator))
	ptr = new FastaIteratorIntern();
	TEST_EXCEPTION (Exception::InvalidIterator,**ptr);
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	ptr->begin();
	FASTAEntry fe = **ptr;
	TEST_EQUAL(fe.first,"Entry 1");
	TEST_EQUAL(fe.second,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
RESULT

CHECK(virtual PepIterator& operator++() throw (Exception::InvalidIterator))
	ptr = new FastaIteratorIntern();
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

CHECK(virtual PepIterator* operator++(int i) throw (Exception::InvalidIterator))
	ptr = new FastaIteratorIntern();
	TEST_EXCEPTION (Exception::InvalidIterator, (*ptr)++);
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	ptr->begin();
	FASTAEntry fe = **ptr;
	PepIterator * pepIt = (*ptr)++;
	TEST_EQUAL ((**pepIt).first,fe.first);
	TEST_EQUAL ((**pepIt).second,fe.second);
RESULT

CHECK(virtual bool begin() throw (Exception::InvalidIterator))
	ptr = new FastaIteratorIntern();
	TEST_EXCEPTION (Exception::InvalidIterator, (*ptr).begin());
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	ptr->begin();
	FASTAEntry fe = **ptr;
	TEST_EQUAL(fe.first,"Entry 1");
	TEST_EQUAL(fe.second,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
RESULT

CHECK(virtual bool isAtEnd())
	ptr = new FastaIteratorIntern();
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	ptr->begin();
	for (int i = 0; i < 5;i++)
	{
		TEST_EQUAL(ptr->isAtEnd(),0);
		++(*ptr);
	}
	TEST_EQUAL(ptr->isAtEnd(),1);
RESULT

CHECK(FastaIteratorIntern(const FastaIteratorIntern &))
	ptr = new FastaIteratorIntern();
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	ptr->begin();
	++*ptr;
	++*ptr;
	FastaIteratorIntern copy (*ptr);
	TEST_EQUAL((**ptr).first,(*copy).first);
	TEST_EQUAL((**ptr).second,(*copy).second);
	TEST_EQUAL((*ptr).getFastaFile(),(copy).getFastaFile());
RESULT


CHECK(virtual void setSpectrum(const std::vector< float > &) throw (Exception::NotImplemented))
	const std::vector<float> spec;
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).setSpectrum(spec));
RESULT

CHECK(virtual const std::vector<float>& getSpectrum() throw (Exception::NotImplemented))
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).getSpectrum());
RESULT

CHECK(virtual void setTolerance(float) throw (Exception::NotImplemented))
	float t = 0.5;
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).setTolerance(t));
RESULT

CHECK(virtual float getTolerance() throw (Exception::NotImplemented))
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).getTolerance());
RESULT

CHECK(static const std::string getName())
	ptr = new FastaIteratorIntern();
	TEST_EQUAL(ptr->getName(),"FastaIteratorIntern");
RESULT

CHECK(static PepIterator* create())
	ptr = new FastaIteratorIntern();
	TEST_NOT_EQUAL(ptr->create(),0);
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
