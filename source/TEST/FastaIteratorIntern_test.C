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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
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

START_SECTION(FastaIteratorIntern())
	ptr = new FastaIteratorIntern();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~FastaIteratorIntern())
	delete ptr;
END_SECTION

START_SECTION(virtual void setFastaFile(const String &f))
	ptr = new FastaIteratorIntern();
	TEST_EXCEPTION (Exception::FileNotFound,ptr->setFastaFile("FileThatNotExists"));
	TEST_EXCEPTION (Exception::FileNotFound,ptr->setFastaFile(""));
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("FastaIterator_test.fasta"));
END_SECTION

START_SECTION(String getFastaFile ())
	ptr = new FastaIteratorIntern();
	TEST_EQUAL (ptr->getFastaFile(),"");
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("FastaIterator_test.fasta"));
	TEST_EQUAL (ptr->getFastaFile(),OPENMS_GET_TEST_DATA_PATH("FastaIterator_test.fasta"));
END_SECTION


START_SECTION(virtual FASTAEntry operator *())
	ptr = new FastaIteratorIntern();
	TEST_EXCEPTION (Exception::InvalidIterator,**ptr);
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("FastaIterator_test.fasta"));
	ptr->begin();
	FASTAEntry fe = **ptr;
	TEST_EQUAL(fe.first,"Entry 1");
	TEST_EQUAL(fe.second,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
END_SECTION

START_SECTION(virtual PepIterator& operator++())
	ptr = new FastaIteratorIntern();
	TEST_EXCEPTION (Exception::InvalidIterator, ++(*ptr));
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("FastaIterator_test.fasta"));
	ptr->begin();
	PepIterator & pepIt = ++(*ptr);
	TEST_EQUAL ((*pepIt).first,(**ptr).first);
	TEST_EQUAL ((*pepIt).second,(**ptr).second);
	pepIt = ++(*ptr);
	TEST_EQUAL ((*pepIt).first,(**ptr).first);
	TEST_EQUAL ((*pepIt).second,(**ptr).second);
END_SECTION

START_SECTION(virtual PepIterator* operator++(int i))
	ptr = new FastaIteratorIntern();
	TEST_EXCEPTION (Exception::InvalidIterator, (*ptr)++);
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("FastaIterator_test.fasta"));
	ptr->begin();
	FASTAEntry fe = **ptr;
	PepIterator * pepIt = (*ptr)++;
	TEST_EQUAL ((**pepIt).first,fe.first);
	TEST_EQUAL ((**pepIt).second,fe.second);
END_SECTION

START_SECTION(virtual bool begin())
	ptr = new FastaIteratorIntern();
	TEST_EXCEPTION (Exception::InvalidIterator, (*ptr).begin());
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("FastaIterator_test.fasta"));
	ptr->begin();
	FASTAEntry fe = **ptr;
	TEST_EQUAL(fe.first,"Entry 1");
	TEST_EQUAL(fe.second,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
END_SECTION

START_SECTION(virtual bool isAtEnd())
	ptr = new FastaIteratorIntern();
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("FastaIterator_test.fasta"));
	ptr->begin();
	for (int i = 0; i < 5;i++)
	{
		TEST_EQUAL(ptr->isAtEnd(), false);
		++(*ptr);
	}
	TEST_EQUAL(ptr->isAtEnd(), true);
END_SECTION

START_SECTION(FastaIteratorIntern(const FastaIteratorIntern &))
	ptr = new FastaIteratorIntern();
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("FastaIterator_test.fasta"));
	ptr->begin();
	++*ptr;
	++*ptr;
	FastaIteratorIntern copy (*ptr);
	TEST_EQUAL((**ptr).first,(*copy).first);
	TEST_EQUAL((**ptr).second,(*copy).second);
	TEST_EQUAL((*ptr).getFastaFile(),(copy).getFastaFile());
END_SECTION


START_SECTION(virtual void setSpectrum(const std::vector< DoubleReal > &))
	const std::vector<DoubleReal> spec;
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).setSpectrum(spec));
END_SECTION

START_SECTION(virtual const std::vector<DoubleReal>& getSpectrum())
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).getSpectrum());
END_SECTION

START_SECTION(virtual void setTolerance(DoubleReal))
	DoubleReal t = 0.5;
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).setTolerance(t));
END_SECTION

START_SECTION(virtual DoubleReal getTolerance())
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).getTolerance());
END_SECTION

START_SECTION(static const std::string getProductName())
	ptr = new FastaIteratorIntern();
	TEST_EQUAL(ptr->getProductName(),"FastaIteratorIntern");
END_SECTION

START_SECTION(static PepIterator* create())
	ptr = new FastaIteratorIntern();
	TEST_NOT_EQUAL(ptr->create(),0);
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
