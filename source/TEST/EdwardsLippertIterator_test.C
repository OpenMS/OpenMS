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
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/EdwardsLippertIterator.h>
#include <OpenMS/CHEMISTRY/PepIterator.h>
#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(EdwardsLippertIterator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
typedef std::pair <String, String> FASTAEntry;

vector<float> spec;
spec.push_back(178.1864);
spec.push_back(441.4806);
const vector<float> specc (spec);

EdwardsLippertIterator* ptr = 0;
CHECK(EdwardsLippertIterator())
        ptr = new EdwardsLippertIterator();
	TEST_REAL_EQUAL(0.5,ptr->getTolerance());
        TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~EdwardsLippertIterator())
        delete ptr;
RESULT

CHECK(EdwardsLippertIterator(const EdwardsLippertIterator &))
	ptr = new EdwardsLippertIterator();
	ptr->setFastaFile("data/EdwardsLippertIterator_test.fasta");
	ptr->setSpectrum(specc);
	ptr->begin();
	++*ptr;
	EdwardsLippertIterator copy (*ptr);
	TEST_EQUAL ((*ptr).getFastaFile(),(copy).getFastaFile());
	TEST_EQUAL ((*ptr).getTolerance(),(copy).getTolerance());
	TEST_EQUAL ((**ptr).first,(*copy).first);
	TEST_EQUAL ((**ptr).second,(*copy).second);
RESULT

CHECK(virtual void setFastaFile(const String &f))
	ptr = new EdwardsLippertIterator();
	TEST_EXCEPTION (Exception::FileNotFound,ptr->setFastaFile("FileThatNotExists"));
	TEST_EXCEPTION (Exception::FileNotFound,ptr->setFastaFile(""));
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	TEST_EQUAL (ptr->getFastaFile(),"data/FastaIterator_test.fasta");
RESULT

CHECK(String getFastaFile ())
	ptr = new EdwardsLippertIterator();
	TEST_EQUAL (ptr->getFastaFile(),"");
	ptr->setFastaFile("data/FastaIterator_test.fasta");
	TEST_EQUAL (ptr->getFastaFile(),"data/FastaIterator_test.fasta");
RESULT



CHECK(virtual FASTAEntry operator *())
	float masse[255];
	ResidueDB* rdb = ResidueDB::getInstance();
		
	char aa[] = "ARNDCEQGHILKMFPSTWYV";
		
	for (unsigned int i = 0; i<255;++i)
	{
		masse[i]=0;
	}
	for (unsigned int i = 0; i<strlen(aa);++i)
	{
		const Residue * r = rdb->getResidue(aa[i]);
		masse[(int)aa[i]]=r->getAverageWeight();
	}
	ptr = new EdwardsLippertIterator();
	TEST_EXCEPTION (Exception::InvalidIterator,**ptr);
	ptr->setFastaFile("data/EdwardsLippertIterator_test.fasta");
	ptr->setSpectrum(specc);
	ptr->begin();
	for (unsigned int i = 0; i < 29;++i)
	{
		FASTAEntry fe = **ptr;
		TEST_EQUAL(fe.first,">Entry 1");
		TEST_EQUAL(fe.second,"AA");
		++(*ptr);
	}
	FASTAEntry fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 2");
	TEST_EQUAL(fe.second,"AA");
	
	++*ptr;
	fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 4");
	TEST_EQUAL(fe.second,"EEE");


	ptr = new EdwardsLippertIterator();
	ptr->setFastaFile("data/EdwardsLippertIterator_test_2.fasta");
	ptr->setSpectrum(specc);
	ptr->begin();
	float tol = 0.2;
	ptr->setTolerance(tol);
	while (!ptr->isAtEnd())
	{
		String seq = (**ptr).second;
		(++*ptr);
		float m = 0;
		for (unsigned int i = 0; i < seq.length();i++)
		{
			m+=masse[(int)seq[i]];
		}
		bool is_in_spec = false;
		for (unsigned int i = 0; i < specc.size();i++)
		{
			is_in_spec |= (m>=specc.at(i)-tol&&m<=specc.at(i)+tol);
		}
		TEST_EQUAL(is_in_spec,1);
	}
RESULT

CHECK(virtual PepIterator& operator++())
	ptr = new EdwardsLippertIterator();
	TEST_EXCEPTION (Exception::InvalidIterator, ++(*ptr));
	ptr->setFastaFile("data/EdwardsLippertIterator_test.fasta");
	ptr->setSpectrum(specc);
	ptr->begin();
	PepIterator & pepIt = ++(*ptr);
	TEST_EQUAL ((*pepIt).first,(**ptr).first);
	TEST_EQUAL ((*pepIt).second,(**ptr).second);
	pepIt = ++(*ptr);
	TEST_EQUAL ((*pepIt).first,(**ptr).first);
	TEST_EQUAL ((*pepIt).second,(**ptr).second);
RESULT

CHECK(virtual PepIterator* operator++(int i))
	ptr = new EdwardsLippertIterator();
	TEST_EXCEPTION (Exception::InvalidIterator, (*ptr)++);
	ptr->setFastaFile("data/EdwardsLippertIterator_test.fasta");
	ptr->setSpectrum(specc);
	ptr->begin();
	FASTAEntry fe = **ptr;
	PepIterator * pepIt = (*ptr)++;
	TEST_EQUAL ((**pepIt).first,fe.first);
	TEST_EQUAL ((**pepIt).second,fe.second);
RESULT

CHECK(virtual bool begin())
	ptr = new EdwardsLippertIterator();
	TEST_EXCEPTION (Exception::InvalidIterator, (*ptr).begin());
	ptr->setFastaFile("data/EdwardsLippertIterator_test.fasta");
	ptr->setSpectrum(specc);
	ptr->begin();
	FASTAEntry fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 1");
	TEST_EQUAL(fe.second,"AA");
RESULT

CHECK(bool isAtEnd ())
	ptr = new EdwardsLippertIterator();
	ptr->setFastaFile("data/EdwardsLippertIterator_test.fasta");
	ptr->setSpectrum(specc);
	ptr->begin();
	for (int i = 0; i < 58;i++)
	{
		TEST_EQUAL(ptr->isAtEnd(),0);
		(*ptr)++;
	}
	TEST_EQUAL(ptr->isAtEnd(),1);
RESULT

CHECK(virtual void setTolerance(float t))
	ptr = new EdwardsLippertIterator();
	ptr->setTolerance(0.4);
	TEST_REAL_EQUAL(0.4,ptr->getTolerance());
	TEST_EXCEPTION (Exception::InvalidValue,ptr->setTolerance(-0.1));
RESULT

CHECK(virtual float getTolerance())
	ptr = new EdwardsLippertIterator();
	TEST_REAL_EQUAL(0.5,ptr->getTolerance());
	ptr->setTolerance(0.4);
	TEST_REAL_EQUAL(0.4,ptr->getTolerance());
RESULT

CHECK(virtual void setSpectrum(const std::vector< float > &s))
	ptr = new EdwardsLippertIterator();
	ptr->setSpectrum(specc);
	vector<float> spec2;
	spec2.push_back(441.4806);
	spec2.push_back(178.1864);
	const vector<float> specc2 (spec2);
	TEST_EXCEPTION (Exception::InvalidValue,ptr->setSpectrum(specc2));
RESULT

CHECK(virtual const std::vector<float>& getSpectrum())
	ptr = new EdwardsLippertIterator();
	ptr->setSpectrum(specc);
	TEST_EQUAL(specc.size(),ptr->getSpectrum().size());
	for (unsigned int i = 0;i < specc.size();++i)
	{
		TEST_EQUAL(specc.at(i),ptr->getSpectrum().at(i));
	}
RESULT

CHECK(virtual bool isDigestingEnd(char, char))
	ptr = new EdwardsLippertIterator();
	TEST_EQUAL(ptr->isDigestingEnd('R','S'),1)
	TEST_EQUAL(ptr->isDigestingEnd('K','S'),1)
	TEST_EQUAL(ptr->isDigestingEnd('R','P'),1)
	TEST_EQUAL(ptr->isDigestingEnd('K','P'),1)
	TEST_EQUAL(ptr->isDigestingEnd('S','S'),1)
RESULT

CHECK(static const String getProductName())
	ptr = new EdwardsLippertIterator();
	TEST_EQUAL(ptr->getProductName(),"EdwardsLippertIterator");
RESULT

CHECK(static PepIterator* create())
	ptr = new EdwardsLippertIterator();
	TEST_NOT_EQUAL(ptr->create(),0);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
