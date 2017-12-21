// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/TrypticIterator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TrypticIterator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
typedef std::pair <String, String> FASTAEntry;

TrypticIterator* ptr = nullptr;
TrypticIterator* nullPointer = nullptr;
START_SECTION(TrypticIterator())
	ptr = new TrypticIterator();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~TrypticIterator())
	delete ptr;
END_SECTION

START_SECTION(TrypticIterator(const TrypticIterator &))
	ptr = new TrypticIterator();
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("TrypticIterator_test.fasta"));
	ptr->begin();
	++*ptr;
	TrypticIterator copy (*ptr);
	TEST_EQUAL ((*ptr).getFastaFile(),(copy).getFastaFile());
	TEST_EQUAL ((**ptr).first,(*copy).first);
	TEST_EQUAL ((**ptr).second,(*copy).second);
END_SECTION

START_SECTION(virtual void setFastaFile(const String &f))
	ptr = new TrypticIterator();
	TEST_EXCEPTION (Exception::FileNotFound,ptr->setFastaFile("FileThatNotExists"));
	TEST_EXCEPTION (Exception::FileNotFound,ptr->setFastaFile(""));
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("TrypticIterator_test.fasta"));
END_SECTION

START_SECTION(virtual String getFastaFile())
	ptr = new TrypticIterator();
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("TrypticIterator_test.fasta"));
	TEST_EQUAL(ptr->getFastaFile(),OPENMS_GET_TEST_DATA_PATH("TrypticIterator_test.fasta"));
END_SECTION

START_SECTION(static const String getProductName())
	ptr = new TrypticIterator();
	TEST_EQUAL(ptr->getProductName(),"TrypticIterator");
END_SECTION

START_SECTION(static PepIterator* create())
	ptr = new TrypticIterator();
  TEST_NOT_EQUAL(ptr->create(),nullPointer);
END_SECTION

START_SECTION(virtual FASTAEntry operator *())
	ptr = new TrypticIterator();
	TEST_EXCEPTION (Exception::InvalidIterator,**ptr);
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("TrypticIterator_test.fasta"));
	ptr->begin();
	FASTAEntry fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 1");
	TEST_EQUAL(fe.second,"AAAAAK");
	++*ptr;
	fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 1");
	TEST_EQUAL(fe.second,"AAAAAKAAAAAAAAAAAAAAAAAAAAAAAA");
	++*ptr;
	fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 1");
	TEST_EQUAL(fe.second,"AAAAAAAAAAAAAAAAAAAAAAAA");
	++*ptr;
	fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 2");
	TEST_EQUAL(fe.second,"K");
	++*ptr;
	fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 2");
	TEST_EQUAL(fe.second,"KCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
	++*ptr;
	fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 2");
	TEST_EQUAL(fe.second,"CCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
	++*ptr;
	fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 3");
	TEST_EQUAL(fe.second,"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDK");
	++*ptr;
	fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 4");
	TEST_EQUAL(fe.second,"EEEEEK");
END_SECTION


START_SECTION(virtual PepIterator& operator++())
	ptr = new TrypticIterator();
	TEST_EXCEPTION (Exception::InvalidIterator, ++(*ptr));
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("TrypticIterator_test.fasta"));
	ptr->begin();
	PepIterator & pepIt = ++(*ptr);
	TEST_EQUAL ((*pepIt).first,(**ptr).first);
	TEST_EQUAL ((*pepIt).second,(**ptr).second);
	pepIt = ++(*ptr);
	TEST_EQUAL ((*pepIt).first,(**ptr).first);
	TEST_EQUAL ((*pepIt).second,(**ptr).second);
END_SECTION

START_SECTION(virtual PepIterator* operator++(int i))
	ptr = new TrypticIterator();
	TEST_EXCEPTION (Exception::InvalidIterator, (*ptr)++);
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("TrypticIterator_test.fasta"));
	ptr->begin();
	FASTAEntry fe = **ptr;
	PepIterator * pepIt = (*ptr)++;
	TEST_EQUAL ((**pepIt).first,fe.first);
	TEST_EQUAL ((**pepIt).second,fe.second);
END_SECTION

START_SECTION(virtual bool begin())
	ptr = new TrypticIterator();
	TEST_EXCEPTION (Exception::InvalidIterator, (*ptr).begin());
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("TrypticIterator_test.fasta"));
	ptr->begin();
	FASTAEntry fe = **ptr;
	TEST_EQUAL(fe.first,">Entry 1");
	TEST_EQUAL(fe.second,"AAAAAK");
END_SECTION

START_SECTION(virtual bool isAtEnd())
	ptr = new TrypticIterator();
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("TrypticIterator_test.fasta"));
	ptr->begin();
	for (int i = 0; i < 13; ++i)
	{
		TEST_EQUAL(ptr->isAtEnd(), false);
		++(*ptr);
	}
	TEST_EQUAL(ptr->isAtEnd(), true);
END_SECTION


START_SECTION(virtual void setSpectrum(const std::vector< double > &))
	ptr = new TrypticIterator();
	const std::vector<double> spec;
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).setSpectrum(spec));
END_SECTION

START_SECTION(virtual const std::vector<double>& getSpectrum())
	ptr = new TrypticIterator();
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).getSpectrum());
END_SECTION

START_SECTION(virtual void setTolerance(double))
	ptr = new TrypticIterator();
	double t = 0.5;
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).setTolerance(t));
END_SECTION

START_SECTION(virtual double getTolerance())
	ptr = new TrypticIterator();
	TEST_EXCEPTION (Exception::NotImplemented, (*ptr).getTolerance());
END_SECTION

START_SECTION(virtual bool isDigestingEnd(char aa1, char aa2))
	ptr = new TrypticIterator();
	TEST_EQUAL (ptr->isDigestingEnd('R','C'),true);
	TEST_EQUAL (ptr->isDigestingEnd('K','C'),true);
	TEST_EQUAL (ptr->isDigestingEnd('R','P'),false);
	TEST_EQUAL (ptr->isDigestingEnd('K','P'),false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



