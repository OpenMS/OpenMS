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
// $Maintainer: Timo Sachsenberg$
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>

using namespace OpenMS;
using namespace std;


///////////////////////////

START_TEST(Bzip2Ifstream_test, "$Id$")

Bzip2Ifstream* ptr = nullptr;
Bzip2Ifstream* nullPointer = nullptr;
START_SECTION((Bzip2Ifstream()))
	ptr = new Bzip2Ifstream;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~Bzip2Ifstream()))
	delete ptr;
END_SECTION

START_SECTION(Bzip2Ifstream(const char * filename))
	TEST_EXCEPTION(Exception::FileNotFound, Bzip2Ifstream bzip2(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	
	Bzip2Ifstream bzip(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	
	TEST_EQUAL(bzip.streamEnd(), false)
	TEST_EQUAL(bzip.isOpen(),true)
	char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	TEST_EQUAL(29, bzip.read(buffer, len))
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))

END_SECTION

START_SECTION(void open(const char *filename))
	Bzip2Ifstream bzip;
	TEST_EXCEPTION(Exception::FileNotFound, bzip.open(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	
	bzip.open(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	
	TEST_EQUAL(bzip.streamEnd(), false)
	TEST_EQUAL(bzip.isOpen(),true)
	char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	TEST_EQUAL(29, bzip.read(buffer, len))
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))
	
END_SECTION

START_SECTION(size_t read(char *s, size_t n))
	//tested in open(const char * filename)
	Bzip2Ifstream bzip(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1_corrupt.bz2"));
		char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	TEST_EXCEPTION(Exception::ParseError, bzip.read(buffer,10))
	
	Bzip2Ifstream bzip2(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	bzip2.read(buffer, len);
	TEST_EQUAL(1, bzip2.read(buffer,10));
	TEST_EQUAL(bzip2.isOpen(), false)
	TEST_EQUAL(bzip2.streamEnd(),true)
	
	bzip2.open(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1_corrupt.bz2"));
	TEST_EXCEPTION(Exception::ParseError, bzip2.read(buffer,10))
	bzip2.close();
	TEST_EQUAL(bzip2.isOpen(), false)
	TEST_EQUAL(bzip2.streamEnd(),true)
	TEST_EXCEPTION(Exception::IllegalArgument, bzip2.read(buffer,10))
	bzip2.close();
	TEST_EQUAL(bzip2.isOpen(), false)
	TEST_EQUAL(bzip2.streamEnd(),true)
	TEST_EXCEPTION(Exception::IllegalArgument, bzip2.read(buffer,10))
	bzip2.open(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	TEST_EQUAL(29, bzip2.read(buffer, len))
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))
END_SECTION

START_SECTION(void close())
	//tested in read
	NOT_TESTABLE
END_SECTION
START_SECTION(bool streamEnd() const )
	//tested in open(const char * filename) and read
	NOT_TESTABLE
END_SECTION
START_SECTION(bool isOpen() const )
	//tested in open(const char * filename) and read
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
