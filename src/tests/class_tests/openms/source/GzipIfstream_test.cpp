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
#include <OpenMS/FORMAT/GzipIfstream.h>
using namespace OpenMS;



///////////////////////////

START_TEST(GzipIfstream, "$Id$")

GzipIfstream* ptr = nullptr;
GzipIfstream* nullPointer = nullptr;
START_SECTION((GzipIfstream()))
	ptr = new GzipIfstream;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~GzipIfstream()))
	delete ptr;
END_SECTION

START_SECTION(GzipIfstream(const char * filename))
	TEST_EXCEPTION(Exception::FileNotFound, GzipIfstream gzip2(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	
	GzipIfstream gzip(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
	
	TEST_EQUAL(gzip.streamEnd(), false)
	TEST_EQUAL(gzip.isOpen(),true)
	char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	TEST_EQUAL(29, gzip.read(buffer, len))
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))

END_SECTION

START_SECTION(void open(const char *filename))
	GzipIfstream gzip;
	TEST_EXCEPTION(Exception::FileNotFound, gzip.open(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	
	gzip.open(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
	
	TEST_EQUAL(gzip.streamEnd(), false)
	TEST_EQUAL(gzip.isOpen(),true)
	char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	TEST_EQUAL(29, gzip.read(buffer, len))
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))
	
END_SECTION

START_SECTION(size_t read(char *s, size_t n))
	//tested in open(const char * filename)
	GzipIfstream gzip(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1_corrupt.gz"));
		char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	//gzip.updateCRC32(buffer,10);
	//Why does that throw a "naked" Exception instead of a ConversionError?
	//~ TEST_EXCEPTION(Exception::BaseException,gzip.read(&buffer[9],10));
//	gzip.updateCRC32(&buffer[9],19);
//	TEST_EQUAL(gzip.isCorrupted(),true)
	
	GzipIfstream gzip2(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
	TEST_EQUAL(gzip2.isOpen(),true)
	gzip2.read(buffer, len);
	TEST_EQUAL(1, gzip2.read(buffer,10));
	TEST_EQUAL(gzip2.isOpen(), false)
	TEST_EQUAL(gzip2.streamEnd(),true)
	
	gzip2.open(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1_corrupt.gz"));

	 //gzip2.updateCRC32(buffer,(size_t)30);
	//TEST_EQUAL(gzip2.isCorrupted(),true )
	gzip2.close();
	TEST_EQUAL(gzip2.isOpen(), false)
	TEST_EQUAL(gzip2.streamEnd(),true)
	TEST_EXCEPTION(Exception::IllegalArgument, gzip2.read(buffer,10))
	gzip2.close();			
	TEST_EQUAL(gzip2.isOpen(), false)
	TEST_EQUAL(gzip2.streamEnd(),true)
	TEST_EXCEPTION(Exception::IllegalArgument, gzip2.read(buffer,10))
	gzip2.open(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
	
		TEST_EQUAL(5, gzip2.read(buffer, 5))
					//			gzip2.updateCRC32(buffer,5);
		TEST_EQUAL(5, gzip2.read(&buffer[5], 5))
						//gzip2.updateCRC32(&buffer[5],5);
		TEST_EQUAL(5, gzip2.read(&buffer[10], 5))
				//		gzip2.updateCRC32(&buffer[10],5);
		TEST_EQUAL(5, gzip2.read(&buffer[15], 5))
					//				gzip2.updateCRC32(&buffer[15],5);
		TEST_EQUAL(5, gzip2.read(&buffer[20], 5))
						//	gzip2.updateCRC32(&buffer[20],5);
					TEST_EQUAL(4, gzip2.read(&buffer[25], 4))
						//			gzip2.updateCRC32(&buffer[25],4);
					char end_of_file[1];
					TEST_EQUAL(1,gzip2.read(end_of_file,2))
		//							gzip2.updateCRC32(end_of_file,1);
					TEST_EQUAL(gzip2.streamEnd(),true)				
					buffer[29]= '\0';				
//	TEST_EQUAL(gzip2.isCorrupted(),false)
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))
END_SECTION

START_SECTION(void close())
	//tested in read
	NOT_TESTABLE
END_SECTION
START_SECTION(bool streamEnd() const )
	//!!!tested in open(const char * filename) and read
	NOT_TESTABLE
END_SECTION
START_SECTION(bool isOpen() const)
	//tested in open(const char * filename) and read
	NOT_TESTABLE
END_SECTION
/*
(updateCRC32(char* s, size_t n))
	//tested in open(const char * filename) and read
	_TESTABLE
_SECTION
_SECTION(isCorrupted())
	//tested in open(const char * filename) and read
	TESTABLE
_SECTION*/

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
