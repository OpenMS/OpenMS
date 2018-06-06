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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/SourceFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SourceFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SourceFile* ptr = nullptr;
SourceFile* nullPointer = nullptr;
START_SECTION((SourceFile()))
	ptr = new SourceFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~SourceFile()))
	delete ptr;
END_SECTION

START_SECTION((float getFileSize() const))
  SourceFile tmp;
  TEST_EQUAL(tmp.getFileSize(),0);
END_SECTION

START_SECTION((void setFileSize(float file_size)))
  SourceFile tmp;
	tmp.setFileSize(1.667f);
  TEST_REAL_SIMILAR(tmp.getFileSize(),1.667f);
END_SECTION

START_SECTION((const String& getFileType() const))
  SourceFile tmp;
  TEST_EQUAL(tmp.getFileType(), "");
END_SECTION

START_SECTION((void setFileType(const String& file_type)))
  SourceFile tmp;
	tmp.setFileType("PEAKDATA");
  TEST_EQUAL(tmp.getFileType(), "PEAKDATA");
END_SECTION

START_SECTION((const String& getNameOfFile() const))
  SourceFile tmp;
  TEST_EQUAL(tmp.getNameOfFile(),"");
END_SECTION

START_SECTION((void setNameOfFile(const String& name_of_file)))
  SourceFile tmp;
  tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
  TEST_EQUAL(tmp.getNameOfFile(),"The White Stripes - Ball and Biscuit");
END_SECTION

START_SECTION((const String& getPathToFile() const))
  SourceFile tmp;
  TEST_EQUAL(tmp.getPathToFile(),"");
END_SECTION

START_SECTION((void setPathToFile(const String& path_path_to_file)))
  SourceFile tmp;
  tmp.setPathToFile("/misc/sturm/mp3/");
  TEST_EQUAL(tmp.getPathToFile(),"/misc/sturm/mp3/");
END_SECTION

START_SECTION((const String& getChecksum() const))
  SourceFile tmp;
  TEST_EQUAL(tmp.getChecksum(), "");
END_SECTION

START_SECTION(ChecksumType getChecksumType() const)
  SourceFile tmp;
  TEST_EQUAL(tmp.getChecksumType(), SourceFile::UNKNOWN_CHECKSUM);
END_SECTION

START_SECTION((void setChecksum(const String& checksum, ChecksumType type)))
  SourceFile tmp;
  tmp.setChecksum("2fd4e1c67a2d28fced849ee1bb76e7391b93eb12",SourceFile::SHA1);
  TEST_EQUAL(tmp.getChecksum(), "2fd4e1c67a2d28fced849ee1bb76e7391b93eb12");
  TEST_EQUAL(tmp.getChecksumType(), SourceFile::SHA1);
END_SECTION

START_SECTION((const String& getNativeIDType() const))
	SourceFile tmp;
	TEST_STRING_EQUAL(tmp.getNativeIDType(), "");
END_SECTION

START_SECTION((void setNativeIDType(const String& type)))
  SourceFile tmp;
  tmp.setNativeIDType("bla");
  TEST_STRING_EQUAL(tmp.getNativeIDType(), "bla");
END_SECTION

START_SECTION((SourceFile(const SourceFile& source)))
	SourceFile tmp;
	tmp.setFileType("CALIBRATIONINFO");
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	tmp.setPathToFile("/misc/sturm/mp3/");
	tmp.setChecksum("2fd4e1c67a2d28fced849ee1bb76e7391b93eb12", SourceFile::MD5);
	tmp.setMetaValue("bla",4.0);
	
	SourceFile tmp2(tmp);
	TEST_EQUAL(tmp2.getFileType(), "CALIBRATIONINFO");
	TEST_EQUAL(tmp2.getNameOfFile(),"The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp2.getPathToFile(),"/misc/sturm/mp3/");
	TEST_EQUAL(tmp2.getChecksum(), "2fd4e1c67a2d28fced849ee1bb76e7391b93eb12");
	TEST_EQUAL(tmp2.getChecksumType(), SourceFile::MD5);
	TEST_REAL_SIMILAR(tmp2.getMetaValue("bla"), 4.0);
END_SECTION

START_SECTION((SourceFile& operator= (const SourceFile& source)))
	SourceFile tmp;
	tmp.setFileType("PUBLICATION");
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	tmp.setPathToFile("/misc/sturm/mp3/");
	tmp.setChecksum("2fd4e1c67a2d28fced849ee1bb76e7391b93eb12", SourceFile::MD5);
	tmp.setMetaValue("bla",4.0);
	
	//normal assignment
	SourceFile tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getFileType(),"PUBLICATION");
	TEST_EQUAL(tmp2.getNameOfFile(),"The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp2.getPathToFile(),"/misc/sturm/mp3/");
	TEST_EQUAL(tmp2.getChecksum(),"2fd4e1c67a2d28fced849ee1bb76e7391b93eb12");
	TEST_EQUAL(tmp2.getChecksumType(), SourceFile::MD5);
	TEST_REAL_SIMILAR(tmp2.getMetaValue("bla"), 4.0);
	
	//assignment of empty object
	tmp2 = SourceFile();
	TEST_EQUAL(tmp2.getFileType(), "");
	TEST_EQUAL(tmp2.getNameOfFile(),"");
	TEST_EQUAL(tmp2.getPathToFile(),"");
	TEST_EQUAL(tmp2.getChecksum(),"");
	TEST_EQUAL(tmp2.getChecksumType(), SourceFile::UNKNOWN_CHECKSUM);
	TEST_EQUAL(tmp2.metaValueExists("bla"), false);
END_SECTION

START_SECTION((bool operator== (const SourceFile& rhs) const))
	SourceFile tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp2.setFileType("PARAMETERSFILE");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setChecksum("", SourceFile::MD5);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
	tmp.setMetaValue("bla",4.0);
	TEST_EQUAL(tmp==tmp2, false);	
	
	tmp2 = tmp;	
	tmp.setPathToFile("/misc/sturm/mp3/");
	TEST_EQUAL(tmp==tmp2, false);
END_SECTION

START_SECTION((bool operator!= (const SourceFile& rhs) const))
	SourceFile tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp2.setFileType("MISC");
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setChecksum("2fd4e1c67a2d28fced849ee1bb76e7391b93eb12", SourceFile::UNKNOWN_CHECKSUM);
	TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;
	tmp.setMetaValue("bla",4.0);
	TEST_EQUAL(tmp!=tmp2, true);	

	tmp2 = tmp;	
	tmp.setPathToFile("/misc/sturm/mp3/");
	TEST_EQUAL(tmp!=tmp2, true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



