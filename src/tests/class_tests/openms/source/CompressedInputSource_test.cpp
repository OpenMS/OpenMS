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
#include <OpenMS/FORMAT/CompressedInputSource.h>
#include <OpenMS/FORMAT/GzipInputStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
using namespace OpenMS;


///////////////////////////

START_TEST(CompressedInputSource, "$Id$")

CompressedInputSource* ptr = nullptr;
CompressedInputSource* nullPointer = nullptr;
START_SECTION(CompressedInputSource(const String& file_path, const char * header,xercesc::MemoryManager* const manager = xercesc::XMLPlatformUtils::fgMemoryManager))
  char header[3];
  header[0] = 'B';
  header[1] = 'Z';
  header[2] = '\0';
  String bz = String(header);
  ptr = new CompressedInputSource(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"), bz);
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~CompressedInputSource()))
  delete ptr;
END_SECTION

START_SECTION(CompressedInputSource(const XMLCh *const file_path, const char *header, xercesc::MemoryManager *const manager=xercesc::XMLPlatformUtils::fgMemoryManager))
  char header[3];
  header[0] = 'B';
  header[1] = 'Z';
  header[2] = '\0';
  String bz = String(header);
  String filename(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
  ptr = new CompressedInputSource(Internal::StringManager().convert(filename.c_str()).c_str(), bz);
  TEST_NOT_EQUAL(ptr, nullPointer)
  delete ptr;
END_SECTION


START_SECTION(virtual xercesc::BinInputStream* makeStream() const)
  char header[3];
  header[0] = 'B';
  header[1] = 'Z';
  header[2] = '\0';
  String bz = String(header);
  CompressedInputSource source(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist"), bz);
  TEST_EXCEPTION(Exception::FileNotFound,source.makeStream())
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
