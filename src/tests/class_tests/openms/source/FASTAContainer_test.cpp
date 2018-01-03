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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FASTAContainer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

typedef FASTAContainer<TFI_Vector> FCVec;
typedef FASTAContainer<TFI_File> FCFile;

FCVec* ptr = nullptr;
FCVec* nullPointer = nullptr;
START_SECTION(FASTAContainer())
{
  ptr = new FCVec(std::vector<FASTAFile::FASTAEntry>());
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~FASTAContainer())
{
  delete ptr;
}
END_SECTION

std::vector<FASTAFile::FASTAEntry> fev = { {"id0", "desc0", "AAAA"},{ "id1", "desc1", "BBBB" },{ "id2", "desc2", "CCCC" },{ "id3", "desc3", "DDDD" } };

START_SECTION(FASTAContainer(const String& FASTA_file))
  FCFile f(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"));
  TEST_EQUAL(f.size(), 0)
END_SECTION

START_SECTION(FASTAContainer(std::vector<FASTAFile::FASTAEntry>& data))
  FCVec fv(fev);
  TEST_EQUAL(fv.size(), 4)
END_SECTION

START_SECTION(size_t getChunkOffset() const)
  // FCFile: tested below
  FCVec fv(fev);
  TEST_EQUAL(fv.getChunkOffset(), 0)
END_SECTION


START_SECTION(bool activateCache())
  // FCFile: tested below
  FCVec fv(fev);
  TEST_EQUAL(fv.activateCache(), 1)
  TEST_EQUAL(fv.activateCache(), 0)
END_SECTION

START_SECTION(bool cacheChunk(int suggested_size))
  // FCFile: tested below
  FCVec fv(fev);
  TEST_EQUAL(fv.cacheChunk(333), 1)
  TEST_EQUAL(fv.cacheChunk(333), 0)
END_SECTION

START_SECTION(size_t chunkSize() const)
  // FCFile: tested below
  FCVec fv(fev);
  TEST_EQUAL(fv.chunkSize(), 4)
END_SECTION

START_SECTION(const FASTAFile::FASTAEntry& chunkAt(size_t pos) const)
  // FCFile: tested below
  FCVec fv(fev);
  FASTAFile::FASTAEntry pe = fv.chunkAt(3);
  TEST_EQUAL(pe.identifier, "id3");
END_SECTION

START_SECTION(bool readAt(FASTAFile::FASTAEntry& protein, size_t pos))
  // FCFile: tested below
  FCVec fv(fev);
  FASTAFile::FASTAEntry pe;
  TEST_EQUAL(fv.readAt(pe, 3), true);
  TEST_EQUAL(pe.identifier, "id3");
END_SECTION

START_SECTION(bool empty() const)
  FCFile f(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"));
  TEST_EQUAL(f.empty(), false)
  FCFile f2(OPENMS_GET_TEST_DATA_PATH("degenerate_cases/empty.fasta"));
  TEST_EQUAL(f2.empty(), true)
  FCVec fv(fev);
  TEST_EQUAL(fv.empty(), false);
  std::vector<FASTAFile::FASTAEntry> feve;
  FCVec fv2(feve);
  TEST_EQUAL(fv2.empty(), true);
END_SECTION

START_SECTION(size_t size() const)
  FCFile f(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"));
  TEST_EQUAL(f.cacheChunk(2), true)
  TEST_EQUAL(f.size(), 2)
  TEST_EQUAL(f.activateCache(), true)
  TEST_EQUAL(f.size(), 2)
  FASTAFile::FASTAEntry pe, pe2;
  TEST_EQUAL(f.readAt(pe, 0), true);
  pe2 = f.chunkAt(0);
  TEST_EQUAL(pe == pe2, true)
  TEST_EQUAL(pe.description, "This is the description of the first protein")
  pe2 = f.chunkAt(1);
  TEST_EQUAL(pe == pe2, false)
  TEST_EQUAL(pe2.description, "This is the description of the second protein")

  // read next chunk, and re-read from disk again, using byte offsets
  TEST_EQUAL(f.cacheChunk(1), true)
  TEST_EQUAL(f.activateCache(), true)
  TEST_EQUAL(f.readAt(pe, 0), true); // third global entry
  TEST_EQUAL(pe.identifier, "P68509|1433F_BOVIN")
  TEST_EQUAL(pe.description, "This is the description of the first protein")
  
  // read until end
  TEST_EQUAL(f.cacheChunk(3), true)  // only 2 can be read, but thats ok
  TEST_EQUAL(f.activateCache(), true)
  TEST_EQUAL(f.chunkSize(), 2)
  pe = f.chunkAt(1);
  TEST_EQUAL(pe.description, " ##0")
  TEST_EQUAL(f.readAt(pe2, 4), true);
  TEST_EQUAL(pe == pe2, true)
  
  // reached the end after 5 entries
  TEST_EQUAL(f.cacheChunk(3), false)
  TEST_EQUAL(f.chunkSize(), 2)
  TEST_EQUAL(f.activateCache(), false)
  TEST_EQUAL(f.chunkSize(), 0)
  TEST_EQUAL(f.cacheChunk(3), false)
  TEST_EQUAL(f.activateCache(), false)

  // read from disk again after reaching EOF, using byte offsets
  TEST_EQUAL(f.readAt(pe, 0), true);
  TEST_EQUAL(pe.identifier, "P68509|1433F_BOVIN")
  TEST_EQUAL(pe.description, "This is the description of the first protein")
  TEST_EQUAL(f.readAt(pe, 4), true);
  TEST_EQUAL(pe.identifier, "test")
  TEST_EQUAL(pe.description, " ##0")

  FCVec fv(fev);
  TEST_EQUAL(fv.size(), 4);
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



