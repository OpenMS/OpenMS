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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
///////////////////////////

#define MULTI_LINE_STRING(...) #__VA_ARGS__ 

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IndexedMzMLDecoder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IndexedMzMLDecoder* ptr = nullptr;
IndexedMzMLDecoder* nullPointer = nullptr;
START_SECTION((IndexedMzMLDecoder()))
  ptr = new IndexedMzMLDecoder;
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~IndexedMzMLDecoder()))
  delete ptr;
END_SECTION

START_SECTION((int parseOffsets(String filename, std::streampos indexoffset, OffsetVector & spectra_offsets, OffsetVector& chromatograms_offsets)))
  // see also IndexedMzMLFile_test.cpp
  std::streampos res = IndexedMzMLDecoder().findIndexListOffset(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_NOT_EQUAL(res, -1)

  IndexedMzMLDecoder::OffsetVector spectra_offsets;
  IndexedMzMLDecoder::OffsetVector chromatograms_offsets;

  int off = IndexedMzMLDecoder().parseOffsets(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), res, spectra_offsets, chromatograms_offsets);

  TEST_EQUAL(off, 0)
  TEST_EQUAL(spectra_offsets.size(), 2)
  TEST_EQUAL(chromatograms_offsets.size(), 1)
END_SECTION

    
START_SECTION((std::streampos findIndexListOffset(String filename, int buffersize = 1023)))
  // see also IndexedMzMLFile_test.cpp
  //
  std::streampos res = IndexedMzMLDecoder().findIndexListOffset(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_NOT_EQUAL(res, -1)

  // A std::streamoff value of -1 is also used to represent error conditions by some of the I/O library functions. 
  std::streampos nonindex = IndexedMzMLDecoder().findIndexListOffset(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(nonindex, -1)

END_SECTION

    

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

