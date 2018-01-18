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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/MS2File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////
///////////////////////////

START_TEST(MS2File, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MS2File * ptr = nullptr;
MS2File* nullPointer = nullptr;
START_SECTION((MS2File()))
ptr = new MS2File;
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MS2File()))
delete ptr;
END_SECTION

TOLERANCE_ABSOLUTE(0.01)

START_SECTION((template <typename MapType> void load(const String &filename, MapType & exp)))
MS2File file;
PeakMap exp;
file.load(OPENMS_GET_TEST_DATA_PATH("MS2File_test_spectra.ms2"), exp);

//test DocumentIdentifier addition
TEST_STRING_EQUAL(exp.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("MS2File_test_spectra.ms2"));
TEST_STRING_EQUAL(FileTypes::typeToName(exp.getLoadedFileType()), "ms2");

TEST_EQUAL(exp.size(), 2)

TEST_EQUAL(exp[0].size(), 4)
TEST_EQUAL(exp[1].size(), 4)

TEST_STRING_EQUAL(exp[0].getNativeID(), "index=0")
TEST_STRING_EQUAL(exp[1].getNativeID(), "index=1")

TEST_REAL_SIMILAR(exp[0].getPrecursors()[0].getMZ(), 444.44)
TEST_REAL_SIMILAR(exp[1].getPrecursors()[0].getMZ(), 555.555)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
