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
// $Authors: Guillaune Belz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/XMassFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

///////////////////////////

START_TEST(XMassFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

XMassFile* ptr = nullptr;
XMassFile* nullPointer = nullptr;
START_SECTION(XMassFile())
	ptr = new XMassFile;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~XMassFile())
	delete ptr;
END_SECTION

START_SECTION(template<typename SpectrumType> void load(const String& filename, MSSpectrum& spectrum) )
	TOLERANCE_ABSOLUTE(0.001)
	MSSpectrum s;
	MSSpectrum::ConstIterator it;
  TextFile::ConstIterator f_it;
	XMassFile f;

	TEST_EXCEPTION(Exception::FileNotFound, f.load("data_Idontexist", s);)

	f.load(OPENMS_GET_TEST_DATA_PATH("XMassFile_test/fid"),s);

	TEST_EQUAL(s.size(), 80478)
	ABORT_IF(s.size() != 80478)

  // read data for comparison
  TextFile file;
  file.load(OPENMS_GET_TEST_DATA_PATH("XMassFile_test_data.txt"));

  TEST_EQUAL((file.end() - file.begin()), 80478)
	ABORT_IF((file.end() - file.begin()) != 80478)

	for(it=s.begin(), f_it = file.begin(); it != s.end() && f_it != file.end(); ++it, ++f_it)
	{
    DoubleList test_values = ListUtils::create<double>(*f_it);
    ABORT_IF(test_values.size() != 2)

	  TEST_REAL_SIMILAR(it->getPosition()[0], test_values[0])
	  TEST_REAL_SIMILAR(it->getIntensity(), test_values[1])
	}

END_SECTION

START_SECTION(template<typename SpectrumType> void store(const String& filename, const MSSpectrum& spectrum) const)
  // not implemented
	TEST_EXCEPTION(Exception::NotImplemented, XMassFile().store(String(), MSSpectrum()))
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

