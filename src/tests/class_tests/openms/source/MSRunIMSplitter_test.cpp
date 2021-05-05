// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/IONMOBILITY/MSRunIMSplitter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MSRunIMSplitter, "$Id$")

/////////////////////////////////////////////////////////////

MSRunIMSplitter* e_ptr = nullptr;
MSRunIMSplitter* e_nullPointer = nullptr;

START_SECTION((MSRunIMSplitter()))
	e_ptr = new MSRunIMSplitter;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~MSRunIMSplitter()))
	delete e_ptr;
END_SECTION

e_ptr = new MSRunIMSplitter();

START_SECTION((std::vector<PeakMap> splitByFAIMSCV(PeakMap& exp)))
	delete e_ptr;
	e_ptr = new MSRunIMSplitter();
  MzMLFile IM_file;
  PeakMap exp;
  IM_file.load(OPENMS_GET_TEST_DATA_PATH("IM_FAIMS_test.mzML"), exp);

  TEST_EQUAL(exp.getSpectra().size(), 19)

  vector<PeakMap> splitPeakMap = e_ptr->splitByFAIMSCV(std::move(exp));
  TEST_EQUAL(splitPeakMap.size(), 3)

	TEST_EQUAL(splitPeakMap[0].size(), 4)
	TEST_EQUAL(splitPeakMap[1].size(), 9)
	TEST_EQUAL(splitPeakMap[2].size(), 6)

	for (PeakMap::Iterator it = splitPeakMap[0].begin(); it != splitPeakMap[0].end(); ++it)
	{
		TEST_EQUAL(it->getDriftTime(), -65.0)
	}
	for (PeakMap::Iterator it = splitPeakMap[1].begin(); it != splitPeakMap[1].end(); ++it)
	{
		TEST_EQUAL(it->getDriftTime(), -55.0)
	}
	for (PeakMap::Iterator it = splitPeakMap[2].begin(); it != splitPeakMap[2].end(); ++it)
	{
		TEST_EQUAL(it->getDriftTime(), -45.0)
	}

	TEST_EQUAL(splitPeakMap[1].getExperimentalSettings().getDateTime().toString(), "2019-09-07T09:40:04")

END_SECTION


delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
