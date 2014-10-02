// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: $
// $Authors: Adam Tenderholt $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include "test_config.h"

///////////////////////////

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/FORMAT/TarFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TarFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TarFile* ptr = 0;
TarFile* nullPointer = 0;
START_SECTION((TarFile()))
  ptr = new TarFile;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~TarFile()))
	delete ptr;
END_SECTION

TarFile file;
MSExperiment<> empty;
MSExperiment<> full;
for (Size i = 0; i < 4; i++)
{
  MSSpectrum<Peak1D> spectrum;
  spectrum.setRT(0.25 * i);
  full.addSpectrum(spectrum);
}

// this tests the expected case
START_SECTION((void load(const String& filename, MSExperiment<Peak1D>& experiment)))
{
  MSExperiment<Peak1D> expr(full);
  file.load(PEAKINVESTIGATOR_GET_TEST_DATA_PATH("TarFile_1.tar.gz"), expr);

  TEST_EQUAL(expr.size(), 4);
  for (Size i = 0; i < 4; i++)
  {
    MSSpectrum<Peak1D> current = expr[i];

    TEST_EQUAL(current.size(), 14);
    for (Size j = 0; j < current.size(); j++)
    {
      TEST_REAL_SIMILAR(current[j].getMZ(), j + 1);
      TEST_REAL_SIMILAR(current[j].getIntensity(), 15.0 - double (j + 1) / double (i + 1));
    }
  }
}
END_SECTION

START_SECTION([EXTRA] load with invalid filename)
{
  MSExperiment<> expr(empty);
  file.load("", expr);
  TEST_EQUAL(expr.size(), 0);
}
END_SECTION

START_SECTION([EXTRA] load with empty experiment)
{
  MSExperiment<> expr(empty);
  file.load(PEAKINVESTIGATOR_GET_TEST_DATA_PATH("TarFile_1.tar.gz"), expr);
  TEST_EQUAL(expr.size(), 0);
}
END_SECTION

START_SECTION([EXTRA] load with no data)
{
  MSExperiment<> expr(full);
  file.load(PEAKINVESTIGATOR_GET_TEST_DATA_PATH("TarFile_2_empty.tar.gz"), expr);

  TEST_EQUAL(expr.size(), 4);
  for (Size i = 0; i < expr.size(); i++)
  {
    TEST_EQUAL(expr[i].size(), 0);
  }
}
END_SECTION

// test load where scan_000003.txt has been replaced with scan_000004.txt
START_SECTION([EXTRA] load with scan mismatch)
{
  MSExperiment<> expr(full);
  file.load(PEAKINVESTIGATOR_GET_TEST_DATA_PATH("TarFile_3_mismatch.tar.gz"), expr);

  TEST_EQUAL(expr.size(), 4);
  for (Size i = 0; i < expr.size(); i++)
  {
    MSSpectrum<Peak1D> current = expr[i];

    if (i < expr.size() - 1)
    {
      TEST_EQUAL(current.size(), 14);
    }
    else
    {
      TEST_EQUAL(current.size(), 0);
    }

    for (Size j = 0; j < current.size(); j++)
    {
      TEST_REAL_SIMILAR(current[j].getMZ(), j + 1);
      TEST_REAL_SIMILAR(current[j].getIntensity(), 15.0 - double (j + 1) / double (i + 1));
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

