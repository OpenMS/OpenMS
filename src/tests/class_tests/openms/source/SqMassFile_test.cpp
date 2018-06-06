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

#include <OpenMS/FORMAT/SqMassFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <QFile>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SqMassFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SqMassFile* ptr = nullptr;
SqMassFile* nullPointer = nullptr;
START_SECTION((SqMassFile()))
  ptr = new SqMassFile;
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~SqMassFile()))
  delete ptr;
END_SECTION

TOLERANCE_RELATIVE(1.0005)

START_SECTION(void load(const String& filename, MapType& map))
{
  SqMassFile file;
  MSExperiment exp;
  file.load(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"), exp);

  MSExperiment exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp2);

  TEST_EQUAL(exp.getNrSpectra(), exp2.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp2.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp2.getSpectra()[0], false) // no exact duplicate

  cout.precision(17);

  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.
  TOLERANCE_ABSOLUTE(1e-4)
  TOLERANCE_RELATIVE(1.001) // 0.1 % error for intensity 

  for (Size i = 0; i < exp.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp.getSpectrum(i).size(), exp2.getSpectra()[i].size())
    for (Size k = 0; k < exp.getSpectrum(i).size(); k++)
    {
      // slof is no good for values smaller than 5
      // if (exp.getSpectrum(i)[k].getIntensity() < 1.0) {continue;} 
      TEST_REAL_SIMILAR(exp.getSpectrum(i)[k].getIntensity(), exp2.getSpectra()[i][k].getIntensity())
    }
  }

  for (Size i = 0; i < exp.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp.getChromatogram(i).size() == exp2.getChromatograms()[i].size(), true)
    for (Size k = 0; k < exp.getChromatogram(i).size(); k++)
    {
      TEST_REAL_SIMILAR(exp.getChromatogram(i)[k].getIntensity(), exp2.getChromatograms()[i][k].getIntensity())
    }
  }

  TOLERANCE_ABSOLUTE(1e-5)
  TOLERANCE_RELATIVE(1.000001) // less than 1ppm error for m/z 
  for (Size i = 0; i < exp.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp.getSpectrum(i).size(), exp2.getSpectra()[i].size())
    for (Size k = 0; k < exp.getSpectrum(i).size(); k++)
    {
      TEST_REAL_SIMILAR(exp.getSpectrum(i)[k].getMZ(), exp2.getSpectra()[i][k].getMZ())
    }

  }
  TOLERANCE_ABSOLUTE(0.05) // max 0.05 seconds error in RT
  for (Size i = 0; i < exp.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp.getChromatogram(i).size() == exp2.getChromatograms()[i].size(), true)
    for (Size k = 0; k < exp.getChromatogram(i).size(); k++)
    {
      TEST_REAL_SIMILAR(exp.getChromatogram(i)[k].getRT(), exp2.getChromatograms()[i][k].getRT())
    }
  }

  // no 1:1 mapping of experimental settings ...
  TEST_EQUAL(exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, false)
}
END_SECTION

// reset error tolerances to default values
TOLERANCE_ABSOLUTE(1e-5)
TOLERANCE_RELATIVE(1+1e-5)

START_SECTION(void store(const String& filename, MapType& map))
{
  MSExperiment exp_orig;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), exp_orig);

  SqMassFile::SqMassConfig config;
  config.use_lossy_numpress = false;
  config.linear_fp_mass_acc = -1;
  config.write_full_meta = false;

  SqMassFile file;
  file.setConfig(config);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  std::cout << "Storing in file " << tmp_filename << std::endl;
  file.store(tmp_filename, exp_orig);

  MSExperiment exp;
  file.load(tmp_filename, exp);

  TEST_EQUAL(exp.getNrSpectra(), exp_orig.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp_orig.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp_orig.getSpectra()[0], false) // no exact duplicate

  MSExperiment exp2 = exp_orig;

  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.
  TOLERANCE_ABSOLUTE(1e-4)
  TOLERANCE_RELATIVE(1.001) // 0.1 % error for intensity 

  // reset error tolerances to default values
  TOLERANCE_ABSOLUTE(1e-5)
  TOLERANCE_RELATIVE(1+1e-5)

  // since we specified no lossy compression, we expect high accuracy
  TOLERANCE_ABSOLUTE(1e-8)
  TOLERANCE_RELATIVE(1.00000001) 

  for (Size i = 0; i < exp.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp.getSpectrum(i).size(), exp2.getSpectra()[i].size())
    for (Size k = 0; k < exp.getSpectrum(i).size(); k++)
    {
      // slof is no good for values smaller than 5
      // if (exp.getSpectrum(i)[k].getIntensity() < 1.0) {continue;} 
      TEST_REAL_SIMILAR(exp.getSpectrum(i)[k].getIntensity(), exp2.getSpectra()[i][k].getIntensity())
    }

  }

  for (Size i = 0; i < exp.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp.getChromatogram(i).size() == exp2.getChromatograms()[i].size(), true)
    for (Size k = 0; k < exp.getChromatogram(i).size(); k++)
    {
      TEST_REAL_SIMILAR(exp.getChromatogram(i)[k].getIntensity(), exp2.getChromatograms()[i][k].getIntensity())
    }
  }

  for (Size i = 0; i < exp.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp.getSpectrum(i).size(), exp2.getSpectra()[i].size())
    for (Size k = 0; k < exp.getSpectrum(i).size(); k++)
    {
      TEST_REAL_SIMILAR(exp.getSpectrum(i)[k].getMZ(), exp2.getSpectra()[i][k].getMZ())
    }
  }

  for (Size i = 0; i < exp.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp.getChromatogram(i).size() == exp2.getChromatograms()[i].size(), true)
    for (Size k = 0; k < exp.getChromatogram(i).size(); k++)
    {
      TEST_REAL_SIMILAR(exp.getChromatogram(i)[k].getRT(), exp2.getChromatograms()[i][k].getRT())
    }
  }

  // no 1:1 mapping of experimental settings ...
  TEST_EQUAL(exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, false)
}
END_SECTION

START_SECTION([EXTRA] void store(const String& filename, MapType& map))
{
  MSExperiment exp_orig;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), exp_orig);

  SqMassFile::SqMassConfig config;
  config.use_lossy_numpress = true;
  config.linear_fp_mass_acc = 0.0001;
  config.write_full_meta = true;

  SqMassFile file;
  file.setConfig(config);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  std::cout << "Storing in file " << tmp_filename << std::endl;
  file.store(tmp_filename, exp_orig);

  MSExperiment exp;
  file.load(tmp_filename, exp);

  TEST_EQUAL(exp.getNrSpectra(), exp_orig.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp_orig.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp_orig.getSpectra()[0], false) // no exact duplicate

  MSExperiment exp2 = exp_orig;

  // using full meta should give 1:1 mapping of experimental settings ...
  TEST_EQUAL(exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, true)

  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.
  TOLERANCE_ABSOLUTE(1e-4)
  TOLERANCE_RELATIVE(1.001) // 0.1 % error for intensity 

  for (Size i = 0; i < exp.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp.getSpectrum(i).size(), exp2.getSpectra()[i].size())
    for (Size k = 0; k < exp.getSpectrum(i).size(); k++)
    {
      // slof is no good for values smaller than 5
      // if (exp.getSpectrum(i)[k].getIntensity() < 1.0) {continue;} 
      TEST_REAL_SIMILAR(exp.getSpectrum(i)[k].getIntensity(), exp2.getSpectra()[i][k].getIntensity())
    }

  }

  for (Size i = 0; i < exp.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp.getChromatogram(i).size() == exp2.getChromatograms()[i].size(), true)
    for (Size k = 0; k < exp.getChromatogram(i).size(); k++)
    {
      TEST_REAL_SIMILAR(exp.getChromatogram(i)[k].getIntensity(), exp2.getChromatograms()[i][k].getIntensity())
    }
  }

  TOLERANCE_ABSOLUTE(1e-5)
  TOLERANCE_RELATIVE(1.000001) // less than 1ppm error for m/z 
  for (Size i = 0; i < exp.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp.getSpectrum(i).size(), exp2.getSpectra()[i].size())
    for (Size k = 0; k < exp.getSpectrum(i).size(); k++)
    {
      TEST_REAL_SIMILAR(exp.getSpectrum(i)[k].getMZ(), exp2.getSpectra()[i][k].getMZ())
    }
  }
  TOLERANCE_ABSOLUTE(0.05) // max 0.05 seconds error in RT
  for (Size i = 0; i < exp.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp.getChromatogram(i).size() == exp2.getChromatograms()[i].size(), true)
    for (Size k = 0; k < exp.getChromatogram(i).size(); k++)
    {
      TEST_REAL_SIMILAR(exp.getChromatogram(i)[k].getRT(), exp2.getChromatograms()[i][k].getRT())
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

