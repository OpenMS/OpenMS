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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>
///////////////////////////

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <QFile>

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

void cmpDataIntensity(MSExperiment& exp1, MSExperiment& exp2, double abs_tol = 1e-5, double rel_tol = 1+1e-5)
{
  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.

  TOLERANCE_ABSOLUTE(abs_tol)
  TOLERANCE_RELATIVE(rel_tol)
  for (Size i = 0; i < exp1.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp1.getSpectrum(i).size(), exp2.getSpectra()[i].size())
    for (Size k = 0; k < exp1.getSpectrum(i).size(); k++)
    {
      // slof is no good for values smaller than 5
      // if (exp.getSpectrum(i)[k].getIntensity() < 1.0) {continue;} 
      TEST_REAL_SIMILAR(exp1.getSpectrum(i)[k].getIntensity(), exp2.getSpectra()[i][k].getIntensity())
    }
  }

  for (Size i = 0; i < exp1.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp1.getChromatogram(i).size() == exp2.getChromatograms()[i].size(), true)
    for (Size k = 0; k < exp1.getChromatogram(i).size(); k++)
    {
      TEST_REAL_SIMILAR(exp1.getChromatogram(i)[k].getIntensity(), exp2.getChromatograms()[i][k].getIntensity())
    }
  }
}

void cmpDataMZ(MSExperiment& exp1, MSExperiment& exp2, double abs_tol = 1e-5, double rel_tol = 1+1e-5)
{
  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.

  TOLERANCE_ABSOLUTE(abs_tol)
  TOLERANCE_RELATIVE(rel_tol)
  for (Size i = 0; i < exp1.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp1.getSpectrum(i).size(), exp2.getSpectra()[i].size())
    for (Size k = 0; k < exp1.getSpectrum(i).size(); k++)
    {
      // slof is no good for values smaller than 5
      // if (exp.getSpectrum(i)[k].getIntensity() < 1.0) {continue;} 
      TEST_REAL_SIMILAR(exp1.getSpectrum(i)[k].getMZ(), exp2.getSpectra()[i][k].getMZ())
    }
  }
}

void cmpDataRT(MSExperiment& exp1, MSExperiment& exp2, double abs_tol = 1e-5, double rel_tol = 1+1e-5)
{
  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.

  TOLERANCE_ABSOLUTE(abs_tol)
  TOLERANCE_RELATIVE(rel_tol)
  for (Size i = 0; i < exp1.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp1.getChromatogram(i).size() == exp2.getChromatograms()[i].size(), true)
    for (Size k = 0; k < exp1.getChromatogram(i).size(); k++)
    {
      TEST_REAL_SIMILAR(exp1.getChromatogram(i)[k].getRT(), exp2.getChromatograms()[i][k].getRT())
    }
  }
}

///////////////////////////

START_TEST(MzMLSqliteHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzMLSqliteHandler* ptr = nullptr;
MzMLSqliteHandler* nullPointer = nullptr;
START_SECTION((MzMLSqliteHandler()))
  ptr = new MzMLSqliteHandler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MzMLSqliteHandler()))
  delete ptr;
END_SECTION

TOLERANCE_RELATIVE(1.0005)

START_SECTION(void readExperiment(MSExperiment & exp, bool meta_only = false) const)
{
  MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));

  MSExperiment exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp2);

  // read in meta data only
  {
    MSExperiment exp;
    handler.readExperiment(exp, true);
    TEST_EQUAL(exp.getNrSpectra(), exp2.getSpectra().size())
    TEST_EQUAL(exp.getNrChromatograms(), exp2.getChromatograms().size())
    TEST_EQUAL(exp.getNrSpectra(), 2)
    TEST_EQUAL(exp.getNrChromatograms(), 1)
    TEST_EQUAL(exp.getSpectrum(0) == exp2.getSpectra()[0], false) // no exact duplicate

    for (Size i = 0; i < exp.getNrSpectra(); i++)
    {
      TEST_EQUAL(exp.getSpectrum(i).size(), 0)
    }

    for (Size i = 0; i < exp.getNrChromatograms(); i++)
    {
      TEST_EQUAL(exp.getChromatogram(i).size(), 0)
    }
    TEST_EQUAL(exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, true)
  }

  MSExperiment exp;
  handler.readExperiment(exp, false);

  TEST_EQUAL(exp.getNrSpectra(), exp2.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp2.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp2.getSpectra()[0], false) // no exact duplicate

  cout.precision(17);

  cmpDataIntensity(exp, exp, 1e-4, 1.001); 
  cmpDataMZ(exp, exp, 1e-5, 1.000001); // less than 1ppm error for m/z 
  cmpDataRT(exp, exp, 0.05, 1.000001); // max 0.05 seconds error in RT

  // 1:1 mapping of experimental settings ...
  TEST_EQUAL(exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, true)
}
END_SECTION

START_SECTION( Size getNrSpectra() const )
{
  MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));
  TEST_EQUAL(handler.getNrSpectra(), 2)
}
END_SECTION

START_SECTION( Size getNrChromatograms() const )
{
  MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));
  TEST_EQUAL(handler.getNrChromatograms(), 1)
}
END_SECTION

START_SECTION( void readSpectra(std::vector<MSSpectrum> & exp, const std::vector<int> & indices, bool meta_only = false) const)
{
  MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));

  MSExperiment exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), exp2);

  // read in meta data only
  {
    std::vector<MSSpectrum> exp;
    std::vector<int> indices = {1};
    handler.readSpectra(exp, indices, true);
    TEST_EQUAL(exp.size(), 1)
    TEST_EQUAL(exp[0].size(), 0)
    TEST_REAL_SIMILAR(exp[0].getRT(), 0.4738)
  }

  {
    std::vector<MSSpectrum> exp;
    std::vector<int> indices = {1};
    handler.readSpectra(exp, indices, false);
    TEST_EQUAL(exp.size(), 1)
    TEST_EQUAL(exp[0].size(), 19800)
    TEST_REAL_SIMILAR(exp[0].getRT(), 0.4738)
  }

  {
    std::vector<MSSpectrum> exp;
    std::vector<int> indices = {0};
    handler.readSpectra(exp, indices, false);
    TEST_EQUAL(exp.size(), 1)
    TEST_EQUAL(exp[0].size(), 19914)
    TEST_REAL_SIMILAR(exp[0].getRT(), 0.2961)
  }

  {
    std::vector<MSSpectrum> exp;
    std::vector<int> indices = {0, 1};
    handler.readSpectra(exp, indices, true);
    TEST_EQUAL(exp.size(), 2)
    TEST_EQUAL(exp[0].size(), 0)
    TEST_EQUAL(exp[1].size(), 0)
    TEST_REAL_SIMILAR(exp[0].getRT(), 0.2961)
    TEST_REAL_SIMILAR(exp[1].getRT(), 0.4738)
  }

  {
    std::vector<MSSpectrum> exp;
    std::vector<int> indices = {0, 1};
    handler.readSpectra(exp, indices, false);
    TEST_EQUAL(exp.size(), 2)
    TEST_EQUAL(exp[0].size(), 19914)
    TEST_EQUAL(exp[1].size(), 19800)
    TEST_REAL_SIMILAR(exp[0].getRT(), 0.2961)
    TEST_REAL_SIMILAR(exp[1].getRT(), 0.4738)
  }

  {
    std::vector<MSSpectrum> exp;
    std::vector<int> indices = {0, 1, 2};
    TEST_EXCEPTION(Exception::IllegalArgument, handler.readSpectra(exp, indices, false));
  }

  {
    std::vector<MSSpectrum> exp;
    std::vector<int> indices = {5};
    TEST_EXCEPTION(Exception::IllegalArgument, handler.readSpectra(exp, indices, false));
  }
}
END_SECTION

START_SECTION(void readChromatograms(std::vector<MSChromatogram> & exp, const std::vector<int> & indices, bool meta_only = false) const)
{
  MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));

  MSExperiment exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), exp2);

  // read in meta data only
  {
    std::vector<MSChromatogram> exp;
    std::vector<int> indices = {0};
    handler.readChromatograms(exp, indices, true);
    TEST_EQUAL(exp.size(), 1)
    TEST_EQUAL(exp[0].size(), 0)
    TEST_STRING_EQUAL(exp[0].getNativeID(), "TIC")
  }

  {
    std::vector<MSChromatogram> exp;
    std::vector<int> indices = {0, 1};
    TEST_EXCEPTION(Exception::IllegalArgument, handler.readChromatograms(exp, indices, false));
  }

  {
    std::vector<MSChromatogram> exp;
    std::vector<int> indices = {5};
    TEST_EXCEPTION(Exception::IllegalArgument, handler.readChromatograms(exp, indices, false));
  }

  {

    MSExperiment exp_orig;
    MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), exp_orig);

    std::string tmp_filename;
    NEW_TMP_FILE(tmp_filename);

    // delete file if present
    QFile file (String(tmp_filename).toQString());
    file.remove();

    auto chroms = exp_orig.getChromatograms();
    chroms.push_back(exp_orig.getChromatograms()[0]);
    chroms.back().setNativeID("second");

    {
      MzMLSqliteHandler handler(tmp_filename);
      handler.setConfig(true, false, 0.0001);
      handler.createTables();
      handler.writeChromatograms(chroms);
    }

    MzMLSqliteHandler handler(tmp_filename);
    {
      std::vector<MSChromatogram> exp;
      std::vector<int> indices = {0};
      handler.readChromatograms(exp, indices, true);
      TEST_EQUAL(exp.size(), 1)
      TEST_EQUAL(exp[0].size(), 0)
      TEST_STRING_EQUAL(exp[0].getNativeID(), "TIC")
    }

    {
      std::vector<MSChromatogram> exp;
      std::vector<int> indices = {1};
      handler.readChromatograms(exp, indices, true);
      TEST_EQUAL(exp.size(), 1)
      TEST_EQUAL(exp[0].size(), 0)
      TEST_STRING_EQUAL(exp[0].getNativeID(), "second")
    }

    {
      std::vector<MSChromatogram> exp;
      std::vector<int> indices = {0, 1};
      handler.readChromatograms(exp, indices, true);
      TEST_EQUAL(exp.size(), 2)
      TEST_EQUAL(exp[0].size(), 0)
      TEST_STRING_EQUAL(exp[0].getNativeID(), "TIC")
      TEST_STRING_EQUAL(exp[1].getNativeID(), "second")
    }

  }

}
END_SECTION

START_SECTION(std::vector<size_t> getSpectraIndicesbyRT(double RT, double deltaRT, const std::vector<int> & indices) const)
{
  MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));

  {
    std::vector<int> indices = {};
    auto res = handler.getSpectraIndicesbyRT(0.4738, 0.1, indices);
    TEST_EQUAL(res.size(), 1)
    TEST_EQUAL(res[0], 1)
  }

  {
    std::vector<int> indices = {};
    auto res = handler.getSpectraIndicesbyRT(0.296, 0.1, indices);
    TEST_EQUAL(res.size(), 1)
    TEST_EQUAL(res[0], 0)
  }

  {
    std::vector<int> indices = {};
    auto res = handler.getSpectraIndicesbyRT(0.296, 1.1, indices);
    TEST_EQUAL(res.size(), 2)
    TEST_EQUAL(res[0], 0)
    TEST_EQUAL(res[1], 1)
  }

  {
    std::vector<int> indices = {1};
    auto res = handler.getSpectraIndicesbyRT(0.296, 1.1, indices);
    TEST_EQUAL(res.size(), 1)
    TEST_EQUAL(res[0], 1)
  }

  {
    std::vector<int> indices = {0};
    auto res = handler.getSpectraIndicesbyRT(0.296, 1.1, indices);
    TEST_EQUAL(res.size(), 1)
    TEST_EQUAL(res[0], 0)
  }

  {
    std::vector<int> indices = {};
    auto res = handler.getSpectraIndicesbyRT(0.0, 0.1, indices);
    TEST_EQUAL(res.size(), 0)
  }

  // negative deltaRT will simply return the first spectrum
  {
    std::vector<int> indices = {};
    auto res = handler.getSpectraIndicesbyRT(0.3, -0.1, indices);
    TEST_EQUAL(res.size(), 1)
    TEST_EQUAL(res[0], 1)
  }

  {
    std::vector<int> indices = {};
    auto res = handler.getSpectraIndicesbyRT(0.0, -0.1, indices);
    TEST_EQUAL(res.size(), 1)
    TEST_EQUAL(res[0], 0)
  }


}
END_SECTION

START_SECTION(void writeExperiment(const MSExperiment & exp))
{
  MSExperiment exp_orig;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), exp_orig);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  // delete file if present
  QFile file (String(tmp_filename).toQString());
  file.remove();

  {
    MzMLSqliteHandler handler(tmp_filename);
    // writing without creating the tables / indices won't work
    TEST_EXCEPTION(Exception::IllegalArgument, handler.writeExperiment(exp_orig));

    // now it will work
    handler.createTables();
    handler.createTables();
    handler.writeExperiment(exp_orig);

    // you can createTables() twice, but it will delete all your data 
    TEST_EQUAL(handler.getNrSpectra(), 2)
    handler.createTables();
    TEST_EQUAL(handler.getNrSpectra(), 0)
    handler.writeExperiment(exp_orig);
    TEST_EQUAL(handler.getNrSpectra(), 2)
  }

  MzMLSqliteHandler handler(tmp_filename);
  MSExperiment exp2 = exp_orig;

  // read in meta data only
  {
    MSExperiment exp;
    handler.readExperiment(exp, true);
    TEST_EQUAL(exp.getNrSpectra(), exp2.getSpectra().size())
    TEST_EQUAL(exp.getNrChromatograms(), exp2.getChromatograms().size())
    TEST_EQUAL(exp.getNrSpectra(), 2)
    TEST_EQUAL(exp.getNrChromatograms(), 1)
    TEST_EQUAL(exp.getSpectrum(0) == exp2.getSpectra()[0], false) // no exact duplicate

    for (Size i = 0; i < exp.getNrSpectra(); i++)
    {
      TEST_EQUAL(exp.getSpectrum(i).size(), 0)
    }

    for (Size i = 0; i < exp.getNrChromatograms(); i++)
    {
      TEST_EQUAL(exp.getChromatogram(i).size(), 0)
    }
    TEST_EQUAL(exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, true)
  }

  MSExperiment exp;
  handler.readExperiment(exp, false);

  TEST_EQUAL(exp.getNrSpectra(), exp2.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp2.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp2.getSpectra()[0], false) // no exact duplicate

  cmpDataIntensity(exp, exp, 1e-4, 1.001); 
  cmpDataMZ(exp, exp, 1e-5, 1.000001); // less than 1ppm error for m/z 
  cmpDataRT(exp, exp, 0.05, 1.000001); // max 0.05 seconds error in RT

  // 1:1 mapping of experimental settings ...
  TEST_EQUAL(exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, true)
}
END_SECTION

START_SECTION(void writeSpectra(const std::vector<MSSpectrum>& spectra))
{
  MSExperiment exp_orig;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), exp_orig);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  // delete file if present
  QFile file (String(tmp_filename).toQString());
  file.remove();

  {
    MzMLSqliteHandler handler(tmp_filename);
    // writing without creating the tables / indices won't work
    TEST_EXCEPTION(Exception::IllegalArgument, handler.writeSpectra(exp_orig.getSpectra()));

    // now it will work
    handler.createTables();
    handler.createTables();
    handler.writeSpectra(exp_orig.getSpectra());
    TEST_EQUAL(handler.getNrSpectra(), 2)
    handler.writeSpectra(exp_orig.getSpectra());
    TEST_EQUAL(handler.getNrSpectra(), 4)
    handler.writeSpectra(exp_orig.getSpectra());
    TEST_EQUAL(handler.getNrSpectra(), 6)
    MSExperiment tmp;
    handler.readExperiment(tmp, false);
    TEST_EQUAL(tmp.getNrSpectra(), 6)
    TEST_EQUAL(tmp[0].size(), 19914)
    TEST_EQUAL(tmp[1].size(), 19800)
    TEST_EQUAL(tmp[2].size(), 19914)
    TEST_EQUAL(tmp[3].size(), 19800)
    TEST_EQUAL(tmp[4].size(), 19914)
    TEST_EQUAL(tmp[5].size(), 19800)

    TEST_REAL_SIMILAR(tmp.getSpectra()[0][100].getMZ(), 204.817)
    TEST_REAL_SIMILAR(tmp.getSpectra()[0][100].getIntensity(), 3857.86)

    // clear
    handler.createTables();
    handler.writeSpectra(exp_orig.getSpectra());
    TEST_EQUAL(handler.getNrSpectra(), 2)
  }

}
END_SECTION

START_SECTION(void writeChromatograms(const std::vector<MSChromatogram>& chroms))
{
  MSExperiment exp_orig;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"), exp_orig);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  // delete file if present
  QFile file (String(tmp_filename).toQString());
  file.remove();

  {
    MzMLSqliteHandler handler(tmp_filename);
    handler.setConfig(true, false, 0.0001);
    // writing without creating the tables / indices won't work
    TEST_EXCEPTION(Exception::IllegalArgument, handler.writeChromatograms(exp_orig.getChromatograms()));

    // now it will work
    handler.createTables();
    handler.createTables();
    handler.writeChromatograms(exp_orig.getChromatograms());
    TEST_EQUAL(handler.getNrChromatograms(), 1)
    handler.writeChromatograms(exp_orig.getChromatograms());
    TEST_EQUAL(handler.getNrChromatograms(), 2)
    handler.writeChromatograms(exp_orig.getChromatograms());
    TEST_EQUAL(handler.getNrChromatograms(), 3)

    MSExperiment tmp;
    handler.readExperiment(tmp, false);
    TEST_EQUAL(tmp.getNrChromatograms(), 3)
    TEST_EQUAL(tmp.getChromatograms()[0].size(), 48)
    TEST_EQUAL(tmp.getChromatograms()[1].size(), 48)
    TEST_EQUAL(tmp.getChromatograms()[2].size(), 48)

    TEST_REAL_SIMILAR(tmp.getChromatograms()[0][20].getRT(), 0.200695)
    TEST_REAL_SIMILAR(tmp.getChromatograms()[0][20].getIntensity(), 147414.578125)

    // clear
    handler.createTables();
    handler.writeChromatograms(exp_orig.getChromatograms());
    TEST_EQUAL(handler.getNrChromatograms(), 1)
  }

  // now test with numpress (accuracy is lower)
  TOLERANCE_RELATIVE(1+2e-4)
  // delete file if present
  file.remove();
  {
    MzMLSqliteHandler handler(tmp_filename);
    handler.setConfig(true, true, 0.0001);
    // writing without creating the tables / indices won't work
    TEST_EXCEPTION(Exception::IllegalArgument, handler.writeChromatograms(exp_orig.getChromatograms()));

    // now it will work
    handler.createTables();
    handler.createTables();
    handler.writeChromatograms(exp_orig.getChromatograms());
    TEST_EQUAL(handler.getNrChromatograms(), 1)
    handler.writeChromatograms(exp_orig.getChromatograms());
    TEST_EQUAL(handler.getNrChromatograms(), 2)
    handler.writeChromatograms(exp_orig.getChromatograms());
    TEST_EQUAL(handler.getNrChromatograms(), 3)

    MSExperiment tmp;
    handler.readExperiment(tmp, false);
    TEST_EQUAL(tmp.getNrChromatograms(), 3)
    TEST_EQUAL(tmp.getChromatograms()[0].size(), 48)
    TEST_EQUAL(tmp.getChromatograms()[1].size(), 48)
    TEST_EQUAL(tmp.getChromatograms()[2].size(), 48)

    TEST_REAL_SIMILAR(tmp.getChromatograms()[0][20].getRT(), 0.200695)
    TEST_REAL_SIMILAR(tmp.getChromatograms()[0][20].getIntensity(), 147414.578125)

    // clear
    handler.createTables();
    handler.writeChromatograms(exp_orig.getChromatograms());
    TEST_EQUAL(handler.getNrChromatograms(), 1)
  }
}
END_SECTION

// reset error tolerances to default values
TOLERANCE_ABSOLUTE(1e-5)
TOLERANCE_RELATIVE(1+1e-5)

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

