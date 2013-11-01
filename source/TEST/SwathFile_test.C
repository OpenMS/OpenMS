// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

///////////////////////////

#include <OpenMS/FORMAT/SwathFile.h>

///////////////////////////

bool sortSwathMaps(const OpenSwath::SwathMap& left, const OpenSwath::SwathMap& right)
{
  // true if left is smaller
  if (left.ms1) return true;
  if (right.ms1) return false;
  return left.lower < right.lower;
}

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>

using namespace OpenMS;

void storeSwathFile(String filename, int nr_swathes=32)
{
  MSExperiment<> exp;
  {
    MSSpectrum<> s;
    s.setMSLevel(1);
    Peak1D p; p.setMZ(101); p.setIntensity(201);
    s.push_back(p);
    exp.addSpectrum(s);
  }
  for (int i = 0; i< nr_swathes; i++)
  {
    MSSpectrum<> s;
    s.setMSLevel(2);
    std::vector<Precursor> prec(1);
    prec[0].setIsolationWindowLowerOffset(400 + i*25);
    prec[0].setIsolationWindowUpperOffset(425 + i*25);
    prec[0].setMZ(400 + i/2.0 *25);
    s.setPrecursors(prec);
    Peak1D p; p.setMZ(101 + i); p.setIntensity(201 + i);
    s.push_back(p);
    exp.addSpectrum(s);
  }
  MzMLFile().store(filename, exp);
}

void storeSplitSwathFile(std::vector<String> filenames)
{
  {
    MSExperiment<> exp;
    MSSpectrum<> s;
    s.setMSLevel(1);
    Peak1D p; p.setMZ(101); p.setIntensity(201);
    s.push_back(p);
    exp.addSpectrum(s);
    MzMLFile().store(filenames[0], exp);
  }
  for (Size i = 0; i< filenames.size() -1; i++)
  {
    MSExperiment<> exp;
    MSSpectrum<> s;
    s.setMSLevel(2);
    std::vector<Precursor> prec(1);
    prec[0].setIsolationWindowLowerOffset(400 + i*25);
    prec[0].setIsolationWindowUpperOffset(425 + i*25);
    prec[0].setMZ(400 + i/2.0 *25);
    s.setPrecursors(prec);
    Peak1D p; p.setMZ(101 + i); p.setIntensity(201 + i);
    s.push_back(p);
    exp.addSpectrum(s);
    MzMLFile().store(filenames[i+1], exp);
  }
}

START_TEST(SwathFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SwathFile* swath_file_ptr = 0;
SwathFile* swath_file_nullPointer = 0;

START_SECTION((SwathFile()))
  swath_file_ptr = new SwathFile;
  TEST_NOT_EQUAL(swath_file_ptr, swath_file_nullPointer)
END_SECTION

START_SECTION(([EXTRA]virtual ~SwathFile()))
    delete swath_file_ptr;
END_SECTION

START_SECTION(std::vector< OpenSwath::SwathMap > loadMzML(String file, String tmp, 
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="normal") )
{
  storeSwathFile("swathFile_1.tmp");
  boost::shared_ptr<ExperimentalSettings> meta = boost::shared_ptr<ExperimentalSettings>(new ExperimentalSettings());
  std::vector< OpenSwath::SwathMap > maps = SwathFile().loadMzML("swathFile_1.tmp", "./", meta);

  TEST_EQUAL(maps.size(), 33)
  TEST_EQUAL(maps[0].ms1, true)
  for (Size i = 0; i< 32; i++)
  {
    TEST_EQUAL(maps[i+1].ms1, false)
    TEST_EQUAL(maps[i+1].sptr->getNrSpectra(), 1)
    TEST_EQUAL(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i+1].lower, 400+i*25.0)
    TEST_REAL_SIMILAR(maps[i+1].upper, 425+i*25.0)
  }
}
END_SECTION

START_SECTION([EXTRA]std::vector< OpenSwath::SwathMap > loadMzML(String file, String tmp, 
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="cache") )
{
  storeSwathFile("swathFile_1.tmp", 2);
  boost::shared_ptr<ExperimentalSettings> meta = boost::shared_ptr<ExperimentalSettings>(new ExperimentalSettings());
  std::vector< OpenSwath::SwathMap > maps = SwathFile().loadMzML("swathFile_1.tmp", "./", meta, "cache");

  TEST_EQUAL(maps.size(), 2+1)
  TEST_EQUAL(maps[0].ms1, true)
  for (Size i = 0; i< 2; i++)
  {
    TEST_EQUAL(maps[i+1].ms1, false)
    TEST_EQUAL(maps[i+1].sptr->getNrSpectra(), 1)
    TEST_EQUAL(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i+1].lower, 400+i*25.0)
    TEST_REAL_SIMILAR(maps[i+1].upper, 425+i*25.0)
  }
}
END_SECTION

START_SECTION(std::vector< OpenSwath::SwathMap > loadSplit(StringList file_list, String tmp, 
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="normal"))
{
  std::vector<String> swath_filenames;
  swath_filenames.push_back("swathFile_2_ms1.tmp");
  swath_filenames.push_back("swathFile_2_sw1.tmp");
  swath_filenames.push_back("swathFile_2_sw2.tmp");
  swath_filenames.push_back("swathFile_2_sw3.tmp");
  swath_filenames.push_back("swathFile_2_sw4.tmp");
  storeSplitSwathFile(swath_filenames);
  boost::shared_ptr<ExperimentalSettings> meta = boost::shared_ptr<ExperimentalSettings>(new ExperimentalSettings());
  std::vector< OpenSwath::SwathMap > maps = SwathFile().loadSplit(swath_filenames, "./", meta);

  // ensure they are sorted ... 
  std::sort(maps.begin(), maps.end(), sortSwathMaps);

  TEST_EQUAL(maps.size(), 5)
  TEST_EQUAL(maps[0].ms1, true)
  for (Size i = 0; i< maps.size() -1; i++)
  {
    TEST_EQUAL(maps[i+1].ms1, false)
    TEST_EQUAL(maps[i+1].sptr->getNrSpectra(), 1)
    TEST_EQUAL(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i+1].lower, 400+i*25.0)
    TEST_REAL_SIMILAR(maps[i+1].upper, 425+i*25.0)
  }

}
END_SECTION

START_SECTION([EXTRA]std::vector< OpenSwath::SwathMap > loadSplit(StringList file_list, String tmp, 
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="cache"))
{
  std::vector<String> swath_filenames;
  swath_filenames.push_back("swathFile_2_ms1.tmp");
  swath_filenames.push_back("swathFile_2_sw1.tmp");
  swath_filenames.push_back("swathFile_2_sw2.tmp");
  swath_filenames.push_back("swathFile_2_sw3.tmp");
  swath_filenames.push_back("swathFile_2_sw4.tmp");
  storeSplitSwathFile(swath_filenames);
  boost::shared_ptr<ExperimentalSettings> meta = boost::shared_ptr<ExperimentalSettings>(new ExperimentalSettings());
  std::vector< OpenSwath::SwathMap > maps = SwathFile().loadSplit(swath_filenames, "./", meta, "cache");

  // ensure they are sorted ... 
  std::sort(maps.begin(), maps.end(), sortSwathMaps);

  TEST_EQUAL(maps.size(), 5)
  TEST_EQUAL(maps[0].ms1, true)
  for (Size i = 0; i< maps.size() -1; i++)
  {
    TEST_EQUAL(maps[i+1].ms1, false)
    TEST_EQUAL(maps[i+1].sptr->getNrSpectra(), 1)
    TEST_EQUAL(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i+1].lower, 400+i*25.0)
    TEST_REAL_SIMILAR(maps[i+1].upper, 425+i*25.0)
  }

}
END_SECTION

START_SECTION(( std::vector< OpenSwath::SwathMap > loadMzXML(String file, String tmp, 
      boost::shared_ptr<ExperimentalSettings>& exp_meta, String readoptions="normal") ) )
{
  NOT_TESTABLE // mzXML is not supported
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
