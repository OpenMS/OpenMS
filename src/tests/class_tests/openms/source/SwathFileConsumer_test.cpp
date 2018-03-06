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

///////////////////////////

#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>

///////////////////////////

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>

using namespace OpenMS;

static bool SortSwathMapByLower(const OpenSwath::SwathMap left, const OpenSwath::SwathMap right)
{
  if (left.ms1 && right.ms1) return left.center < right.center;
  else if (left.ms1) return true;
  else if (right.ms1) return false;
  return left.lower < right.lower;
}

void getSwathFile(PeakMap& exp, int nr_swathes=32, bool ms1=true)
{
  if (ms1)
  {
    MSSpectrum s;
    s.setMSLevel(1);
    Peak1D p; p.setMZ(100); p.setIntensity(200);
    s.push_back(p);
    exp.addSpectrum(s);
  }
  for (int i = 0; i< nr_swathes; i++)
  {
    MSSpectrum s;
    s.setMSLevel(2);
    std::vector<Precursor> prec(1);
    prec[0].setIsolationWindowLowerOffset(12.5);
    prec[0].setIsolationWindowUpperOffset(12.5);
    prec[0].setMZ(400 + i*25 + 12.5);
    s.setPrecursors(prec);
    Peak1D p; p.setMZ(101 + i); p.setIntensity(201 + i);
    s.push_back(p);
    exp.addSpectrum(s);
  }
}

START_TEST(SwathFileConsumer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Test "regular" / in memory consumer
{

RegularSwathFileConsumer* regular_sfc_ptr = nullptr;
RegularSwathFileConsumer* regular_sfc_nullPointer = nullptr;

START_SECTION(([EXTRA] RegularSwathFileConsumer()))
  regular_sfc_ptr = new RegularSwathFileConsumer;
  TEST_NOT_EQUAL(regular_sfc_ptr, regular_sfc_nullPointer)
END_SECTION

START_SECTION(([EXTRA] virtual ~RegularSwathFileConsumer()))
    delete regular_sfc_ptr;
END_SECTION

START_SECTION(([EXTRA] consumeAndRetrieve))
{
  regular_sfc_ptr = new RegularSwathFileConsumer();
  PeakMap exp;
  getSwathFile(exp);
  // Consume all the spectra
  for (Size i = 0; i < exp.getSpectra().size(); i++)
  {
    regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]);
  }

  std::vector< OpenSwath::SwathMap > maps;
  regular_sfc_ptr->retrieveSwathMaps(maps);

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
  delete regular_sfc_ptr;
}
END_SECTION

START_SECTION(([EXTRA] consumeAndRetrieve_known_boundaries))
{
  // Using the second constructor
  std::vector< OpenSwath::SwathMap > boundaries;
  PeakMap exp;
  getSwathFile(exp);
  for (int i = 0; i< 32; i++)
  {
    OpenSwath::SwathMap m;
    m.center = 400 + i*25 + 12.5;
    // enforce slightly different windows than the one in the file
    m.lower = m.center - 5;
    m.upper = m.center + 5;
    boundaries.push_back(m);
  }
  regular_sfc_ptr = new RegularSwathFileConsumer(boundaries);

  // Consume all the spectra
  for (Size i = 0; i < exp.getSpectra().size(); i++)
  {
    regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]);
  }

  std::vector< OpenSwath::SwathMap > maps;
  regular_sfc_ptr->retrieveSwathMaps(maps);

  TEST_EQUAL(maps.size(), 33)
  TEST_EQUAL(maps[0].ms1, true)
  for (Size i = 0; i< 32; i++)
  {
    TEST_EQUAL(maps[i+1].ms1, false)
    TEST_EQUAL(maps[i+1].sptr->getNrSpectra(), 1)
    TEST_EQUAL(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i+1].lower, 400+i*25.0 + 12.5 - 5.0)
    TEST_REAL_SIMILAR(maps[i+1].upper, 400+i*25.0 + 12.5 + 5.0)
  }
  delete regular_sfc_ptr;
}
END_SECTION

START_SECTION(([EXTRA] consumeAndRetrieve_scrambled))
{
  // Feed the SWATH maps to the consumer in a scrambled fashion:
  // Consume an MS2 spectrum, then an MS1 spectrum, then 5 more MS2, then an MS1, skip 5
  // of the MS2 spectra and consumer another 5 MS2 spectra.

  regular_sfc_ptr = new RegularSwathFileConsumer();
  PeakMap exp;
  getSwathFile(exp);

  regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[6]);   // MS2
  regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[0]);   // MS1
  for (Size i = 1; i < 5; i++)
  {
    regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]); // MS2
  }
  regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[0]);   // MS1
  for (Size i = 10; i < 15; i++)
  {
    regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]); // MS2
  }

  // Consume all spectra again to make sure we have "seen" them all (1 MS1 + 32 MS2)
  for (Size i = 0; i < exp.getSpectra().size(); i++)
  {
    regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]);
  }

  std::vector< OpenSwath::SwathMap > maps;
  regular_sfc_ptr->retrieveSwathMaps(maps);

  TEST_EQUAL(maps.size(), 33)
  std::sort(maps.begin(), maps.end(), SortSwathMapByLower);
  TEST_EQUAL(maps[0].ms1, true)
  TEST_EQUAL(maps[0].sptr->getNrSpectra(), 3)
  for (Size i = 0; i< 32; i++)
  {
    TEST_EQUAL(maps[i+1].ms1, false)
    if (i>15) TEST_EQUAL(maps[i+1].sptr->getNrSpectra(), 1) // some now also have 2 or 3 spectra
    TEST_EQUAL(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i+1].lower, 400+i*25.0)
    TEST_REAL_SIMILAR(maps[i+1].upper, 425+i*25.0)
  }
  delete regular_sfc_ptr;
}
END_SECTION

START_SECTION(([EXTRA] consumeAndRetrieve_scrambled_known_boundaries))
{
  // Using the second constructor
  std::vector< OpenSwath::SwathMap > boundaries;
  PeakMap exp;
  getSwathFile(exp);
  // add some extra windows for confusion
  for (int i = 0; i < 2; i++)
  {
    OpenSwath::SwathMap m;
    m.center = 100 + i*25 + 12.5;
    m.lower = m.center - 5;
    m.upper = m.center + 5;
    boundaries.push_back(m);
  }

  for (int i = 0; i < 32; i++)
  {
    OpenSwath::SwathMap m;
    m.center = 400 + i*25 + 12.5;
    // enforce slightly different windows than the one in the file
    m.lower = m.center - 5;
    m.upper = m.center + 5;
    boundaries.push_back(m);
  }

  // add some extra windows for confusion
  for (int i = 0; i < 2; i++)
  {
    OpenSwath::SwathMap m;
    m.center = 5000 + i*25 + 12.5;
    m.lower = m.center - 5;
    m.upper = m.center + 5;
    boundaries.push_back(m);
  }
  regular_sfc_ptr = new RegularSwathFileConsumer(boundaries);


  // Feed the SWATH maps to the consumer in a scrambled fashion:
  // Consume an MS2 spectrum, then an MS1 spectrum, then 5 more MS2, then an MS1, skip 5
  // of the MS2 spectra and consumer another 5 MS2 spectra.

  regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[6]);   // MS2
  regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[0]);   // MS1
  for (Size i = 1; i < 5; i++)
  {
    regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]); // MS2
  }
  regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[0]);   // MS1
  for (Size i = 10; i < 15; i++)
  {
    regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]); // MS2
  }

  // Consume all spectra again to make sure we have "seen" them all (1 MS1 + 32 MS2)
  for (Size i = 0; i < exp.getSpectra().size(); i++)
  {
    regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]);
  }

  std::vector< OpenSwath::SwathMap > maps;
  regular_sfc_ptr->retrieveSwathMaps(maps);

  TEST_EQUAL(maps.size(), 33+2)
  std::sort(maps.begin(), maps.end(), SortSwathMapByLower);
  TEST_EQUAL(maps[0].ms1, true)
  TEST_EQUAL(maps[0].sptr->getNrSpectra(), 3)

  // two empty ones
  TEST_EQUAL(maps[1].ms1, false)
  TEST_EQUAL(maps[1].sptr->getNrSpectra(), 0)
  TEST_EQUAL(maps[2].ms1, false)
  TEST_EQUAL(maps[2].sptr->getNrSpectra(), 0)

  for (Size i = 2; i< 32+2; i++)
  {
    TEST_EQUAL(maps[i+1].ms1, false)
    TEST_EQUAL(maps[i+1].sptr->getNrSpectra() > 0, true)
    TEST_EQUAL(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+(i-2))
    TEST_REAL_SIMILAR(maps[i+1].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+(i-2))
    TEST_REAL_SIMILAR(maps[i+1].lower, 400+(i-2)*25.0 + 12.5 - 5.0)
    TEST_REAL_SIMILAR(maps[i+1].upper, 400+(i-2)*25.0 + 12.5 + 5.0)
  }
  delete regular_sfc_ptr;
}
END_SECTION

START_SECTION(([EXTRA] consumeAndRetrieve_noMS1))
{
  regular_sfc_ptr = new RegularSwathFileConsumer();
  int nr_swath = 32;
  PeakMap exp;
  getSwathFile(exp, nr_swath, false);
  // Consume all the spectra
  for (Size i = 0; i < exp.getSpectra().size(); i++)
  {
    regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]);
  }

  std::vector< OpenSwath::SwathMap > maps;
  regular_sfc_ptr->retrieveSwathMaps(maps);

  TEST_EQUAL(maps.size(), nr_swath) // Swath number
  TEST_EQUAL(maps[0].ms1, false)
  for (int i = 0; i< nr_swath; i++)
  {
    TEST_EQUAL(maps[i].ms1, false)
    TEST_EQUAL(maps[i].sptr->getNrSpectra(), 1)
    TEST_EQUAL(maps[i].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i].lower, 400+i*25.0)
    TEST_REAL_SIMILAR(maps[i].upper, 425+i*25.0)
  }
  delete regular_sfc_ptr;
}
END_SECTION

START_SECTION(([EXTRA] consumeAndRetrieve_noMS2))
{
  int nr_swath = 0;
  regular_sfc_ptr = new RegularSwathFileConsumer();
  PeakMap exp;
  getSwathFile(exp, nr_swath, true);
  // Consume all the spectra
  for (Size i = 0; i < exp.getSpectra().size(); i++)
  {
    regular_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]);
  }

  std::vector< OpenSwath::SwathMap > maps;
  regular_sfc_ptr->retrieveSwathMaps(maps);

  TEST_EQUAL(maps.size(), 1) // Only MS1
  TEST_EQUAL(maps[0].ms1, true)
  TEST_EQUAL(maps[0].sptr->getNrSpectra(), 1)
  TEST_EQUAL(maps[0].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
  TEST_REAL_SIMILAR(maps[0].sptr->getSpectrumById(0)->getMZArray()->data[0], 100.0)
  TEST_REAL_SIMILAR(maps[0].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 200.0)
  delete regular_sfc_ptr;
}
END_SECTION

START_SECTION(([EXTRA] void retrieveSwathMaps(std::vector< OpenSwath::SwathMap > & maps))) 
{
  NOT_TESTABLE // already tested consumeAndRetrieve
}
END_SECTION

START_SECTION(([EXTRA] void consumeChromatogram(MapType::ChromatogramType &) )) 
{
  regular_sfc_ptr = new RegularSwathFileConsumer();
  MSChromatogram c;
  regular_sfc_ptr->consumeChromatogram(c);
  TEST_EQUAL(true, true)
  delete regular_sfc_ptr;
}
END_SECTION

START_SECTION(([EXTRA] void consumeSpectrum(MapType::SpectrumType & s))) 
{
  regular_sfc_ptr = new RegularSwathFileConsumer();
  MSSpectrum s;
  s.setMSLevel(1);
  regular_sfc_ptr->consumeSpectrum(s);

  s.setMSLevel(2);
  TEST_EXCEPTION(Exception::InvalidParameter, regular_sfc_ptr->consumeSpectrum(s))

  std::vector<Precursor> prec(1);
  s.setPrecursors(prec);
  TEST_EXCEPTION(Exception::InvalidParameter, regular_sfc_ptr->consumeSpectrum(s))

  prec[0].setIsolationWindowLowerOffset(12.5);
  prec[0].setIsolationWindowUpperOffset(12.5);
  s.setPrecursors(prec);
  TEST_EXCEPTION(Exception::InvalidParameter, regular_sfc_ptr->consumeSpectrum(s))

  prec[0].setMZ(100);
  s.setPrecursors(prec);
  regular_sfc_ptr->consumeSpectrum(s);
  delete regular_sfc_ptr;
}
END_SECTION
}

// Test cached consumer
// - shared functions in the base class are already tested, only test I/O here
{

CachedSwathFileConsumer* cached_sfc_ptr = nullptr;
CachedSwathFileConsumer* cached_sfc_nullPointer = nullptr;

START_SECTION(([EXTRA] CachedSwathFileConsumer()))
  cached_sfc_ptr = new CachedSwathFileConsumer("./", "tmp_osw_cached", 0, std::vector<int>());
  TEST_NOT_EQUAL(cached_sfc_ptr, cached_sfc_nullPointer)
END_SECTION

START_SECTION(([EXTRA] virtual ~CachedSwathFileConsumer()))
    delete cached_sfc_ptr;
END_SECTION

START_SECTION(([EXTRA] consumeAndRetrieve))
{
  // 2 SWATH should be sufficient for the test
  //int nr_swath = 1;
  int nr_swath = 2;
  std::vector<int> nr_ms2_spectra(nr_swath,1);
  cached_sfc_ptr = new CachedSwathFileConsumer("./", "tmp_osw_cached", 1, nr_ms2_spectra);
  PeakMap exp;
  getSwathFile(exp, nr_swath);
  // Consume all the spectra
  for (Size i = 0; i < exp.getSpectra().size(); i++)
  {
    cached_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]);
  }

  std::vector< OpenSwath::SwathMap > maps;
  cached_sfc_ptr->retrieveSwathMaps(maps);

  TEST_EQUAL(maps.size(), nr_swath+1) // Swath number + MS1 
  TEST_EQUAL(maps[0].ms1, true)
  TEST_EQUAL(maps[0].sptr->getNrSpectra(), 1)
  TEST_EQUAL(maps[0].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
  TEST_REAL_SIMILAR(maps[0].sptr->getSpectrumById(0)->getMZArray()->data[0], 100.0)
  TEST_REAL_SIMILAR(maps[0].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 200.0)

  for (int i = 0; i< nr_swath; i++)
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

START_SECTION(([EXTRA] consumeAndRetrieve_noMS1))
{
  // 2 SWATH should be sufficient for the test
  int nr_swath = 2;
  std::vector<int> nr_ms2_spectra(nr_swath,1);
  cached_sfc_ptr = new CachedSwathFileConsumer("./", "tmp_osw_cached", 1, nr_ms2_spectra);
  PeakMap exp;
  getSwathFile(exp, nr_swath, false);
  // Consume all the spectra
  for (Size i = 0; i < exp.getSpectra().size(); i++)
  {
    cached_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]);
  }

  std::vector< OpenSwath::SwathMap > maps;
  cached_sfc_ptr->retrieveSwathMaps(maps);

  TEST_EQUAL(maps.size(), nr_swath) // Swath number
  TEST_EQUAL(maps[0].ms1, false)
  for (int i = 0; i< nr_swath; i++)
  {
    TEST_EQUAL(maps[i].ms1, false)
    TEST_EQUAL(maps[i].sptr->getNrSpectra(), 1)
    TEST_EQUAL(maps[i].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
    TEST_REAL_SIMILAR(maps[i].sptr->getSpectrumById(0)->getMZArray()->data[0], 101.0+i)
    TEST_REAL_SIMILAR(maps[i].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 201.0+i)
    TEST_REAL_SIMILAR(maps[i].lower, 400+i*25.0)
    TEST_REAL_SIMILAR(maps[i].upper, 425+i*25.0)
  }

}
END_SECTION

START_SECTION(([EXTRA] consumeAndRetrieve_noMS2))
{
  int nr_swath = 0;
  std::vector<int> nr_ms2_spectra(nr_swath,1);
  cached_sfc_ptr = new CachedSwathFileConsumer("./", "tmp_osw_cached", 1, nr_ms2_spectra);
  PeakMap exp;
  getSwathFile(exp, nr_swath, true);
  // Consume all the spectra
  for (Size i = 0; i < exp.getSpectra().size(); i++)
  {
    cached_sfc_ptr->consumeSpectrum(exp.getSpectra()[i]);
  }

  std::vector< OpenSwath::SwathMap > maps;
  cached_sfc_ptr->retrieveSwathMaps(maps);

  TEST_EQUAL(maps.size(), 1) // Only MS1
  TEST_EQUAL(maps[0].ms1, true)
  TEST_EQUAL(maps[0].sptr->getNrSpectra(), 1)
  TEST_EQUAL(maps[0].sptr->getSpectrumById(0)->getMZArray()->data.size(), 1)
  TEST_REAL_SIMILAR(maps[0].sptr->getSpectrumById(0)->getMZArray()->data[0], 100.0)
  TEST_REAL_SIMILAR(maps[0].sptr->getSpectrumById(0)->getIntensityArray()->data[0], 200.0)
}
END_SECTION

START_SECTION(([EXTRA] void retrieveSwathMaps(std::vector< OpenSwath::SwathMap > & maps))) 
{
  NOT_TESTABLE // already tested consumeAndRetrieve
}
END_SECTION

START_SECTION(([EXTRA] void consumeChromatogram(MapType::ChromatogramType &) )) 
{
  NOT_TESTABLE // already tested consumeAndRetrieve
}
END_SECTION

START_SECTION(([EXTRA] void consumeSpectrum(MapType::SpectrumType & s))) 
{
  NOT_TESTABLE // already tested consumeAndRetrieve
}
END_SECTION
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
