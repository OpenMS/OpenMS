// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak2D.h>

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(double a, double b)
{
  DPosition<1> pa(a), pb(b);
  return DRange<1>(pa, pb);
}

///////////////////////////

START_TEST(DTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DTA2DFile* ptr = nullptr;
DTA2DFile* nullPointer = nullptr;
START_SECTION((DTA2DFile()))
  ptr = new DTA2DFile;
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~DTA2DFile()))
  delete ptr;
END_SECTION

START_SECTION(const PeakFileOptions& getOptions() const)
  DTA2DFile file;
  TEST_EQUAL(file.getOptions().hasMSLevels(),false)
END_SECTION

START_SECTION(PeakFileOptions& getOptions())
  DTA2DFile file;
  file.getOptions().addMSLevel(1);
  TEST_EQUAL(file.getOptions().hasMSLevels(),true);
END_SECTION

START_SECTION((template<typename MapType> void load(const String& filename, MapType& map) ))
  TOLERANCE_ABSOLUTE(0.01)

  PeakMap e;
  DTA2DFile file;

  //test exception
  TEST_EXCEPTION( Exception::FileNotFound , file.load("dummy/dummy.dta2d",e) )

  // real test
  file.load(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"),e);

  //test DocumentIdentifier addition
  TEST_STRING_EQUAL(e.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"));
  TEST_STRING_EQUAL(FileTypes::typeToName(e.getLoadedFileType()),"dta2d");

  TEST_EQUAL(e.size(), 9);
  ABORT_IF(e.size() != 9)

  TEST_STRING_EQUAL(e[0].getNativeID(),"index=0")
  TEST_STRING_EQUAL(e[1].getNativeID(),"index=1")
  TEST_STRING_EQUAL(e[2].getNativeID(),"index=2")
  TEST_STRING_EQUAL(e[3].getNativeID(),"index=3")
  TEST_STRING_EQUAL(e[4].getNativeID(),"index=4")
  TEST_STRING_EQUAL(e[5].getNativeID(),"index=5")
  TEST_STRING_EQUAL(e[6].getNativeID(),"index=6")
  TEST_STRING_EQUAL(e[7].getNativeID(),"index=7")
  TEST_STRING_EQUAL(e[8].getNativeID(),"index=8")

  PeakMap::const_iterator it(e.begin());
  TEST_REAL_SIMILAR((*it)[0].getPosition()[0], 230.02)
  TEST_REAL_SIMILAR(it->getRT(), 4711.1)
  TEST_REAL_SIMILAR((*it)[0].getIntensity(), 47218.89)
  ++it;

  TEST_REAL_SIMILAR((*it)[0].getPosition()[0], 231.51)
  TEST_REAL_SIMILAR(it->getRT(), 4711.2)
  TEST_REAL_SIMILAR((*it)[0].getIntensity(), 89935.22)
  ++it;

  TEST_REAL_SIMILAR((*it)[0].getPosition()[0], 139.42)
  TEST_REAL_SIMILAR(it->getRT(), 4711.3)
  TEST_REAL_SIMILAR((*it)[0].getIntensity(), 318.52)
  ++it;

  TEST_REAL_SIMILAR((*it)[0].getPosition()[0], 149.93)
  TEST_REAL_SIMILAR(it->getRT(), 4711.4)
  TEST_REAL_SIMILAR((*it)[0].getIntensity(), 61870.99)
  ++it;

  TEST_REAL_SIMILAR((*it)[0].getPosition()[0], 169.65)
  TEST_REAL_SIMILAR(it->getRT(), 4711.5)
  TEST_REAL_SIMILAR((*it)[0].getIntensity(), 62074.22)
  ++it;

  TEST_REAL_SIMILAR((*it)[0].getPosition()[0], 189.30)
  TEST_REAL_SIMILAR(it->getRT(), 4711.6)
  TEST_REAL_SIMILAR((*it)[0].getIntensity(), 53737.85)
  ++it;

  TEST_REAL_SIMILAR((*it)[0].getPosition()[0], 202.28)
  TEST_REAL_SIMILAR(it->getRT(), 4711.7)
  TEST_REAL_SIMILAR((*it)[0].getIntensity(), 49410.25)
  ++it;

  TEST_REAL_SIMILAR((*it)[0].getPosition()[0], 207.82)
  TEST_REAL_SIMILAR(it->getRT(), 4711.8)
  TEST_REAL_SIMILAR((*it)[0].getIntensity(), 17038.71)
  ++it;

  TEST_REAL_SIMILAR((*it)[0].getPosition()[0], 219.72)
  TEST_REAL_SIMILAR(it->getRT(), 4711.9)
  TEST_REAL_SIMILAR((*it)[0].getIntensity(), 73629.98)


  //test with header
  file.load(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_2.dta2d"),e);
  std::vector<Peak2D> array;
  e.get2DData(array);
  TEST_EQUAL(array.size(), 11);
  ABORT_IF(array.size() != 11)

  std::vector<Peak2D>::const_iterator it2 = array.begin();

  TEST_REAL_SIMILAR(it2->getMZ(), 230.02)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.1)
  TEST_REAL_SIMILAR(it2->getIntensity(), 47218.89)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 430.02)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.1)
  TEST_REAL_SIMILAR(it2->getIntensity(), 47219.89)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 630.02)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.1)
  TEST_REAL_SIMILAR(it2->getIntensity(), 47210.89)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 231.51)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.2)
  TEST_REAL_SIMILAR(it2->getIntensity(), 89935.22)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 139.42)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.3)
  TEST_REAL_SIMILAR(it2->getIntensity(), 318.52)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 149.93)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.4)
  TEST_REAL_SIMILAR(it2->getIntensity(), 61870.99)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 169.65)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.5)
  TEST_REAL_SIMILAR(it2->getIntensity(), 62074.22)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 189.30)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.6)
  TEST_REAL_SIMILAR(it2->getIntensity(), 53737.85)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 202.28)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.7)
  TEST_REAL_SIMILAR(it2->getIntensity(), 49410.25)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 207.82)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.8)
  TEST_REAL_SIMILAR(it2->getIntensity(), 17038.71)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 219.72)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.9)
  TEST_REAL_SIMILAR(it2->getIntensity(), 73629.98)


  PeakMap e3;
  file.load(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"),e3);
  TEST_EQUAL(e3.size(), 9);
  ABORT_IF(e3.size() != 9)

  PeakMap::const_iterator it3(e3.begin());
  TEST_EQUAL(it3->size(), 3);
  ABORT_IF(it3->size() != 3)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.1)
  TEST_REAL_SIMILAR((*it3)[0].getPosition()[0], 230.02)
  TEST_REAL_SIMILAR((*it3)[0].getIntensity(), 47218.89)
  TEST_REAL_SIMILAR((*it3)[1].getPosition()[0], 430.02)
  TEST_REAL_SIMILAR((*it3)[1].getIntensity(), 47219.89)
  TEST_REAL_SIMILAR((*it3)[2].getPosition()[0], 630.02)
  TEST_REAL_SIMILAR((*it3)[2].getIntensity(), 47210.89)
  ++it3;

  TEST_REAL_SIMILAR((*it3)[0].getPosition()[0], 231.51)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.2)
  TEST_REAL_SIMILAR((*it3)[0].getIntensity(), 89935.22)
  ++it3;

  TEST_REAL_SIMILAR((*it3)[0].getPosition()[0], 139.42)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.3)
  TEST_REAL_SIMILAR((*it3)[0].getIntensity(), 318.52)
  ++it3;

  TEST_REAL_SIMILAR((*it3)[0].getPosition()[0], 149.93)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.4)
  TEST_REAL_SIMILAR((*it3)[0].getIntensity(), 61870.99)
  ++it3;

  TEST_REAL_SIMILAR((*it3)[0].getPosition()[0], 169.65)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.5)
  TEST_REAL_SIMILAR((*it3)[0].getIntensity(), 62074.22)
  ++it3;

  TEST_REAL_SIMILAR((*it3)[0].getPosition()[0], 189.30)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.6)
  TEST_REAL_SIMILAR((*it3)[0].getIntensity(), 53737.85)
  ++it3;

  TEST_REAL_SIMILAR((*it3)[0].getPosition()[0], 202.28)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.7)
  TEST_REAL_SIMILAR((*it3)[0].getIntensity(), 49410.25)
  ++it3;

  TEST_REAL_SIMILAR((*it3)[0].getPosition()[0], 207.82)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.8)
  TEST_REAL_SIMILAR((*it3)[0].getIntensity(), 17038.71)
  ++it3;

  TEST_REAL_SIMILAR((*it3)[0].getPosition()[0], 219.72)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.9)
  TEST_REAL_SIMILAR((*it3)[0].getIntensity(), 73629.98)

  //test with header and minutes instead of seconds
  PeakMap e4;
  file.load(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_3.dta2d"),e4);
  TEST_EQUAL(e4.size(),9)
  TEST_REAL_SIMILAR(e4[0].getRT(), 282666)
  TEST_REAL_SIMILAR(e4[1].getRT(), 282672)
  TEST_REAL_SIMILAR(e4[2].getRT(), 282678)
  TEST_REAL_SIMILAR(e4[3].getRT(), 282684)
  TEST_REAL_SIMILAR(e4[4].getRT(), 282690)

END_SECTION

START_SECTION((template<typename MapType> void store(const String& filename, const MapType& map) const ))
  TOLERANCE_ABSOLUTE(0.1)
  std::string tmp_filename;
  PeakMap e;
  DTA2DFile f;

  NEW_TMP_FILE(tmp_filename);
  f.load(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"),e);
  f.store(tmp_filename,e);

  PeakMap e2;
  f.load(tmp_filename,e2);
  std::vector<Peak2D> array;
  e2.get2DData(array);
  TEST_EQUAL(array.size(), 11);
  ABORT_IF(array.size() != 11)

  std::vector<Peak2D>::const_iterator it2 = array.begin();

  TEST_REAL_SIMILAR(it2->getMZ(), 230.02)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.1)
  TEST_REAL_SIMILAR(it2->getIntensity(), 47218.89)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 430.02)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.1)
  TEST_REAL_SIMILAR(it2->getIntensity(), 47219.89)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 630.02)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.1)
  TEST_REAL_SIMILAR(it2->getIntensity(), 47210.89)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 231.51)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.2)
  TEST_REAL_SIMILAR(it2->getIntensity(), 89935.22)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 139.42)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.3)
  TEST_REAL_SIMILAR(it2->getIntensity(), 318.52)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 149.93)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.4)
  TEST_REAL_SIMILAR(it2->getIntensity(), 61870.99)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 169.65)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.5)
  TEST_REAL_SIMILAR(it2->getIntensity(), 62074.22)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 189.30)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.6)
  TEST_REAL_SIMILAR(it2->getIntensity(), 53737.85)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 202.28)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.7)
  TEST_REAL_SIMILAR(it2->getIntensity(), 49410.25)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 207.82)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.8)
  TEST_REAL_SIMILAR(it2->getIntensity(), 17038.71)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 219.72)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.9)
  TEST_REAL_SIMILAR(it2->getIntensity(), 73629.98)

  PeakMap e3;
  f.load(tmp_filename,e3);
  std::vector<Peak2D > array2;
  e2.get2DData(array2);
  TEST_EQUAL(array2.size(), 11);
  ABORT_IF(array2.size() != 11)

  std::vector<Peak2D >::const_iterator it3 = array2.begin();

  TEST_REAL_SIMILAR(it3->getMZ(), 230.02)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.1)
  TEST_REAL_SIMILAR(it3->getIntensity(), 47218.89)
  ++it3;

  TEST_REAL_SIMILAR(it3->getMZ(), 430.02)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.1)
  TEST_REAL_SIMILAR(it3->getIntensity(), 47219.89)
  ++it3;

  TEST_REAL_SIMILAR(it3->getMZ(), 630.02)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.1)
  TEST_REAL_SIMILAR(it3->getIntensity(), 47210.89)
  ++it3;

  TEST_REAL_SIMILAR(it3->getMZ(), 231.51)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.2)
  TEST_REAL_SIMILAR(it3->getIntensity(), 89935.22)
  ++it3;

  TEST_REAL_SIMILAR(it3->getMZ(), 139.42)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.3)
  TEST_REAL_SIMILAR(it3->getIntensity(), 318.52)
  ++it3;

  TEST_REAL_SIMILAR(it3->getMZ(), 149.93)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.4)
  TEST_REAL_SIMILAR(it3->getIntensity(), 61870.99)
  ++it3;

  TEST_REAL_SIMILAR(it3->getMZ(), 169.65)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.5)
  TEST_REAL_SIMILAR(it3->getIntensity(), 62074.22)
  ++it3;

  TEST_REAL_SIMILAR(it3->getMZ(), 189.30)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.6)
  TEST_REAL_SIMILAR(it3->getIntensity(), 53737.85)
  ++it3;

  TEST_REAL_SIMILAR(it3->getMZ(), 202.28)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.7)
  TEST_REAL_SIMILAR(it3->getIntensity(), 49410.25)
  ++it3;

  TEST_REAL_SIMILAR(it3->getMZ(), 207.82)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.8)
  TEST_REAL_SIMILAR(it3->getIntensity(), 17038.71)
  ++it3;

  TEST_REAL_SIMILAR(it3->getMZ(), 219.72)
  TEST_REAL_SIMILAR(it3->getRT(), 4711.9)
  TEST_REAL_SIMILAR(it3->getIntensity(), 73629.98)


END_SECTION

START_SECTION((template<typename MapType> void storeTIC(const String& filename, const MapType& map) const ))
  TOLERANCE_ABSOLUTE(0.1)
  std::string tmp_filename;
  PeakMap e;
  DTA2DFile f;

  NEW_TMP_FILE(tmp_filename);
  f.load(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"),e);
  f.storeTIC(tmp_filename,e);

  PeakMap e2;
  f.load(tmp_filename,e2);
  std::vector<Peak2D> array;
  e2.get2DData(array);
  TEST_EQUAL(array.size(), 9);
  ABORT_IF(array.size() != 9)

  std::vector<Peak2D>::const_iterator it2 = array.begin();

  TEST_REAL_SIMILAR(it2->getMZ(), 0)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.1)
  TEST_REAL_SIMILAR(it2->getIntensity(), 141650)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 0)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.2)
  TEST_REAL_SIMILAR(it2->getIntensity(), 89935.22)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 0)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.3)
  TEST_REAL_SIMILAR(it2->getIntensity(), 318.52)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 0)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.4)
  TEST_REAL_SIMILAR(it2->getIntensity(), 61870.99)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 0)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.5)
  TEST_REAL_SIMILAR(it2->getIntensity(), 62074.22)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 0)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.6)
  TEST_REAL_SIMILAR(it2->getIntensity(), 53737.85)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 0)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.7)
  TEST_REAL_SIMILAR(it2->getIntensity(), 49410.25)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 0)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.8)
  TEST_REAL_SIMILAR(it2->getIntensity(), 17038.71)
  ++it2;

  TEST_REAL_SIMILAR(it2->getMZ(), 0)
  TEST_REAL_SIMILAR(it2->getRT(), 4711.9)
  TEST_REAL_SIMILAR(it2->getIntensity(), 73629.98)
END_SECTION

START_SECTION(([EXTRA] load with RT range))
  TOLERANCE_ABSOLUTE(0.01)

  PeakMap e;
  DTA2DFile file;

  file.getOptions().setRTRange(makeRange(4711.15, 4711.45));
  file.load(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"),e);

  TEST_EQUAL(e.size(), 3)

  TEST_REAL_SIMILAR(e[0].getRT(), 4711.2)
  TEST_EQUAL(e[0].size(), 1)
  TEST_REAL_SIMILAR(e[0][0].getMZ(), 231.51)
  TEST_STRING_EQUAL(e[0].getNativeID(),"index=1")

  TEST_REAL_SIMILAR(e[1].getRT(),  4711.3)
  TEST_EQUAL(e[1].size(), 1)
  TEST_REAL_SIMILAR(e[1][0].getMZ(), 139.42)
  TEST_STRING_EQUAL(e[1].getNativeID(),"index=2")

  TEST_REAL_SIMILAR(e[2].getRT(),  4711.4)
  TEST_EQUAL(e[2].size(), 1)
  TEST_REAL_SIMILAR(e[2][0].getMZ(), 149.93)
  TEST_STRING_EQUAL(e[2].getNativeID(),"index=3")

END_SECTION

START_SECTION(([EXTRA] load with MZ range))
  TOLERANCE_ABSOLUTE(0.01)

  PeakMap e;
  DTA2DFile file;

  file.getOptions().setMZRange(makeRange(150, 220));
  file.load(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"),e);

  TEST_EQUAL(e.size(), 5)

  TEST_REAL_SIMILAR(e[0].getRT(), 4711.5)
  TEST_EQUAL(e[0].size(), 1)
  TEST_REAL_SIMILAR(e[0][0].getMZ(), 169.65)
  TEST_STRING_EQUAL(e[0].getNativeID(),"index=4")

  TEST_REAL_SIMILAR(e[1].getRT(), 4711.6)
  TEST_EQUAL(e[1].size(), 1)
  TEST_REAL_SIMILAR(e[1][0].getMZ(), 189.30)
  TEST_STRING_EQUAL(e[1].getNativeID(),"index=5")

  TEST_REAL_SIMILAR(e[2].getRT(), 4711.7)
  TEST_EQUAL(e[2].size(), 1)
  TEST_REAL_SIMILAR(e[2][0].getMZ(), 202.28)
  TEST_STRING_EQUAL(e[2].getNativeID(),"index=6")

  TEST_REAL_SIMILAR(e[3].getRT(), 4711.8)
  TEST_EQUAL(e[3].size(), 1)
  TEST_REAL_SIMILAR(e[3][0].getMZ(), 207.82)
  TEST_STRING_EQUAL(e[3].getNativeID(),"index=7")

  TEST_REAL_SIMILAR(e[4].getRT(), 4711.9)
  TEST_EQUAL(e[4].size(), 1)
  TEST_REAL_SIMILAR(e[4][0].getMZ(), 219.72)
  TEST_STRING_EQUAL(e[4].getNativeID(),"index=8")

END_SECTION

START_SECTION(([EXTRA] load with intensity range))
  TOLERANCE_ABSOLUTE(0.01)

  PeakMap e;
  DTA2DFile file;

  file.getOptions().setIntensityRange(makeRange(30000, 70000));
  file.load(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"),e);

  TEST_EQUAL(e.size(), 5)

  TEST_REAL_SIMILAR(e[0].getRT(), 4711.1)
  TEST_EQUAL(e[0].size(), 3)
  TEST_REAL_SIMILAR(e[0][0].getMZ(), 230.02)
  TEST_REAL_SIMILAR(e[0][1].getMZ(), 430.02)
  TEST_REAL_SIMILAR(e[0][2].getMZ(), 630.02)
  TEST_STRING_EQUAL(e[0].getNativeID(),"index=0")

  TEST_REAL_SIMILAR(e[1].getRT(), 4711.4)
  TEST_EQUAL(e[1].size(), 1)
  TEST_REAL_SIMILAR(e[1][0].getMZ(), 149.93)
  TEST_STRING_EQUAL(e[1].getNativeID(),"index=3")

  TEST_REAL_SIMILAR(e[2].getRT(), 4711.5)
  TEST_EQUAL(e[2].size(), 1)
  TEST_REAL_SIMILAR(e[2][0].getMZ(), 169.65)
  TEST_STRING_EQUAL(e[2].getNativeID(),"index=4")

  TEST_REAL_SIMILAR(e[3].getRT(), 4711.6)
  TEST_EQUAL(e[3].size(), 1)
  TEST_REAL_SIMILAR(e[3][0].getMZ(), 189.30)
  TEST_STRING_EQUAL(e[3].getNativeID(),"index=5")

  TEST_REAL_SIMILAR(e[4].getRT(), 4711.7)
  TEST_EQUAL(e[4].size(), 1)
  TEST_REAL_SIMILAR(e[4][0].getMZ(), 202.28)
  TEST_STRING_EQUAL(e[4].getNativeID(),"index=6")

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
