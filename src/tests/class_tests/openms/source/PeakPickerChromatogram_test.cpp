// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <boost/assign/std/vector.hpp>

#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerChromatogram.h>

using namespace OpenMS;
using namespace std;

typedef MSChromatogram RichPeakChromatogram;

RichPeakChromatogram get_chrom(int i)
{
  // this is a simulated SRM experiment where the two traces are not sampled at
  // the exact same time points, thus a resampling is necessary before applying
  // the algorithm.
  static const double rtdata_1[] = {1474.34, 1477.11,  1479.88, 1482.64, 1485.41, 1488.19, 1490.95, 1493.72, 1496.48, 1499.25, 1502.03, 1504.8 , 1507.56, 1510.33, 1513.09, 1515.87, 1518.64, 1521.42};
  static const double rtdata_2[] = {1473.55, 1476.31,  1479.08, 1481.84, 1484.61, 1487.39, 1490.15, 1492.92, 1495.69, 1498.45, 1501.23, 1504   , 1506.76, 1509.53, 1512.29, 1515.07, 1517.84, 1520.62};

  static const double intdata_1[] = {3.26958, 3.74189, 3.31075, 86.1901, 3.47528, 387.864, 13281  , 6375.84, 39852.6, 2.66726, 612.747, 3.34313, 793.12 , 3.29156, 4.00586, 4.1591 , 3.23035, 3.90591};
  static const double intdata_2[] = {3.44054, 2142.31, 3.58763, 3076.97, 6663.55, 45681 ,  157694 , 122844 , 86034.7, 85391.1, 15992.8, 2293.94, 6934.85, 2735.18, 459.413, 3.93863, 3.36564, 3.44005};

  RichPeakChromatogram chromatogram;
  for (int k = 0; k < 18; k++)
  {   
    ChromatogramPeak peak;
    if (i == 0)
    {
      peak.setRT(rtdata_1[k]);
      peak.setIntensity(intdata_1[k]);
    }
    else if (i == 1)
    {
      peak.setRT(rtdata_2[k]);
      peak.setIntensity(intdata_2[k]);
    }
    chromatogram.push_back(peak);
  } 
  return chromatogram;
  
}

START_TEST(PeakPickerChromatogram, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerChromatogram* ptr = nullptr;
PeakPickerChromatogram* nullPointer = nullptr;

START_SECTION(PeakPickerChromatogram())
{
	ptr = new PeakPickerChromatogram();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~PeakPickerChromatogram())
{
  delete ptr;
}
END_SECTION

START_SECTION(void pickChromatogram(const RichPeakChromatogram &chromatogram, RichPeakChromatogram &picked_chrom))
{
  RichPeakChromatogram picked_chrom, smoothed_chrom, chrom;

  //RichPeakChromatogram chrom = transition_group.getChromatograms()[0];
  chrom = get_chrom(0);
  PeakPickerChromatogram picker;
  Param picker_param = picker.getDefaults();
  picker_param.setValue("method", "legacy"); // old parameters
  picker_param.setValue("peak_width", 40.0); // old parameters
  picker.setParameters(picker_param);
  picker.pickChromatogram(chrom, picked_chrom);

  TEST_EQUAL( picked_chrom.size(), 1);
  TEST_EQUAL( picked_chrom.getFloatDataArrays().size(), PeakPickerChromatogram::SIZE_OF_FLOATINDICES);

  // Peak picking is done by cubic spline interpolation and searching for the
  // point with zero derivative.
  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 9981.93933103869);
  TEST_REAL_SIMILAR( picked_chrom[0].getRT(), 1495.11321013749);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE][0], 60124.9); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][0], 1490.95); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][0], 1502.03); // rightWidth

  // chrom = transition_group.getChromatograms()[1];
  chrom = get_chrom(1);
  picker.pickChromatogram(chrom, picked_chrom);

  TEST_EQUAL( picked_chrom.size(), 1);
  TEST_EQUAL( picked_chrom.getFloatDataArrays().size(), PeakPickerChromatogram::SIZE_OF_FLOATINDICES);

  // Peak picking is done by cubic spline interpolation and searching for the
  // point with zero derivative.
  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 78719.134569503);
  TEST_REAL_SIMILAR( picked_chrom[0].getRT(), 1492.830608593);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE][0], 523378); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][0], 1481.84); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][0], 1501.23); // rightWidth

  ///////////////////////////////////////////////////////////////////////////
  // New method: Peak picking is done on the smoothed data and no minimal peak
  // width is set.
  chrom = get_chrom(0);
  picker_param.setValue("method", "corrected");
  picker_param.setValue("peak_width", -1.0);
  picker.setParameters(picker_param);
  picker.pickChromatogram(chrom, picked_chrom);
  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 9981.93933103869);
  TEST_REAL_SIMILAR( picked_chrom[0].getRT(), 1495.11368082583);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE][0], 60605.7); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][0], 1482.64); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][0], 1504.8);  // rightWidth

  chrom = get_chrom(1);
  picker.pickChromatogram(chrom, picked_chrom);
  TEST_EQUAL( picked_chrom.size(), 1);
  TEST_EQUAL( picked_chrom.getFloatDataArrays().size(), PeakPickerChromatogram::SIZE_OF_FLOATINDICES);

  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 78719.1346);
  TEST_REAL_SIMILAR( picked_chrom[0].getRT(), 1492.8305);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE][0], 525672.0); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][0], 1481.84);  // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][0], 1504.0);   // rightWidth

#ifdef WITH_CRAWDAD
  chrom = get_chrom(0);
  picker_param.setValue("method", "crawdad");
  picker_param.setValue("peak_width", 40.0); // old parameters
  picker.setParameters(picker_param);
  picker.pickChromatogram(chrom, picked_chrom);
  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 61366.56640625);
  TEST_REAL_SIMILAR( picked_chrom[0].getRT(), 1496.48);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE][0], 61366.6); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][0], 1479.88); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][0], 1510.33); // rightWidth

  chrom = get_chrom(1);
  picker.pickChromatogram(chrom, picked_chrom);
  TEST_EQUAL( picked_chrom.size(), 1);
  TEST_EQUAL( picked_chrom.getFloatDataArrays().size(), 3);

  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 533936.875);
  TEST_REAL_SIMILAR( picked_chrom[0].getRT(), 1490.15);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE][0], 533936.875); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][0], 1479.08); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][0], 1509.53); // rightWidth
#endif
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

