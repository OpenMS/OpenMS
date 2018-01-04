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
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <boost/assign/std/vector.hpp>

#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerMRM.h>

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
      peak.setMZ(rtdata_1[k]);
      peak.setIntensity(intdata_1[k]);
    }
    else if (i == 1)
    {
      peak.setMZ(rtdata_2[k]);
      peak.setIntensity(intdata_2[k]);
    }
    chromatogram.push_back(peak);
  } 
  return chromatogram;
  
}

START_TEST(PeakPickerMRM, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerMRM* ptr = nullptr;
PeakPickerMRM* nullPointer = nullptr;

START_SECTION(PeakPickerMRM())
{
	ptr = new PeakPickerMRM();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~PeakPickerMRM())
{
  delete ptr;
}
END_SECTION

START_SECTION(void pickChromatogram(const RichPeakChromatogram &chromatogram, RichPeakChromatogram &picked_chrom))
{
  RichPeakChromatogram picked_chrom, smoothed_chrom, chrom;

  //RichPeakChromatogram chrom = transition_group.getChromatograms()[0];
  chrom = get_chrom(0);
  PeakPickerMRM picker;
  Param picker_param = picker.getDefaults();
  picker_param.setValue("method", "legacy"); // old parameters
  picker_param.setValue("peak_width", 40.0); // old parameters
  picker.setParameters(picker_param);
  picker.pickChromatogram(chrom, picked_chrom);

  TEST_EQUAL( picked_chrom.size(), 1);
  TEST_EQUAL( picked_chrom.getFloatDataArrays().size(), 3);

  // Peak picking is done by cubic spline interpolation and searching for the
  // point with zero derivative.
  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 9981.93933103869);
  TEST_REAL_SIMILAR( picked_chrom[0].getMZ(), 1495.11321013749);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[0][0], 60124.9); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[1][0], 1490.95); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[2][0], 1502.03); // rightWidth

  // chrom = transition_group.getChromatograms()[1];
  chrom = get_chrom(1);
  picker.pickChromatogram(chrom, picked_chrom);

  TEST_EQUAL( picked_chrom.size(), 1);
  TEST_EQUAL( picked_chrom.getFloatDataArrays().size(), 3);

  // Peak picking is done by cubic spline interpolation and searching for the
  // point with zero derivative.
  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 78719.134569503);
  TEST_REAL_SIMILAR( picked_chrom[0].getMZ(), 1492.830608593);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[0][0], 523378); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[1][0], 1481.84); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[2][0], 1501.23); // rightWidth

  ///////////////////////////////////////////////////////////////////////////
  // New method: Peak picking is done on the smoothed data and no minimal peak
  // width is set.
  chrom = get_chrom(0);
  picker_param.setValue("method", "corrected");
  picker_param.setValue("peak_width", -1.0);
  picker.setParameters(picker_param);
  picker.pickChromatogram(chrom, picked_chrom);
  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 9981.93933103869);
  TEST_REAL_SIMILAR( picked_chrom[0].getMZ(), 1495.11368082583);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[0][0], 60605.7); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[1][0], 1482.64); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[2][0], 1504.8);  // rightWidth

  chrom = get_chrom(1);
  picker.pickChromatogram(chrom, picked_chrom);
  TEST_EQUAL( picked_chrom.size(), 1);
  TEST_EQUAL( picked_chrom.getFloatDataArrays().size(), 3);

  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 78719.1346);
  TEST_REAL_SIMILAR( picked_chrom[0].getMZ(), 1492.8305);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[0][0], 525672.0); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[1][0], 1481.84);  // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[2][0], 1504.0);   // rightWidth

#ifdef WITH_CRAWDAD
  chrom = get_chrom(0);
  picker_param.setValue("method", "crawdad");
  picker_param.setValue("peak_width", 40.0); // old parameters
  picker.setParameters(picker_param);
  picker.pickChromatogram(chrom, picked_chrom);
  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 61366.56640625);
  TEST_REAL_SIMILAR( picked_chrom[0].getMZ(), 1496.48);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[0][0], 61366.6); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[1][0], 1479.88); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[2][0], 1510.33); // rightWidth

  chrom = get_chrom(1);
  picker.pickChromatogram(chrom, picked_chrom);
  TEST_EQUAL( picked_chrom.size(), 1);
  TEST_EQUAL( picked_chrom.getFloatDataArrays().size(), 3);

  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 533936.875);
  TEST_REAL_SIMILAR( picked_chrom[0].getMZ(), 1490.15);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[0][0], 533936.875); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[1][0], 1479.08); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[2][0], 1509.53); // rightWidth
#endif
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

