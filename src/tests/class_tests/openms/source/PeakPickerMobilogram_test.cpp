// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <boost/assign/std/vector.hpp>

#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerMobilogram.h>

using namespace OpenMS;
using namespace std;

typedef Mobilogram RichPeakMobilogram;

RichPeakMobilogram get_mobilogram(int i)
{
  static const double im_data1[] = {0.93734956, 0.93835497, 0.9393905,  0.9404265,  0.9414631, 0.9425003,  0.9435083,  0.9445466,  0.9455854,  0.94662476, 0.947635,
       0.94867545, 0.9496867,  0.95072824, 0.95177037, 0.9527832, 0.9538264,  0.95484036, 0.9558847,  0.9568997,  0.9579451,  0.9589612,
       0.96000767, 0.9610248,  0.96207243, 0.96309066, 0.9641094, 0.96515864, 0.9661785,  0.96719885, 0.96824974, 0.9692712,  0.9702931,
       0.97131556, 0.9723687,  0.97339225, 0.9744164,  0.975441,  0.9764661,  0.977522,   0.9785482,  0.979575,   0.98060226, 0.98163015,
       0.9826585,  0.9836874,  0.98471683, 0.9857468,  0.9867773, 0.98780835, 0.98883986, 0.989872,   0.9909046,  0.9919378,  0.99297154,
       0.99397534, 0.99501014, 0.9960454,  0.9970813,  0.9981177, 1.0125301};

  static const double int_data1[] = {341.995757, 200.000741, 325.0076900000001, 79.002655, 227.00538, 487.00108, 576.0032550000001, 526.996605, 715.9864020000001, 778.010235, 686.003534, 548.0023570000001, 481.00395999999995, 547.00563, 950.006949, 1453.9987189999997, 831.026847, 1839.009966, 2306.021277, 2588.9953999999993, 2126.0273610000004, 3254.0281579999996, 3961.0108749999995, 4708.020719000002, 4991.072095, 5230.084562000001, 6559.103413000001, 6486.9137040000005, 9512.035332, 14461.755780000001, 16570.804681, 19887.687351000004, 22368.220470000004, 31683.985819, 37679.199531000006, 46487.67955, 48958.673168000016, 52820.41709300001, 58070.474128999995, 57655.34146700001, 58849.449743000005, 56964.34470800001, 56100.919934, 48964.276839, 43019.251419, 39263.106511000005, 28199.275315000003, 25018.205981999996, 20248.883416, 14727.949587000001, 10372.823519, 6309.068338999999, 5370.119231000001, 2253.037892, 2025.998686, 868.9924269999999, 209.0077, 44.99693, 317.999914, 57.00195, 78.0025};

  RichPeakMobilogram mobilogram;
  for (int k = 0; k < 60; k++)
  {
    MobilityPeak1D peak;
    if (i == 0)
    {
      peak.setPos(im_data1[k]);
      peak.setIntensity(int_data1[k]);
    }
    mobilogram.push_back(peak);
  }
  return mobilogram;
}

START_TEST(PeakPickerMobilogram, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerMobilogram* ptr = nullptr;
PeakPickerMobilogram* nullPointer = nullptr;

START_SECTION(PeakPickerMobilogram())
{
	ptr = new PeakPickerMobilogram();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~PeakPickerMobilogram())
{
  delete ptr;
}
END_SECTION

START_SECTION(void pickMobilogram(const RichPeakMobilogram& mobilogram, RichPeakMobilogram& picked_mobilogram))
{
  RichPeakMobilogram mobilogram, picked_mobilogram, smoothed_mobilogram;
  
  mobilogram = get_mobilogram(0);
  PeakPickerMobilogram picker;
  Param picker_param = picker.getParameters();
  picker_param.setValue("method", "corrected");
  picker_param.setValue("use_gauss", "false");
  picker.setParameters(picker_param);
  picker.pickMobilogram(mobilogram, picked_mobilogram, smoothed_mobilogram);

  TEST_REAL_SIMILAR(picked_mobilogram[0].getIntensity(), 58956.1);
  TEST_REAL_SIMILAR(picked_mobilogram[0].getPos(), 0.978364);
  TEST_REAL_SIMILAR(picked_mobilogram.getFloatDataArrays()[PeakPickerMobilogram::IDX_ABUNDANCE][0], 884145); // IntegratedIntensity
  TEST_REAL_SIMILAR(picked_mobilogram.getFloatDataArrays()[PeakPickerMobilogram::IDX_LEFTBORDER][0], 0.948675); // leftWidth
  TEST_REAL_SIMILAR(picked_mobilogram.getFloatDataArrays()[PeakPickerMobilogram::IDX_RIGHTBORDER][0], 0.997081); // rightWidth
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST