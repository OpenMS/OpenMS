// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumHelper, "$Id$")

START_SECTION((MSSpectrum::FloatDataArrays::iterator getDataArrayByName(MSSpectrum::FloatDataArrays& a, const String& name)))
	MSSpectrum ds;
	MSSpectrum::FloatDataArray float_array;
        float_array.push_back(56);  
        float_array.push_back(201); 
        float_array.push_back(31);  
        float_array.push_back(201); 
        float_array.push_back(201);  
	
        ds.getFloatDataArrays() = std::vector<MSSpectrum::FloatDataArray>(3, float_array);
	ds.getFloatDataArrays()[0].setName("f1");
	ds.getFloatDataArrays()[1].setName("f2");
	ds.getFloatDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f2") == ds.getFloatDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "NOT_THERE") == ds.getFloatDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f1") - ds.getFloatDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f2") - ds.getFloatDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f3") - ds.getFloatDataArrays().begin(), 2);

        // test const version
        const MSSpectrum::FloatDataArrays cfa(ds.getFloatDataArrays());
        TEST_EQUAL(getDataArrayByName(cfa, "f2") == cfa.end(), false);
        TEST_EQUAL(getDataArrayByName(cfa, "f1") - cfa.begin(), 0);
        TEST_EQUAL(getDataArrayByName(cfa, "NOT_THERE") == cfa.end(), true);

END_SECTION

START_SECTION((MSSpectrum::StringDataArrays::iterator getDataArrayByName(MSSpectrum::StringDataArrays& a, const String& name)))
	MSSpectrum ds;
	MSSpectrum::StringDataArray string_array;
        string_array.push_back("56");  
        string_array.push_back("201"); 
        string_array.push_back("31");  
        string_array.push_back("201"); 
        string_array.push_back("201");  
	ds.getStringDataArrays() = std::vector<MSSpectrum::StringDataArray>(3, string_array);
	ds.getStringDataArrays()[0].setName("f1");
	ds.getStringDataArrays()[1].setName("f2");
	ds.getStringDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f2") == ds.getStringDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "NOT_THERE") == ds.getStringDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f1") - ds.getStringDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f2") - ds.getStringDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f3") - ds.getStringDataArrays().begin(), 2);
END_SECTION

START_SECTION((MSSpectrum::IntegerDataArrays::iterator getDataArrayByName(MSSpectrum::IntegerDataArrays& a, const String& name)))
	MSSpectrum ds;
	MSSpectrum::IntegerDataArray int_array;
        int_array.push_back(56);  
        int_array.push_back(201); 
        int_array.push_back(31);  
        int_array.push_back(201); 
        int_array.push_back(201);  

	ds.getIntegerDataArrays() = std::vector<MSSpectrum::IntegerDataArray>(3, int_array);
	ds.getIntegerDataArrays()[0].setName("f1");
	ds.getIntegerDataArrays()[1].setName("f2");
	ds.getIntegerDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f2") == ds.getIntegerDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "NOT_THERE") == ds.getIntegerDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f1") - ds.getIntegerDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f2") - ds.getIntegerDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f3") - ds.getIntegerDataArrays().begin(), 2);
END_SECTION

START_SECTION((MSChromatogram::FloatDataArrays::iterator getDataArrayByName(const MSChromatogram::FloatDataArrays& a, const String& name)))
	MSChromatogram ds;
	MSChromatogram::FloatDataArray float_array;
        float_array.push_back(56);  
        float_array.push_back(201); 
        float_array.push_back(31);  
        float_array.push_back(201); 
        float_array.push_back(201);  
	
        ds.getFloatDataArrays() = std::vector<MSChromatogram::FloatDataArray>(3, float_array);
	ds.getFloatDataArrays()[0].setName("f1");
	ds.getFloatDataArrays()[1].setName("f2");
	ds.getFloatDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f2") == ds.getFloatDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "NOT_THERE") == ds.getFloatDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f1") - ds.getFloatDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f2") - ds.getFloatDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f3") - ds.getFloatDataArrays().begin(), 2);
END_SECTION

START_SECTION((MSChromatogram::StringDataArrays::iterator getDataArrayByName(MSChromatogram::StringDataArrays& a, const String& name)))
	MSChromatogram ds;
	MSChromatogram::StringDataArray string_array;
        string_array.push_back("56");  
        string_array.push_back("201"); 
        string_array.push_back("31");  
        string_array.push_back("201"); 
        string_array.push_back("201");  
	ds.getStringDataArrays() = std::vector<MSChromatogram::StringDataArray>(3, string_array);
	ds.getStringDataArrays()[0].setName("f1");
	ds.getStringDataArrays()[1].setName("f2");
	ds.getStringDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f2") == ds.getStringDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "NOT_THERE") == ds.getStringDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f1") - ds.getStringDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f2") - ds.getStringDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f3") - ds.getStringDataArrays().begin(), 2);
END_SECTION

START_SECTION((MSChromatogram::IntegerDataArrays::iterator getDataArrayByName(MSChromatogram::IntegerDataArrays& a, const String& name)))
	MSChromatogram ds;
	MSChromatogram::IntegerDataArray int_array;
        int_array.push_back(56);  
        int_array.push_back(201); 
        int_array.push_back(31);  
        int_array.push_back(201); 
        int_array.push_back(201);  

	ds.getIntegerDataArrays() = std::vector<MSChromatogram::IntegerDataArray>(3, int_array);
	ds.getIntegerDataArrays()[0].setName("f1");
	ds.getIntegerDataArrays()[1].setName("f2");
	ds.getIntegerDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f2") == ds.getIntegerDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "NOT_THERE") == ds.getIntegerDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f1") - ds.getIntegerDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f2") - ds.getIntegerDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f3") - ds.getIntegerDataArrays().begin(), 2);
END_SECTION

START_SECTION(void removePeaks(PeakContainerT& p, const double pos_start, const double pos_end))
{
  MSSpectrum s;
  MSChromatogram c;
  DataArrays::IntegerDataArray ida;

  for (Size i = 5; i < 15; ++i) // RTs: [5 14]
  {
    s.push_back(Peak1D(i, 0));
    c.push_back(ChromatogramPeak(i, 0));
    ida.push_back(i);
  }

  s.getIntegerDataArrays().push_back(ida);
  c.getIntegerDataArrays().push_back(ida);

  MSSpectrum s1 = s;
  removePeaks(s1, 3, 6); // start rt (3) is lower than the minimum (5) within the spectrum
  TEST_EQUAL(s1.size(), 2)
  TEST_EQUAL(s1.getIntegerDataArrays()[0].size(), 2)
  TEST_REAL_SIMILAR(s1[0].getPos(), 5)
  TEST_REAL_SIMILAR(s1[1].getPos(), 6)

  MSSpectrum s2 = s;
  removePeaks(s2, 0, 4); // no peak within the requested range
  TEST_EQUAL(s2.size(), 0)
  TEST_EQUAL(s2.getIntegerDataArrays()[0].size(), 0)

  MSChromatogram c1 = c;
  removePeaks(c1, 12, 16); // end rt (16) is higher than the maximum (14) within the chromatogram
  TEST_EQUAL(c1.size(), 3)
  TEST_EQUAL(c1.getIntegerDataArrays()[0].size(), 3)
  TEST_REAL_SIMILAR(c1[0].getPos(), 12)
  TEST_REAL_SIMILAR(c1[1].getPos(), 13)
  TEST_REAL_SIMILAR(c1[2].getPos(), 14)

  MSChromatogram c2 = c;
  removePeaks(c2, 9, 12); // all within the range
  TEST_EQUAL(c2.size(), 4)
  TEST_EQUAL(c2.getIntegerDataArrays()[0].size(), 4)
  TEST_REAL_SIMILAR(c2[0].getPos(), 9)
  TEST_REAL_SIMILAR(c2[1].getPos(), 10)
  TEST_REAL_SIMILAR(c2[2].getPos(), 11)
  TEST_REAL_SIMILAR(c2[3].getPos(), 12)

  MSSpectrum s_empty;
  removePeaks(s_empty, 9, 12);
  TEST_EQUAL(s_empty.size(), 0)
}
END_SECTION

START_SECTION(void subtractMinimumIntensity(PeakContainerT& p))
{
  MSSpectrum s;
  MSChromatogram c;
  
  for (Int i = -5; i < 5; ++i) // Intensities: [-5 4]
  {
    s.push_back(Peak1D(0, i));
    c.push_back(ChromatogramPeak(0, i));
  }
  
  subtractMinimumIntensity(s);
  TEST_REAL_SIMILAR(s[0].getIntensity(), 0)
  TEST_REAL_SIMILAR(s[1].getIntensity(), 1)
  TEST_REAL_SIMILAR(s[9].getIntensity(), 9)
  
  subtractMinimumIntensity(c);
  TEST_REAL_SIMILAR(c[0].getIntensity(), 0)
  TEST_REAL_SIMILAR(c[1].getIntensity(), 1)
  TEST_REAL_SIMILAR(c[9].getIntensity(), 9)
  
  MSChromatogram c_empty;
  subtractMinimumIntensity(c_empty);
  TEST_EQUAL(c_empty.size(), 0)
  
  s.clear(true);
  c.clear(true);
  
  for (Int i = 5; i < 15; ++i) // Intensities: [5 14]
  {
    s.push_back(Peak1D(0, i));
    c.push_back(ChromatogramPeak(0, i));
  }
  
  subtractMinimumIntensity(s);
  TEST_REAL_SIMILAR(s[0].getIntensity(), 0)
  TEST_REAL_SIMILAR(s[1].getIntensity(), 1)
  TEST_REAL_SIMILAR(s[9].getIntensity(), 9)
  
  subtractMinimumIntensity(c);
  TEST_REAL_SIMILAR(c[0].getIntensity(), 0)
  TEST_REAL_SIMILAR(c[1].getIntensity(), 1)
  TEST_REAL_SIMILAR(c[9].getIntensity(), 9)
}
END_SECTION

START_SECTION(void makePeakPositionUnique(PeakContainerT& p, const int m) )
{
  MSSpectrum s;
  MSChromatogram c;
  
  MSSpectrum s_template;
  MSChromatogram c_template;
  
  s_template.push_back(Peak1D(1, 1));
  s_template.push_back(Peak1D(2, 4));
  s_template.push_back(Peak1D(2, 8));
  s_template.push_back(Peak1D(3, 9));
  s_template.push_back(Peak1D(4, 7));
  s_template.push_back(Peak1D(2, 10));
  
  c_template.push_back(ChromatogramPeak(1, 1));
  c_template.push_back(ChromatogramPeak(2, 4));
  c_template.push_back(ChromatogramPeak(2, 8));
  c_template.push_back(ChromatogramPeak(3, 9));
  c_template.push_back(ChromatogramPeak(4, 7));
  c_template.push_back(ChromatogramPeak(2, 10));

  s = s_template;
  makePeakPositionUnique(s, IntensityAveragingMethod::MEDIAN);
  TEST_REAL_SIMILAR(s[0].getIntensity(), 1)
  TEST_REAL_SIMILAR(s[1].getIntensity(), 8)
  TEST_REAL_SIMILAR(s[2].getIntensity(), 9)
  TEST_REAL_SIMILAR(s[3].getIntensity(), 7)
  
  c = c_template;
  makePeakPositionUnique(c, IntensityAveragingMethod::MEDIAN);
  TEST_REAL_SIMILAR(c[0].getIntensity(), 1)
  TEST_REAL_SIMILAR(c[1].getIntensity(), 8)
  TEST_REAL_SIMILAR(c[2].getIntensity(), 9)
  TEST_REAL_SIMILAR(c[3].getIntensity(), 7)
  
  s = s_template;
  makePeakPositionUnique(s, IntensityAveragingMethod::MEAN);
  TEST_REAL_SIMILAR(s[0].getIntensity(), 1)
  TEST_REAL_SIMILAR(s[1].getIntensity(), (4+8+10)/3.0)
  TEST_REAL_SIMILAR(s[2].getIntensity(), 9)
  TEST_REAL_SIMILAR(s[3].getIntensity(), 7)
  
  c = c_template;
  makePeakPositionUnique(c, IntensityAveragingMethod::MEAN);
  TEST_REAL_SIMILAR(c[0].getIntensity(), 1)
  TEST_REAL_SIMILAR(c[1].getIntensity(), (4+8+10)/3.0)
  TEST_REAL_SIMILAR(c[2].getIntensity(), 9)
  TEST_REAL_SIMILAR(c[3].getIntensity(), 7)
  
  s = s_template;
  makePeakPositionUnique(s, IntensityAveragingMethod::SUM);
  TEST_REAL_SIMILAR(s[0].getIntensity(), 1)
  TEST_REAL_SIMILAR(s[1].getIntensity(), (4+8+10))
  TEST_REAL_SIMILAR(s[2].getIntensity(), 9)
  TEST_REAL_SIMILAR(s[3].getIntensity(), 7)
  
  c = c_template;
  makePeakPositionUnique(c, IntensityAveragingMethod::SUM);
  TEST_REAL_SIMILAR(c[0].getIntensity(), 1)
  TEST_REAL_SIMILAR(c[1].getIntensity(), (4+8+10))
  TEST_REAL_SIMILAR(c[2].getIntensity(), 9)
  TEST_REAL_SIMILAR(c[3].getIntensity(), 7)
  
  s = s_template;
  makePeakPositionUnique(s, IntensityAveragingMethod::MIN);
  TEST_REAL_SIMILAR(s[0].getIntensity(), 1)
  TEST_REAL_SIMILAR(s[1].getIntensity(), 4)
  TEST_REAL_SIMILAR(s[2].getIntensity(), 9)
  TEST_REAL_SIMILAR(s[3].getIntensity(), 7)
  
  c = c_template;
  makePeakPositionUnique(c, IntensityAveragingMethod::MIN);
  TEST_REAL_SIMILAR(c[0].getIntensity(), 1)
  TEST_REAL_SIMILAR(c[1].getIntensity(), 4)
  TEST_REAL_SIMILAR(c[2].getIntensity(), 9)
  TEST_REAL_SIMILAR(c[3].getIntensity(), 7)
  
  s = s_template;
  makePeakPositionUnique(s, IntensityAveragingMethod::MAX);
  TEST_REAL_SIMILAR(s[0].getIntensity(), 1)
  TEST_REAL_SIMILAR(s[1].getIntensity(), 10)
  TEST_REAL_SIMILAR(s[2].getIntensity(), 9)
  TEST_REAL_SIMILAR(s[3].getIntensity(), 7)
  
  c = c_template;
  makePeakPositionUnique(c, IntensityAveragingMethod::MAX);
  TEST_REAL_SIMILAR(c[0].getIntensity(), 1)
  TEST_REAL_SIMILAR(c[1].getIntensity(), 10)
  TEST_REAL_SIMILAR(c[2].getIntensity(), 9)
  TEST_REAL_SIMILAR(c[3].getIntensity(), 7)
  
  MSSpectrum s_empty;
  makePeakPositionUnique(s_empty);
  TEST_EQUAL(s_empty.size(), 0)
  
  MSChromatogram c_empty;
  makePeakPositionUnique(c_empty);
  TEST_EQUAL(c_empty.size(), 0)
}
END_SECTION

END_TEST

