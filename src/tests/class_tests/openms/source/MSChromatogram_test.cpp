// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_WINDOWSPLATFORM
#pragma clang diagnostic push
// Ignore -Wpessimizing-move, becuase it's intentional
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/MSChromatogram.h>
///////////////////////////

#include <sstream>

using namespace OpenMS;
using namespace std;

START_TEST(MSChromatogram, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSChromatogram* ptr = nullptr;
MSChromatogram* nullPointer = nullptr;
START_SECTION(MSChromatogram())
{
  ptr = new MSChromatogram();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~MSChromatogram())
{
  delete ptr;
}
END_SECTION

ptr = new MSChromatogram();
ChromatogramPeak p1;
p1.setIntensity(1.0f);
p1.setRT(2.0);

ChromatogramPeak p2;
p2.setIntensity(2.0f);
p2.setRT(10.0);

ChromatogramPeak p3;
p3.setIntensity(3.0f);
p3.setRT(30.0);

ChromatogramPeak p4;
p4.setIntensity(6.0f);
p4.setRT(30);

ChromatogramPeak p5;
p5.setIntensity(3.0f);
p5.setRT(30.0001);


START_SECTION((const String& getName() const ))
{
  TEST_STRING_EQUAL(ptr->getName(), "")
  ptr->setName("my_fancy_name");
  TEST_STRING_EQUAL(ptr->getName(), "my_fancy_name")
  delete ptr;
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((const FloatDataArrays& getFloatDataArrays() const ))
{
  MSChromatogram chrom;
  TEST_EQUAL(chrom.getFloatDataArrays().size(),0)
}
END_SECTION

START_SECTION((FloatDataArrays& getFloatDataArrays()))
{
  MSChromatogram chrom;
  chrom.getFloatDataArrays().resize(2);
  TEST_EQUAL(chrom.getFloatDataArrays().size(), 2)
}
END_SECTION

START_SECTION((const StringDataArrays& getStringDataArrays() const ))
{
  MSChromatogram chrom;
  TEST_EQUAL(chrom.getStringDataArrays().size(),0)
}
END_SECTION

START_SECTION((StringDataArrays& getStringDataArrays()))
{
  MSChromatogram chrom;
  chrom.getStringDataArrays().resize(2);
  TEST_EQUAL(chrom.getStringDataArrays().size(), 2)
}
END_SECTION

START_SECTION((const IntegerDataArrays& getIntegerDataArrays() const ))
{
  MSChromatogram chrom;
  TEST_EQUAL(chrom.getIntegerDataArrays().size(), 0)
}
END_SECTION

START_SECTION((IntegerDataArrays& getIntegerDataArrays()))
{
  MSChromatogram chrom;
  chrom.getIntegerDataArrays().resize(2);
  TEST_EQUAL(chrom.getIntegerDataArrays().size(), 2)
}
END_SECTION

START_SECTION((void sortByIntensity(bool reverse=false)))
{
  MSChromatogram ds;
  ChromatogramPeak p;
  MSChromatogram::FloatDataArray float_array;
  MSChromatogram::StringDataArray string_array;
  MSChromatogram::IntegerDataArray int_array;
  std::vector<double> rts, intensities;
  MSChromatogram::IntegerDataArray in_array;
  intensities.push_back(201); rts.push_back(420.130); float_array.push_back(420.130f); string_array.push_back("420.13"); int_array.push_back(420);
  intensities.push_back(60);  rts.push_back(412.824); float_array.push_back(412.824f); string_array.push_back("412.82"); int_array.push_back(412);
  intensities.push_back(56);  rts.push_back(423.269); float_array.push_back(423.269f); string_array.push_back("423.27"); int_array.push_back(423);
  intensities.push_back(37);  rts.push_back(415.287); float_array.push_back(415.287f); string_array.push_back("415.29"); int_array.push_back(415);
  intensities.push_back(34);  rts.push_back(413.800); float_array.push_back(413.800f); string_array.push_back("413.80"); int_array.push_back(413);
  intensities.push_back(31);  rts.push_back(419.113); float_array.push_back(419.113f); string_array.push_back("419.11"); int_array.push_back(419);
  intensities.push_back(31);  rts.push_back(416.293); float_array.push_back(416.293f); string_array.push_back("416.29"); int_array.push_back(416);
  intensities.push_back(31);  rts.push_back(418.232); float_array.push_back(418.232f); string_array.push_back("418.23"); int_array.push_back(418);
  intensities.push_back(29);  rts.push_back(414.301); float_array.push_back(414.301f); string_array.push_back("414.30"); int_array.push_back(414);
  intensities.push_back(29);  rts.push_back(412.321); float_array.push_back(412.321f); string_array.push_back("412.32"); int_array.push_back(412);

  for (Size i = 0; i < rts.size(); ++i)
  {
    p.setIntensity(intensities[i]); 
    p.setRT(rts[i]);
    ds.push_back(p);
  }
  ds.sortByIntensity();
  std::vector<double> intensities_copy(intensities);
  std::sort(intensities_copy.begin(),intensities_copy.end());
  MSChromatogram::iterator it_ds = ds.begin();
  for(std::vector<double>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
  {
    if(it_ds == ds.end())
    {
      TEST_EQUAL(true,false)
    }
    TEST_EQUAL(it_ds->getIntensity(), *it);
    ++it_ds;
  }
  
  ds.clear(true);
  for (Size i = 0; i < rts.size(); ++i)
  {
    p.setIntensity(intensities[i]); 
    p.setRT(rts[i]);
    ds.push_back(p);
  }
  intensities_copy = intensities;
  std::sort(intensities_copy.begin(),intensities_copy.end());

  ds.getFloatDataArrays() = std::vector<MSChromatogram::FloatDataArray>(3,float_array);
  ds.getFloatDataArrays()[0].setName("f1");
  ds.getFloatDataArrays()[1].setName("f2");
  ds.getFloatDataArrays()[2].setName("f3");

  ds.getStringDataArrays() = std::vector<MSChromatogram::StringDataArray>(2, string_array);
  ds.getStringDataArrays()[0].setName("s1");
  ds.getStringDataArrays()[1].setName("s2");

  ds.getIntegerDataArrays() = std::vector<MSChromatogram::IntegerDataArray>(1, int_array);
  ds.getIntegerDataArrays()[0].setName("i1");

  ds.sortByIntensity();

  TEST_STRING_EQUAL(ds.getFloatDataArrays()[0].getName(),"f1")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[1].getName(),"f2")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[2].getName(),"f3")

  TEST_STRING_EQUAL(ds.getStringDataArrays()[0].getName(),"s1")
  TEST_STRING_EQUAL(ds.getStringDataArrays()[1].getName(),"s2")

  TEST_STRING_EQUAL(ds.getIntegerDataArrays()[0].getName(),"i1")

  MSChromatogram::iterator it1 = ds.begin();
  MSChromatogram::FloatDataArray::iterator it2 = ds.getFloatDataArrays()[1].begin();
  MSChromatogram::StringDataArray::iterator it3 = ds.getStringDataArrays()[0].begin();
  MSChromatogram::IntegerDataArray::iterator it4 = ds.getIntegerDataArrays()[0].begin();
  TOLERANCE_ABSOLUTE(0.0001)
  for(std::vector<double>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
  {
    if(it1 != ds.end() && it2 != ds.getFloatDataArrays()[1].end() && it3 != ds.getStringDataArrays()[0].end() && it4 != ds.getIntegerDataArrays()[0].end())
    {
      //metadataarray values == mz values
      TEST_REAL_SIMILAR(it1->getIntensity(), *it);
      TEST_REAL_SIMILAR(*it2 , it1->getRT());
      TEST_STRING_EQUAL(*it3 , String::number(it1->getRT(),2));
      TEST_EQUAL(*it4 , (Int)floor(it1->getRT()));
      ++it1;
      ++it2;
      ++it3;
      ++it4;
    }
    else
    {
      TEST_EQUAL(true,false)
    }
  }
}
END_SECTION

START_SECTION((void sortByPosition()))
{
  MSChromatogram ds;
  ChromatogramPeak p;
  MSChromatogram::FloatDataArray float_array;
  MSChromatogram::StringDataArray string_array;
  MSChromatogram::IntegerDataArray int_array;
  std::vector<double> rts, intensities;
  intensities.push_back(56);  rts.push_back(423.269); float_array.push_back(56);  string_array.push_back("56");  int_array.push_back(56);
  intensities.push_back(201); rts.push_back(420.130); float_array.push_back(201); string_array.push_back("201"); int_array.push_back(201);
  intensities.push_back(31);  rts.push_back(419.113); float_array.push_back(31);  string_array.push_back("31");  int_array.push_back(31);
  intensities.push_back(31);  rts.push_back(418.232); float_array.push_back(31);  string_array.push_back("31");  int_array.push_back(31);
  intensities.push_back(31);  rts.push_back(416.293); float_array.push_back(31);  string_array.push_back("31");  int_array.push_back(31);
  intensities.push_back(37);  rts.push_back(415.287); float_array.push_back(37);  string_array.push_back("37");  int_array.push_back(37);
  intensities.push_back(29);  rts.push_back(414.301); float_array.push_back(29);  string_array.push_back("29");  int_array.push_back(29);
  intensities.push_back(34);  rts.push_back(413.800); float_array.push_back(34);  string_array.push_back("34");  int_array.push_back(34);
  intensities.push_back(60);  rts.push_back(412.824); float_array.push_back(60);  string_array.push_back("60");  int_array.push_back(60);
  intensities.push_back(29);  rts.push_back(412.321); float_array.push_back(29);  string_array.push_back("29");  int_array.push_back(29);

  for (Size i = 0; i < rts.size(); ++i)
  {
    p.setIntensity(intensities[i]); p.setRT(rts[i]);
    ds.push_back(p);
  }
  ds.sortByPosition();
  MSChromatogram::iterator it = ds.begin();
  for(std::vector<double>::reverse_iterator rit = intensities.rbegin(); rit != intensities.rend(); ++rit)
  {
    if(it == ds.end())
    {
      TEST_EQUAL(true,false)
    }
    TEST_EQUAL(it->getIntensity(), *rit);
    ++it;
  }
  ds.clear(true);
  for (Size i = 0; i < rts.size(); ++i)
  {
    p.setIntensity(intensities[i]); p.setRT(rts[i]);
    ds.push_back(p);
  }
  ds.getFloatDataArrays() = std::vector<MSChromatogram::FloatDataArray>(3,float_array);
  ds.getFloatDataArrays()[0].setName("f1");
  ds.getFloatDataArrays()[1].setName("f2");
  ds.getFloatDataArrays()[2].setName("f3");

  ds.getStringDataArrays() = std::vector<MSChromatogram::StringDataArray>(2, string_array);
  ds.getStringDataArrays()[0].setName("s1");
  ds.getStringDataArrays()[1].setName("s2");

  ds.getIntegerDataArrays() = std::vector<MSChromatogram::IntegerDataArray>(2, int_array);
  ds.getIntegerDataArrays()[0].setName("i1");

  ds.sortByPosition();

  TEST_STRING_EQUAL(ds.getFloatDataArrays()[0].getName(),"f1")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[1].getName(),"f2")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[2].getName(),"f3")

  TEST_STRING_EQUAL(ds.getStringDataArrays()[0].getName(),"s1")
  TEST_STRING_EQUAL(ds.getStringDataArrays()[1].getName(),"s2")

  TEST_STRING_EQUAL(ds.getIntegerDataArrays()[0].getName(),"i1")

  MSChromatogram::iterator it1 = ds.begin();
  MSChromatogram::FloatDataArray::iterator it2 = ds.getFloatDataArrays()[1].begin();
  MSChromatogram::StringDataArray::iterator it3 = ds.getStringDataArrays()[0].begin();
  MSChromatogram::IntegerDataArray::iterator it4 = ds.getIntegerDataArrays()[0].begin();
  for(std::vector<double>::reverse_iterator rit = intensities.rbegin(); rit != intensities.rend(); ++rit)
  {
    if(it1 != ds.end() && it2 != ds.getFloatDataArrays()[1].end() && it3 != ds.getStringDataArrays()[0].end())
    {
      //metadataarray values == intensity values
      TEST_REAL_SIMILAR(it1->getIntensity(), *rit);
      TEST_REAL_SIMILAR(*it2 , *rit);
      TEST_STRING_EQUAL(*it3 , String::number(*rit,0));
      TEST_EQUAL(*it4 , (Int)floor(*rit));
      ++it1;
      ++it2;
      ++it3;
      ++it4;
    }
    else
    {
      TEST_EQUAL(true,false)
    }
  }

}
END_SECTION

START_SECTION((bool isSorted() const ))
{
  //make test dataset
  MSChromatogram spec;
  ChromatogramPeak p;
  p.setIntensity(1.0);
  p.setRT(1000.0);
  spec.push_back(p);

  p.setIntensity(1.0);
  p.setRT(1001.0);
  spec.push_back(p);

  p.setIntensity(1.0);
  p.setRT(1002.0);
  spec.push_back(p);

  TEST_EQUAL(spec.isSorted(),true)

  reverse(spec.begin(), spec.end());
  TEST_EQUAL(spec.isSorted(),false)

}
END_SECTION

START_SECTION((Size findNearest(CoordinateType rt) const ))
{
  MSChromatogram tmp;
  ChromatogramPeak p;
  p.setIntensity(29.0f); p.setRT(412.321); tmp.push_back(p); //0
  p.setIntensity(60.0f); p.setRT(412.824); tmp.push_back(p); //1
  p.setIntensity(34.0f); p.setRT(413.8); tmp.push_back(p); //2
  p.setIntensity(29.0f); p.setRT(414.301); tmp.push_back(p); //3
  p.setIntensity(37.0f); p.setRT(415.287); tmp.push_back(p); //4
  p.setIntensity(31.0f); p.setRT(416.293); tmp.push_back(p); //5
  p.setIntensity(31.0f); p.setRT(418.232); tmp.push_back(p); //6
  p.setIntensity(31.0f); p.setRT(419.113); tmp.push_back(p); //7
  p.setIntensity(201.0f); p.setRT(420.13); tmp.push_back(p); //8
  p.setIntensity(56.0f); p.setRT(423.269); tmp.push_back(p); //9
  p.setIntensity(34.0f); p.setRT(426.292); tmp.push_back(p); //10
  p.setIntensity(82.0f); p.setRT(427.28); tmp.push_back(p); //11
  p.setIntensity(87.0f); p.setRT(428.322); tmp.push_back(p); //12
  p.setIntensity(30.0f); p.setRT(430.269); tmp.push_back(p); //13
  p.setIntensity(29.0f); p.setRT(431.246); tmp.push_back(p); //14
  p.setIntensity(42.0f); p.setRT(432.289); tmp.push_back(p); //15
  p.setIntensity(32.0f); p.setRT(436.161); tmp.push_back(p); //16
  p.setIntensity(54.0f); p.setRT(437.219); tmp.push_back(p); //17
  p.setIntensity(40.0f); p.setRT(439.186); tmp.push_back(p); //18
  p.setIntensity(40); p.setRT(440.27); tmp.push_back(p); //19
  p.setIntensity(23.0f); p.setRT(441.224); tmp.push_back(p); //20

  //test outside mass range
  TEST_EQUAL(tmp.findNearest(400.0),0);
  TEST_EQUAL(tmp.findNearest(500.0),20);
  //test mass range borders
  TEST_EQUAL(tmp.findNearest(412.4),0);
  TEST_EQUAL(tmp.findNearest(441.224),20);
  //test inside scan
  TEST_EQUAL(tmp.findNearest(426.29),10);
  TEST_EQUAL(tmp.findNearest(426.3),10);
  TEST_EQUAL(tmp.findNearest(427.2),11);
  TEST_EQUAL(tmp.findNearest(427.3),11);

  //empty spectrum
  MSChromatogram tmp2;
  TEST_PRECONDITION_VIOLATED(tmp2.findNearest(427.3));
}
END_SECTION


START_SECTION((Iterator RTBegin(CoordinateType rt)))
{
  MSChromatogram tmp;
  MSChromatogram::PeakType rdp;
  rdp.getPosition()[0] = 1.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 2.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 3.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 4.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 5.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 6.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 7.0;
  tmp.push_back(rdp);

  MSChromatogram::Iterator it;

  it = tmp.RTBegin(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.RTBegin(5.0);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.RTBegin(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)
}
END_SECTION

START_SECTION((Iterator RTBegin(Iterator begin, CoordinateType rt, Iterator end)))
{
  MSChromatogram tmp;
  MSChromatogram::PeakType rdp;
  rdp.getPosition()[0] = 1.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 2.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 3.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 4.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 5.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 6.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 7.0;
  tmp.push_back(rdp);

  MSChromatogram::Iterator it;

  it = tmp.RTBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.RTBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.RTBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])

}
END_SECTION

START_SECTION((Iterator RTEnd(CoordinateType rt)))
{
  MSChromatogram tmp;
  MSChromatogram::PeakType rdp;
  rdp.getPosition()[0] = 1.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 2.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 3.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 4.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 5.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 6.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 7.0;
  tmp.push_back(rdp);

  MSChromatogram::Iterator it;

  it = tmp.RTEnd(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.RTEnd(5.0);
  TEST_EQUAL(it->getPosition()[0],6.0)
  it = tmp.RTEnd(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)
}
END_SECTION

START_SECTION((Iterator RTEnd(Iterator begin, CoordinateType rt, Iterator end)))
{
  MSChromatogram tmp;
  MSChromatogram::PeakType rdp;
  rdp.getPosition()[0] = 1.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 2.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 3.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 4.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 5.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 6.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 7.0;
  tmp.push_back(rdp);

  MSChromatogram::Iterator it;

  it = tmp.RTEnd(tmp.begin(), 3.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],4.0)
  it = tmp.RTEnd(tmp.begin(), 5.0, tmp.end());
  TEST_EQUAL(it->getPosition()[0],6.0)
  it = tmp.RTEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
}
END_SECTION

START_SECTION((ConstIterator RTBegin(CoordinateType rt) const ))
{
  MSChromatogram tmp;
  MSChromatogram::PeakType rdp;
  rdp.getPosition()[0] = 1.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 2.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 3.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 4.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 5.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 6.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 7.0;
  tmp.push_back(rdp);

  MSChromatogram::ConstIterator it;

  it = tmp.RTBegin(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.RTBegin(5.0);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.RTBegin(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)

}
END_SECTION

START_SECTION((ConstIterator RTBegin(ConstIterator begin, CoordinateType rt, ConstIterator end) const ))
{
  MSChromatogram tmp;
  MSChromatogram::PeakType rdp;
  rdp.getPosition()[0] = 1.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 2.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 3.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 4.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 5.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 6.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 7.0;
  tmp.push_back(rdp);

  MSChromatogram::ConstIterator it;

  it = tmp.RTBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.RTBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.RTBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
}
END_SECTION

START_SECTION((ConstIterator RTEnd(CoordinateType rt) const ))
{
  MSChromatogram tmp;
  MSChromatogram::PeakType rdp;
  rdp.getPosition()[0] = 1.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 2.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 3.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 4.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 5.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 6.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 7.0;
  tmp.push_back(rdp);

  MSChromatogram::ConstIterator it;

  it = tmp.RTEnd(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.RTEnd(5.0);
  TEST_EQUAL(it->getPosition()[0],6.0)
  it = tmp.RTEnd(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)

}
END_SECTION

START_SECTION((ConstIterator RTEnd(ConstIterator begin, CoordinateType rt, ConstIterator end) const ))
{
  MSChromatogram tmp;
  MSChromatogram::PeakType rdp;
  rdp.getPosition()[0] = 1.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 2.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 3.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 4.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 5.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 6.0;
  tmp.push_back(rdp);
  rdp.getPosition()[0] = 7.0;
  tmp.push_back(rdp);

  MSChromatogram::ConstIterator it;

  it = tmp.RTEnd(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.RTEnd(tmp.begin(), 5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],6.0)
  it = tmp.RTEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
}
END_SECTION

MSChromatogram tmp;
vector<double> position = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
for (Size i=0; i<position.size(); ++i)
{
  ChromatogramPeak cp;
  cp.setPos(position[i]);
  tmp.push_back(cp);
}

START_SECTION((Iterator PosBegin(CoordinateType rt)))
{
  MSChromatogram::Iterator it;
  it = tmp.PosBegin(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.0);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((Iterator PosBegin(Iterator begin, CoordinateType rt, Iterator end)))
{
  MSChromatogram::Iterator it;
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(tmp.begin(), 5.5, tmp.end());
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it-1)->getPos(), (tmp.end()-1)->getPos())
}
END_SECTION

START_SECTION((ConstIterator PosBegin(CoordinateType rt) const ))
{
  MSChromatogram::ConstIterator it;
  it = tmp.PosBegin(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.0);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((ConstIterator PosBegin(ConstIterator begin, CoordinateType rt, ConstIterator end) const ))
{
  MSChromatogram::ConstIterator it;
  it = tmp.PosBegin(tmp.begin(), 3.5, tmp.end());
  TEST_EQUAL(it->getPos(), 4.0)
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it-1)->getPos(), (tmp.end()-1)->getPos())
}
END_SECTION

START_SECTION((Iterator PosEnd(CoordinateType rt)))
{
  MSChromatogram::Iterator it;
  it = tmp.PosEnd(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(5.0);
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosEnd(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((Iterator PosEnd(Iterator begin, CoordinateType rt, Iterator end)))
{
  MSChromatogram::Iterator it;
  it = tmp.PosEnd(tmp.begin(), 3.5, tmp.end());
  TEST_EQUAL(it->getPos(), 4.0)
  it = tmp.PosEnd(tmp.begin(), 4.0, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it-1)->getPos(), (tmp.end()-1)->getPos())
}
END_SECTION

START_SECTION((ConstIterator PosEnd(CoordinateType rt) const ))
{
  MSChromatogram::ConstIterator it;
  it = tmp.PosEnd(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(5.0);
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosEnd(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((ConstIterator PosEnd(ConstIterator begin, CoordinateType rt, ConstIterator end) const ))
{
  MSChromatogram::ConstIterator it;
  it = tmp.PosEnd(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(tmp.begin(), 5.0, tmp.end());
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it-1)->getPos(), (tmp.end()-1)->getPos())
}
END_SECTION

/////////////////////////////////////////////////////////////
// Copy constructor, move constructor, assignment operator, move assignment operator, equality

START_SECTION((MSChromatogram(const MSChromatogram &source)))
{
  MSChromatogram tmp;
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  tmp.setMetaValue("label",5.0);
  Product prod;
  prod.setMZ(7.0);
  tmp.setProduct(prod);
  tmp.setName("bla");
  //peaks
  MSChromatogram::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);

  MSChromatogram tmp2(tmp);
  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),1);
  TEST_REAL_SIMILAR(tmp2.getMetaValue("label"), 5.0)
  TEST_REAL_SIMILAR(tmp2.getMZ(), 7.0)
  TEST_EQUAL(tmp2.getName(),"bla")
  //peaks
  TEST_EQUAL(tmp2.size(),1);
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0],47.11);
  
}
END_SECTION

START_SECTION((MSChromatogram(const MSChromatogram &&source)))
{
  // Ensure that MSChromatogram has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(MSChromatogram(std::declval<MSChromatogram&&>())), true)

  MSChromatogram tmp;
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  tmp.setMetaValue("label",5.0);
  Product prod;
  prod.setMZ(7.0);
  tmp.setProduct(prod);
  tmp.setName("bla");
  //peaks
  MSChromatogram::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);

  //copy tmp so we can move one of them
  MSChromatogram orig = tmp;
  MSChromatogram tmp2(std::move(tmp));

  TEST_EQUAL(tmp2, orig); // should be equal to the original

  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),1);
  TEST_REAL_SIMILAR(tmp2.getMetaValue("label"), 5.0)
  TEST_REAL_SIMILAR(tmp2.getMZ(), 7.0)
  TEST_EQUAL(tmp2.getName(),"bla")
  //peaks
  TEST_EQUAL(tmp2.size(),1);
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0],47.11);

  // test move
  TEST_EQUAL(tmp.size(),0);
  TEST_EQUAL(tmp.metaValueExists("label"), false);
}
END_SECTION

START_SECTION((MSChromatogram& operator=(const MSChromatogram &source)))
{
  MSChromatogram tmp;
  tmp.setMetaValue("label",5.0);
  Product prod;
  prod.setMZ(7.0);
  tmp.setProduct(prod);
  tmp.setName("bla");
  //peaks
  MSChromatogram::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);

  //normal assignment
  MSChromatogram tmp2;
  tmp2 = tmp;
  TEST_REAL_SIMILAR(tmp2.getMetaValue("label"), 5.0)
  TEST_REAL_SIMILAR(tmp2.getMZ(), 7.0)
  TEST_EQUAL(tmp2.getName(),"bla")
  TEST_EQUAL(tmp2.size(),1);
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0],47.11);

  //Assignment of empty object
  //normal assignment
  tmp2 = MSChromatogram();
  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),0);
  TEST_EQUAL(tmp2.metaValueExists("label"), false)
  TEST_REAL_SIMILAR(tmp2.getMZ(), 0.0)
  TEST_EQUAL(tmp2.getName(),"")
  TEST_EQUAL(tmp2.size(),0);
}
END_SECTION

START_SECTION((MSChromatogram& operator=(const MSChromatogram &&source)))
{
  MSChromatogram tmp;
  tmp.setMetaValue("label",5.0);
  Product prod;
  prod.setMZ(7.0);
  tmp.setProduct(prod);
  tmp.setName("bla");
  //peaks
  MSChromatogram::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);

  //copy tmp so we can move one of them
  MSChromatogram orig = tmp;

  //move assignment
  MSChromatogram tmp2;
  tmp2 = std::move(tmp);

  TEST_EQUAL(tmp2, orig); // should be equal to the original

  //normal assignment
  TEST_REAL_SIMILAR(tmp2.getMetaValue("label"), 5.0)
  TEST_REAL_SIMILAR(tmp2.getMZ(), 7.0)
  TEST_EQUAL(tmp2.getName(),"bla")
  TEST_EQUAL(tmp2.size(),1);
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0],47.11);

  // test move
  TEST_EQUAL(tmp.size(),0);
  TEST_EQUAL(tmp.metaValueExists("label"), false);

  //Assignment of empty object
  //normal assignment
#ifndef OPENMS_WINDOWSPLATFORM
#pragma clang diagnostic push
// Ignore -Wpessimizing-move, because we want to test the move assignment operator.
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
  tmp2 = std::move(MSChromatogram());
#ifndef OPENMS_WINDOWSPLATFORM
#pragma clang diagnostic pop
#endif
  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),0);
  TEST_EQUAL(tmp2.metaValueExists("label"), false)
  TEST_REAL_SIMILAR(tmp2.getMZ(), 0.0)
  TEST_EQUAL(tmp2.getName(),"")
  TEST_EQUAL(tmp2.size(),0);
}
END_SECTION

START_SECTION((bool operator==(const MSChromatogram &rhs) const ))
{
  MSChromatogram edit, empty;

  TEST_EQUAL(edit==empty,true);

  edit = empty;
  edit.resize(1);
  TEST_EQUAL(edit==empty,false);

  edit = empty;
  edit.setMetaValue("label",String("bla"));
  TEST_EQUAL(empty==edit, false);

  Product prod;
  prod.setMZ(5);
  edit.setProduct(prod);
  TEST_EQUAL(empty==edit, false);

  edit = empty;
  edit.getFloatDataArrays().resize(5);
  TEST_EQUAL(empty==edit, false);

  edit = empty;
  edit.getStringDataArrays().resize(5);
  TEST_EQUAL(empty==edit, false);

  edit = empty;
  edit.getIntegerDataArrays().resize(5);
  TEST_EQUAL(empty==edit, false);

  //name is not checked => no change
  edit = empty;
  edit.setName("bla");
  TEST_TRUE(empty == edit);

  edit = empty;
  edit.push_back(p1);
  edit.push_back(p2);
  edit.updateRanges();
  edit.clear(false);
  TEST_EQUAL(empty==edit, false);

}
END_SECTION

START_SECTION((bool operator!=(const MSChromatogram &rhs) const ))
{
  MSChromatogram edit, empty;

  TEST_EQUAL(edit!=empty,false);

  edit = empty;
  edit.getInstrumentSettings().getScanWindows().resize(1);
  TEST_EQUAL(edit!=empty,true);

  edit = empty;
  edit.resize(1);
  TEST_EQUAL(edit!=empty,true);

  edit = empty;
  edit.setMetaValue("label",String("bla"));
  TEST_EQUAL(edit!=empty,true);

  Product prod;
  prod.setMZ(5);
  edit.setProduct(prod);
  TEST_EQUAL(edit!=empty,true);

  edit = empty;
  edit.getFloatDataArrays().resize(5);
  TEST_EQUAL(edit!=empty,true);

  edit = empty;
  edit.getIntegerDataArrays().resize(5);
  TEST_EQUAL(edit!=empty,true);

  edit = empty;
  edit.getStringDataArrays().resize(5);
  TEST_EQUAL(edit!=empty,true);

  //name is not checked => no change
  edit = empty;
  edit.setName("bla");
  TEST_EQUAL(edit!=empty,false);

  edit = empty;
  edit.push_back(p1);
  edit.push_back(p2);
  edit.updateRanges();
  edit.clear(false);
  TEST_EQUAL(edit!=empty,true);
}
END_SECTION

START_SECTION((virtual void updateRanges()))
{
  MSChromatogram s;
  s.push_back(p1);
  s.push_back(p2);
  s.push_back(p1);

  s.updateRanges();
  s.updateRanges(); //second time to check the initialization

  TEST_REAL_SIMILAR(s.getMaxIntensity(),2)
  TEST_REAL_SIMILAR(s.getMinIntensity(), 1)
  TEST_REAL_SIMILAR(s.getMaxRT(),10)
  TEST_REAL_SIMILAR(s.getMinRT(),2)

  //test with only one peak

  s.clear(true);
  s.push_back(p1);
  s.updateRanges();
  TEST_REAL_SIMILAR(s.getMaxIntensity(), 1)
  TEST_REAL_SIMILAR(s.getMinIntensity(), 1)
  TEST_REAL_SIMILAR(s.getMaxRT(),2)
  TEST_REAL_SIMILAR(s.getMinRT(),2)
}
END_SECTION

START_SECTION(void clear(bool clear_meta_data))
{
  MSChromatogram edit;
  edit.getInstrumentSettings().getScanWindows().resize(1);
  edit.resize(1);
  edit.setMetaValue("label",String("bla"));
  edit.getProduct().setMZ(5);
  edit.getFloatDataArrays().resize(5);
  edit.getIntegerDataArrays().resize(5);
  edit.getStringDataArrays().resize(5);

  edit.clear(false);
  TEST_EQUAL(edit.size(),0)
  TEST_EQUAL(edit.empty(),true)
  TEST_EQUAL(edit == MSChromatogram(),false)

  edit.clear(true);
  TEST_EQUAL(edit.size(),0)
  TEST_EQUAL(edit.empty(),true)
  TEST_EQUAL(edit == MSChromatogram(),true)
}
END_SECTION

START_SECTION((double getMZ() const))
{
  MSChromatogram tmp;
  Product prod;
  prod.setMZ(0.1);
  TEST_REAL_SIMILAR(tmp.getMZ(), 0.0)
    tmp.setProduct(prod);
  TEST_REAL_SIMILAR(tmp.getMZ(), 0.1)
}
END_SECTION

START_SECTION(([MSChromatogram::MZLess] bool operator()(const MSChromatogram &a, const MSChromatogram &b) const))
{
    MSChromatogram a;
    Product pa;
    pa.setMZ(1000.0);
    a.setProduct(pa);

    MSChromatogram b;
    Product pb;
    pb.setMZ(1000.1);
    b.setProduct(pb);

    TEST_EQUAL(MSChromatogram::MZLess().operator ()(a,b), true)
    TEST_EQUAL(MSChromatogram::MZLess().operator ()(b,a), false)

    TEST_EQUAL(MSChromatogram::MZLess().operator ()(a,a), false)
}
END_SECTION

START_SECTION(void mergePeaks(MSChromatogram & other) )
{
  MSChromatogram a, b, c;
  a.push_back(p1);
  b.push_back(p2);
  a.push_back(p3);
  b.push_back(p5);
  c.push_back(p1);
  c.push_back(p2);
  c.push_back(p4);
  a.sortByPosition();
  b.sortByPosition();
  c.sortByPosition();
  a.mergePeaks(b, true);
  DoubleList dl = { b.getMZ() };
  c.setMetaValue(Constants::UserParam::MERGED_CHROMATOGRAM_MZS, dl);
  TEST_EQUAL(a, c)
}
END_SECTION

START_SECTION( std::ostream& operator<<(std::ostream& os, const MSChromatogram& chrom)) 
{
  MSChromatogram a;
  Product pa;
  pa.setMZ(1000.0);
  a.setProduct(pa);

  MSChromatogram::PeakType peak;
  peak.getPosition()[0] = 47.11;
  a.push_back(peak);

  std::ostringstream os;
  os << a;

  TEST_EQUAL(String(os.str()).hasSubstring("MSCHROMATOGRAM BEGIN"), true);
  TEST_EQUAL(String(os.str()).hasSubstring("47.11"), true);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#ifndef OPENMS_WINDOWSPLATFORM
#pragma clang diagnostic pop
#endif
