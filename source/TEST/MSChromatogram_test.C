// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/KERNEL/MSChromatogram.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSChromatogram, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSChromatogram<ChromatogramPeak>* ptr = 0;
START_SECTION(MSChromatogram())
{
	ptr = new MSChromatogram<ChromatogramPeak>();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(virtual ~MSChromatogram())
{
	delete ptr;
}
END_SECTION

ptr = new MSChromatogram<ChromatogramPeak>();

START_SECTION((const String& getName() const ))
{
  TEST_STRING_EQUAL(ptr->getName(), "")
	ptr->setName("my_fancy_name");
	TEST_STRING_EQUAL(ptr->getName(), "my_fancy_name")
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((const FloatDataArrays& getFloatDataArrays() const ))
{
  MSChromatogram<> chrom;
	TEST_EQUAL(chrom.getFloatDataArrays().size(),0)
}
END_SECTION

START_SECTION((FloatDataArrays& getFloatDataArrays()))
{
	MSChromatogram<> chrom;
	chrom.getFloatDataArrays().resize(2);
	TEST_EQUAL(chrom.getFloatDataArrays().size(), 2)
}
END_SECTION

START_SECTION((const StringDataArrays& getStringDataArrays() const ))
{
  MSChromatogram<> chrom;
  TEST_EQUAL(chrom.getStringDataArrays().size(),0)
}
END_SECTION

START_SECTION((StringDataArrays& getStringDataArrays()))
{
  MSChromatogram<> chrom;
  chrom.getStringDataArrays().resize(2);
  TEST_EQUAL(chrom.getStringDataArrays().size(), 2)
}
END_SECTION

START_SECTION((const IntegerDataArrays& getIntegerDataArrays() const ))
{
	MSChromatogram<> chrom;
  TEST_EQUAL(chrom.getIntegerDataArrays().size(), 0)
}
END_SECTION

START_SECTION((IntegerDataArrays& getIntegerDataArrays()))
{
	MSChromatogram<> chrom;
  chrom.getIntegerDataArrays().resize(2);
  TEST_EQUAL(chrom.getIntegerDataArrays().size(), 2)
}
END_SECTION

START_SECTION((void sortByIntensity(bool reverse=false)))
{
  MSChromatogram<> ds;
  ChromatogramPeak p;
  MSChromatogram<>::FloatDataArray float_array;
  MSChromatogram<>::StringDataArray string_array;
  MSChromatogram<>::IntegerDataArray int_array;
  std::vector<DoubleReal> rts, intensities;
  MSChromatogram<>::IntegerDataArray in_array;
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
  std::vector<DoubleReal> intensities_copy(intensities);
  std::sort(intensities_copy.begin(),intensities_copy.end());
  MSChromatogram<>::iterator it_ds = ds.begin();
  for(std::vector<DoubleReal>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
  {
    if(it_ds == ds.end())
    {
      TEST_EQUAL(true,false)
    }
    TEST_EQUAL(it_ds->getIntensity(), *it);
    ++it_ds;
  }
  
	ds = MSChromatogram<>();
  for (Size i = 0; i < rts.size(); ++i)
  {
    p.setIntensity(intensities[i]); 
		p.setRT(rts[i]);
    ds.push_back(p);
  }
  intensities_copy = intensities;
  std::sort(intensities_copy.begin(),intensities_copy.end());

  ds.getFloatDataArrays() = std::vector<MSChromatogram<>::FloatDataArray>(3,float_array);
  ds.getFloatDataArrays()[0].setName("f1");
  ds.getFloatDataArrays()[1].setName("f2");
  ds.getFloatDataArrays()[2].setName("f3");

  ds.getStringDataArrays() = std::vector<MSChromatogram<>::StringDataArray>(2, string_array);
  ds.getStringDataArrays()[0].setName("s1");
  ds.getStringDataArrays()[1].setName("s2");

  ds.getIntegerDataArrays() = std::vector<MSChromatogram<>::IntegerDataArray>(1, int_array);
  ds.getIntegerDataArrays()[0].setName("i1");

  ds.sortByIntensity();

  TEST_STRING_EQUAL(ds.getFloatDataArrays()[0].getName(),"f1")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[1].getName(),"f2")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[2].getName(),"f3")

  TEST_STRING_EQUAL(ds.getStringDataArrays()[0].getName(),"s1")
  TEST_STRING_EQUAL(ds.getStringDataArrays()[1].getName(),"s2")

  TEST_STRING_EQUAL(ds.getIntegerDataArrays()[0].getName(),"i1")

  MSChromatogram<>::iterator it1 = ds.begin();
  MSChromatogram<>::FloatDataArray::iterator it2 = ds.getFloatDataArrays()[1].begin();
  MSChromatogram<>::StringDataArray::iterator it3 = ds.getStringDataArrays()[0].begin();
  MSChromatogram<>::IntegerDataArray::iterator it4 = ds.getIntegerDataArrays()[0].begin();
  TOLERANCE_ABSOLUTE(0.0001)
  for(std::vector<DoubleReal>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
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
  // TODO
}
END_SECTION

START_SECTION((bool isSorted() const ))
{
  // TODO
}
END_SECTION

START_SECTION((Size findNearest(CoordinateType mz) const ))
{
  // TODO
}
END_SECTION

START_SECTION((Iterator MZBegin(CoordinateType mz)))
{
  // TODO
}
END_SECTION

START_SECTION((Iterator MZBegin(Iterator begin, CoordinateType mz, Iterator end)))
{
  // TODO
}
END_SECTION

START_SECTION((Iterator MZEnd(CoordinateType mz)))
{
  // TODO
}
END_SECTION

START_SECTION((Iterator MZEnd(Iterator begin, CoordinateType mz, Iterator end)))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator MZBegin(CoordinateType mz) const ))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator MZBegin(ConstIterator begin, CoordinateType mz, ConstIterator end) const ))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator MZEnd(CoordinateType mz) const ))
{
  // TODO
}
END_SECTION

START_SECTION((ConstIterator MZEnd(ConstIterator begin, CoordinateType mz, ConstIterator end) const ))
{
  // TODO
}
END_SECTION

START_SECTION((MSChromatogram(const MSChromatogram &source)))
{
  // TODO
}
END_SECTION

START_SECTION((MSChromatogram& operator=(const MSChromatogram &source)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const MSChromatogram &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const MSChromatogram &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((virtual void updateRanges()))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



