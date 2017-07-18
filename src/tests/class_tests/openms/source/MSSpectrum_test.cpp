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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/MSSpectrum.h>
///////////////////////////

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

START_TEST(MSSpectrum, "$Id$")

/////////////////////////////////////////////////////////////
// Dummy peak data

Peak1D p1;
p1.setIntensity(1.0f);
p1.setMZ(2.0);

Peak1D p2;
p2.setIntensity(2.0f);
p2.setMZ(10.0);

Peak1D p3;
p3.setIntensity(3.0f);
p3.setMZ(30.0);


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSSpectrum<>* ptr = 0;
MSSpectrum<>* nullPointer = 0;
START_SECTION((MSSpectrum()))
  ptr = new MSSpectrum<>();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MSSpectrum()))
  delete ptr;
END_SECTION

START_SECTION(([EXTRA] MSSpectrum<>()))
  MSSpectrum<Peak1D > tmp;
  Peak1D peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);
  TEST_EQUAL(tmp.size(),1);
  TEST_REAL_SIMILAR(tmp[0].getMZ(),47.11);  
END_SECTION

/////////////////////////////////////////////////////////////
// Member accessors

START_SECTION((UInt getMSLevel() const))
  MSSpectrum<> spec;
  TEST_EQUAL(spec.getMSLevel(),1)
END_SECTION

START_SECTION((void setMSLevel(UInt ms_level)))
  MSSpectrum<> spec;
  spec.setMSLevel(17);
  TEST_EQUAL(spec.getMSLevel(),17)
END_SECTION

START_SECTION((const String& getName() const))
  MSSpectrum<> s;
  TEST_STRING_EQUAL(s.getName(),"")
END_SECTION

START_SECTION((void setName(const String &name)))
  MSSpectrum<> s;
  s.setName("bla");
  TEST_STRING_EQUAL(s.getName(),"bla")
END_SECTION

START_SECTION((double getRT() const ))
  MSSpectrum<> s;
  TEST_REAL_SIMILAR(s.getRT(),-1.0)
END_SECTION

START_SECTION((void setRT(double rt)))
  MSSpectrum<> s;
  s.setRT(0.451);
  TEST_REAL_SIMILAR(s.getRT(),0.451)
END_SECTION

START_SECTION((double getDriftTime() const ))
  MSSpectrum<> s;
  TEST_REAL_SIMILAR(s.getDriftTime(),-1.0)
END_SECTION

START_SECTION((void setDriftTime(double dt)))
  MSSpectrum<> s;
  s.setDriftTime(0.451);
  TEST_REAL_SIMILAR(s.getDriftTime(),0.451)
END_SECTION

START_SECTION((const FloatDataArrays& getFloatDataArrays() const))
  MSSpectrum<> s;
  TEST_EQUAL(s.getFloatDataArrays().size(),0)
END_SECTION

START_SECTION((FloatDataArrays& getFloatDataArrays()))
  MSSpectrum<> s;
  s.getFloatDataArrays().resize(2);
  TEST_EQUAL(s.getFloatDataArrays().size(),2)
END_SECTION

START_SECTION((const StringDataArrays& getStringDataArrays() const))
  MSSpectrum<> s;
  TEST_EQUAL(s.getStringDataArrays().size(),0)
END_SECTION

START_SECTION((StringDataArrays& getStringDataArrays()))
  MSSpectrum<> s;
  s.getStringDataArrays().resize(2);
  TEST_EQUAL(s.getStringDataArrays().size(),2)
END_SECTION

START_SECTION((const IntegerDataArrays& getIntegerDataArrays() const))
  MSSpectrum<> s;
  TEST_EQUAL(s.getIntegerDataArrays().size(),0)
END_SECTION

START_SECTION((IntegerDataArrays& getIntegerDataArrays()))
  MSSpectrum<> s;
  s.getIntegerDataArrays().resize(2);
  TEST_EQUAL(s.getIntegerDataArrays().size(),2)
END_SECTION

START_SECTION((MSSpectrum& select(const std::vector<Size>& indices)))
  MSSpectrum<> s;
  s.push_back(p1);
  s.push_back(p2);
  s.push_back(p3);
  s.push_back(p3);
  s.push_back(p2);

  int air[] = {1, 2, 3, 4, 5};
  float afr[] = {1.0, 2.0, 3.0, 4.0, 5.0};
  String asr[] = {"1", "2" , "3", "4", "5"};
  std::vector<int> ai(&air[0], &air[5]);
  MSSpectrum<>::IntegerDataArray aia;
  swap(aia, ai);
  std::vector<float> af(&afr[0], &afr[5]);
  MSSpectrum<>::FloatDataArray afa;
  swap(afa, af);
  std::vector<String> as(&asr[0], &asr[5]);
  MSSpectrum<>::StringDataArray asa;
  swap(asa, as);
  //MSSpectrum<>::IntegerDataArray
  s.getFloatDataArrays().push_back(afa);
  s.getIntegerDataArrays().push_back(aia);
  s.getStringDataArrays().push_back(asa);
  s.getFloatDataArrays().push_back(afa);
  s.getIntegerDataArrays().push_back(aia);
  s.getStringDataArrays().push_back(asa);

  TEST_REAL_SIMILAR(s[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(s[4].getIntensity(), 2.0)
  TEST_EQUAL(s.getFloatDataArrays().size(), 2)
  TEST_EQUAL(s.getFloatDataArrays()[0].size(), 5)
  TEST_EQUAL(s.getIntegerDataArrays().size(), 2)
  TEST_EQUAL(s.getIntegerDataArrays()[0].size(), 5)
  TEST_EQUAL(s.getStringDataArrays().size(), 2)
  TEST_EQUAL(s.getStringDataArrays()[0].size(), 5)

  // re-order
  MSSpectrum<> s2 = s;
  Size order[] = {4, 2, 3, 1, 0};
  s2.select(std::vector<Size>(&order[0], &order[5]));
  TEST_REAL_SIMILAR(s2[0].getIntensity(), 2.0)
  TEST_REAL_SIMILAR(s2[4].getIntensity(), 1.0)
  TEST_EQUAL(s2.getFloatDataArrays().size(), 2)
  TEST_EQUAL(s2.getFloatDataArrays()[0].size(), 5)
  TEST_EQUAL(s2.getIntegerDataArrays().size(), 2)
  TEST_EQUAL(s2.getIntegerDataArrays()[0].size(), 5)
  TEST_EQUAL(s2.getStringDataArrays().size(), 2)
  TEST_EQUAL(s2.getStringDataArrays()[0].size(), 5)
  
  TEST_REAL_SIMILAR(s2.getFloatDataArrays()[0][1], 3.0)
  TEST_EQUAL(s2.getIntegerDataArrays()[0][1], 3)
  TEST_EQUAL(s2.getStringDataArrays()[0][1], "3")

  // subset
  s2 = s;
  Size subset[] = {4, 2, 3};
  // --> new values in Meta arrays are: 
  //     5, 3, 4
  s2.select(std::vector<Size>(&subset[0], &subset[3]));
  TEST_REAL_SIMILAR(s2[0].getIntensity(), 2.0)
  TEST_REAL_SIMILAR(s2[1].getIntensity(), 3.0)
  TEST_REAL_SIMILAR(s2[2].getIntensity(), 3.0)
  TEST_EQUAL(s2.getFloatDataArrays().size(), 2)
  TEST_EQUAL(s2.getFloatDataArrays()[0].size(), 3)
  TEST_EQUAL(s2.getIntegerDataArrays().size(), 2)
  TEST_EQUAL(s2.getIntegerDataArrays()[0].size(), 3)
  TEST_EQUAL(s2.getStringDataArrays().size(), 2)
  TEST_EQUAL(s2.getStringDataArrays()[0].size(), 3)

  TEST_REAL_SIMILAR(s2.getFloatDataArrays()[0][1], 3.0)
  TEST_EQUAL(s2.getIntegerDataArrays()[0][1], 3)
  TEST_EQUAL(s2.getStringDataArrays()[0][1], "3")


END_SECTION

/////////////////////////////////////////////////////////////
// RangeManager

START_SECTION((virtual void updateRanges()))
  MSSpectrum<> s;
  s.push_back(p1);
  s.push_back(p2);
  s.push_back(p1);

  s.updateRanges();
  s.updateRanges(); //second time to check the initialization

  TEST_REAL_SIMILAR(s.getMaxInt(),2)
  TEST_REAL_SIMILAR(s.getMinInt(),1)
  TEST_REAL_SIMILAR(s.getMax()[0],10)
  TEST_REAL_SIMILAR(s.getMin()[0],2)

  //test with only one peak

  s.clear(true);
  s.push_back(p1);
  s.updateRanges();
  TEST_REAL_SIMILAR(s.getMaxInt(),1)
  TEST_REAL_SIMILAR(s.getMinInt(),1)
  TEST_REAL_SIMILAR(s.getMax()[0],2)
  TEST_REAL_SIMILAR(s.getMin()[0],2)
END_SECTION


/////////////////////////////////////////////////////////////
// Copy constructor, assignement operator, equality

START_SECTION((MSSpectrum(const MSSpectrum& source)))
  MSSpectrum<> tmp;
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  tmp.setMetaValue("label",5.0);
  tmp.setMSLevel(17);
  tmp.setRT(7.0);
  tmp.setDriftTime(8.0);
  tmp.setName("bla");
  //peaks
  MSSpectrum<>::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);
  
  MSSpectrum<> tmp2(tmp);
  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),1);
  TEST_REAL_SIMILAR(tmp2.getMetaValue("label"), 5.0)
  TEST_EQUAL(tmp2.getMSLevel(), 17)
  TEST_REAL_SIMILAR(tmp2.getRT(), 7.0)
  TEST_REAL_SIMILAR(tmp2.getDriftTime(), 8.0)
  TEST_EQUAL(tmp2.getName(),"bla")
  //peaks
  TEST_EQUAL(tmp2.size(),1);
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0],47.11);
END_SECTION


START_SECTION((MSSpectrum& operator= (const MSSpectrum& source)))
  MSSpectrum<> tmp;
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  tmp.setMetaValue("label",5.0);
  tmp.setMSLevel(17);
  tmp.setRT(7.0);
  tmp.setDriftTime(8.0);
  tmp.setName("bla");
  //peaks
  MSSpectrum<>::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);
  
  //normal assignment
  MSSpectrum<> tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),1);
  TEST_REAL_SIMILAR(tmp2.getMetaValue("label"), 5.0)
  TEST_EQUAL(tmp2.getMSLevel(), 17)
  TEST_REAL_SIMILAR(tmp2.getRT(), 7.0)
  TEST_REAL_SIMILAR(tmp2.getDriftTime(), 8.0)
  TEST_EQUAL(tmp2.getName(),"bla")
  TEST_EQUAL(tmp2.size(),1);
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0],47.11);
  
  //Assignment of empty object
  //normal assignment
  tmp2 = MSSpectrum<>();
  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),0);
  TEST_EQUAL(tmp2.metaValueExists("label"), false)
  TEST_EQUAL(tmp2.getMSLevel(),1)
  TEST_REAL_SIMILAR(tmp2.getRT(), -1.0)
  TEST_REAL_SIMILAR(tmp2.getDriftTime(), -1.0)
  TEST_EQUAL(tmp2.getName(),"")
  TEST_EQUAL(tmp2.size(),0);
END_SECTION

START_SECTION((bool operator== (const MSSpectrum& rhs) const))
  MSSpectrum<> edit, empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit = empty;
  edit.getInstrumentSettings().getScanWindows().resize(1);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.resize(1);
  TEST_EQUAL(edit==empty,false);

  edit = empty;
  edit.setMetaValue("label",String("bla"));
  TEST_EQUAL(empty==edit, false);

  edit = empty;
  edit.setDriftTime(5);
  TEST_EQUAL(empty==edit, false);

  edit = empty;
  edit.setRT(5);
  TEST_EQUAL(empty==edit, false);

  edit = empty;
  edit.setMSLevel(5);
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
  TEST_EQUAL(empty==edit, true);

  edit = empty;
  edit.push_back(p1);
  edit.push_back(p2);
  edit.updateRanges();
  edit.clear(false);
  TEST_EQUAL(empty==edit, false);
END_SECTION

START_SECTION((bool operator!= (const MSSpectrum& rhs) const))
  MSSpectrum<> edit, empty;
  
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

  edit = empty;
  edit.setDriftTime(5);
  TEST_EQUAL(edit!=empty,true);

  edit = empty;
  edit.setRT(5);
  TEST_EQUAL(edit!=empty,true);

  edit = empty;
  edit.setMSLevel(5);
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

END_SECTION


/////////////////////////////////////////////////////////////
// Sorting


START_SECTION((void sortByIntensity(bool reverse=false)))
  MSSpectrum<> ds;
  Peak1D p;
  MSSpectrum<>::FloatDataArray float_array;
  MSSpectrum<>::StringDataArray string_array;
  MSSpectrum<>::IntegerDataArray int_array;
  std::vector<double> mzs, intensities;
  MSSpectrum<>::IntegerDataArray in_array;
  intensities.push_back(201); mzs.push_back(420.130); float_array.push_back(420.130f); string_array.push_back("420.13"); int_array.push_back(420);
  intensities.push_back(60);  mzs.push_back(412.824); float_array.push_back(412.824f); string_array.push_back("412.82"); int_array.push_back(412);
  intensities.push_back(56);  mzs.push_back(423.269); float_array.push_back(423.269f); string_array.push_back("423.27"); int_array.push_back(423);
  intensities.push_back(37);  mzs.push_back(415.287); float_array.push_back(415.287f); string_array.push_back("415.29"); int_array.push_back(415);
  intensities.push_back(34);  mzs.push_back(413.800); float_array.push_back(413.800f); string_array.push_back("413.80"); int_array.push_back(413);
  intensities.push_back(31);  mzs.push_back(419.113); float_array.push_back(419.113f); string_array.push_back("419.11"); int_array.push_back(419);
  intensities.push_back(31);  mzs.push_back(416.293); float_array.push_back(416.293f); string_array.push_back("416.29"); int_array.push_back(416);
  intensities.push_back(31);  mzs.push_back(418.232); float_array.push_back(418.232f); string_array.push_back("418.23"); int_array.push_back(418);
  intensities.push_back(29);  mzs.push_back(414.301); float_array.push_back(414.301f); string_array.push_back("414.30"); int_array.push_back(414);
  intensities.push_back(29);  mzs.push_back(412.321); float_array.push_back(412.321f); string_array.push_back("412.32"); int_array.push_back(412);

  for (Size i = 0; i < mzs.size(); ++i)
  {
    p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
    ds.push_back(p);
  }
  ds.sortByIntensity();
  std::vector<double> intensities_copy(intensities);
  std::sort(intensities_copy.begin(),intensities_copy.end());
  MSSpectrum<>::iterator it_ds = ds.begin();
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
  for (Size i = 0; i < mzs.size(); ++i)
  {
    p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
    ds.push_back(p);
  }
  intensities_copy = intensities;
  std::sort(intensities_copy.begin(),intensities_copy.end());

  ds.getFloatDataArrays() = std::vector<MSSpectrum<>::FloatDataArray>(3,float_array);
  ds.getFloatDataArrays()[0].setName("f1");
  ds.getFloatDataArrays()[1].setName("f2");
  ds.getFloatDataArrays()[2].setName("f3");

  ds.getStringDataArrays() = std::vector<MSSpectrum<>::StringDataArray>(2, string_array);
  ds.getStringDataArrays()[0].setName("s1");
  ds.getStringDataArrays()[1].setName("s2");

  ds.getIntegerDataArrays() = std::vector<MSSpectrum<>::IntegerDataArray>(1, int_array);
  ds.getIntegerDataArrays()[0].setName("i1");

  ds.sortByIntensity();

  TEST_STRING_EQUAL(ds.getFloatDataArrays()[0].getName(),"f1")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[1].getName(),"f2")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[2].getName(),"f3")

  TEST_STRING_EQUAL(ds.getStringDataArrays()[0].getName(),"s1")
  TEST_STRING_EQUAL(ds.getStringDataArrays()[1].getName(),"s2")

  TEST_STRING_EQUAL(ds.getIntegerDataArrays()[0].getName(),"i1")
  
  MSSpectrum<>::iterator it1 = ds.begin();
  MSSpectrum<>::FloatDataArray::iterator it2 = ds.getFloatDataArrays()[1].begin();
  MSSpectrum<>::StringDataArray::iterator it3 = ds.getStringDataArrays()[0].begin();
  MSSpectrum<>::IntegerDataArray::iterator it4 = ds.getIntegerDataArrays()[0].begin();
  TOLERANCE_ABSOLUTE(0.0001)
  for(std::vector<double>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
  {
    if(it1 != ds.end() && it2 != ds.getFloatDataArrays()[1].end() && it3 != ds.getStringDataArrays()[0].end() && it4 != ds.getIntegerDataArrays()[0].end())
    {
      //metadataarray values == mz values
      TEST_REAL_SIMILAR(it1->getIntensity(), *it);
      TEST_REAL_SIMILAR(*it2 , it1->getMZ());
      TEST_STRING_EQUAL(*it3 , String::number(it1->getMZ(),2));
      TEST_EQUAL(*it4 , (Int)floor(it1->getMZ()));
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
END_SECTION

START_SECTION((void sortByPosition()))
  MSSpectrum<> ds;
  Peak1D p;
  MSSpectrum<>::FloatDataArray float_array;
  MSSpectrum<>::StringDataArray string_array;
  MSSpectrum<>::IntegerDataArray int_array;
  std::vector<double> mzs, intensities;
  intensities.push_back(56);  mzs.push_back(423.269); float_array.push_back(56);  string_array.push_back("56");  int_array.push_back(56);  
  intensities.push_back(201); mzs.push_back(420.130); float_array.push_back(201); string_array.push_back("201"); int_array.push_back(201); 
  intensities.push_back(31);  mzs.push_back(419.113); float_array.push_back(31);  string_array.push_back("31");  int_array.push_back(31);  
  intensities.push_back(31);  mzs.push_back(418.232); float_array.push_back(31);  string_array.push_back("31");  int_array.push_back(31);  
  intensities.push_back(31);  mzs.push_back(416.293); float_array.push_back(31);  string_array.push_back("31");  int_array.push_back(31);  
  intensities.push_back(37);  mzs.push_back(415.287); float_array.push_back(37);  string_array.push_back("37");  int_array.push_back(37);  
  intensities.push_back(29);  mzs.push_back(414.301); float_array.push_back(29);  string_array.push_back("29");  int_array.push_back(29);  
  intensities.push_back(34);  mzs.push_back(413.800); float_array.push_back(34);  string_array.push_back("34");  int_array.push_back(34);  
  intensities.push_back(60);  mzs.push_back(412.824); float_array.push_back(60);  string_array.push_back("60");  int_array.push_back(60);  
  intensities.push_back(29);  mzs.push_back(412.321); float_array.push_back(29);  string_array.push_back("29");  int_array.push_back(29);  

  for (Size i = 0; i < mzs.size(); ++i)
  {
    p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
    ds.push_back(p);
  }
  ds.sortByPosition();
  MSSpectrum<>::iterator it = ds.begin();
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
  for (Size i = 0; i < mzs.size(); ++i)
  {
    p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
    ds.push_back(p);
  }
  ds.getFloatDataArrays() = std::vector<MSSpectrum<>::FloatDataArray>(3,float_array);
  ds.getFloatDataArrays()[0].setName("f1");
  ds.getFloatDataArrays()[1].setName("f2");
  ds.getFloatDataArrays()[2].setName("f3");
  
  ds.getStringDataArrays() = std::vector<MSSpectrum<>::StringDataArray>(2, string_array);
  ds.getStringDataArrays()[0].setName("s1");
  ds.getStringDataArrays()[1].setName("s2");

  ds.getIntegerDataArrays() = std::vector<MSSpectrum<>::IntegerDataArray>(2, int_array);
  ds.getIntegerDataArrays()[0].setName("i1");
  
  ds.sortByPosition();

  TEST_STRING_EQUAL(ds.getFloatDataArrays()[0].getName(),"f1")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[1].getName(),"f2")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[2].getName(),"f3")

  TEST_STRING_EQUAL(ds.getStringDataArrays()[0].getName(),"s1")
  TEST_STRING_EQUAL(ds.getStringDataArrays()[1].getName(),"s2")
  
  TEST_STRING_EQUAL(ds.getIntegerDataArrays()[0].getName(),"i1")

  MSSpectrum<>::iterator it1 = ds.begin();
  MSSpectrum<>::FloatDataArray::iterator it2 = ds.getFloatDataArrays()[1].begin();
  MSSpectrum<>::StringDataArray::iterator it3 = ds.getStringDataArrays()[0].begin();
  MSSpectrum<>::IntegerDataArray::iterator it4 = ds.getIntegerDataArrays()[0].begin();
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

END_SECTION

START_SECTION(bool isSorted() const)
  //make test dataset
  MSSpectrum<> spec;
  Peak1D p;
  p.setIntensity(1.0);
  p.setMZ(1000.0);
  spec.push_back(p);
  
  p.setIntensity(1.0);
  p.setMZ(1001.0);
  spec.push_back(p);
  
  p.setIntensity(1.0);
  p.setMZ(1002.0);
  spec.push_back(p);
  
  TEST_EQUAL(spec.isSorted(),true)
  
  reverse(spec.begin(), spec.end());
  TEST_EQUAL(spec.isSorted(),false)
END_SECTION

/////////////////////////////////////////////////////////////
// Finding peaks or peak ranges

START_SECTION((Iterator MZEnd(CoordinateType mz)))
  MSSpectrum<> tmp;
  MSSpectrum<>::PeakType rdp;
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

  MSSpectrum<>::Iterator it;

  it = tmp.MZBegin(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZBegin(5.0);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZBegin(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION((Iterator MZBegin(CoordinateType mz)))
  MSSpectrum<> tmp;
  MSSpectrum<>::PeakType rdp;
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

  MSSpectrum<>::Iterator it;

  it = tmp.MZEnd(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZEnd(5.0);
  TEST_EQUAL(it->getPosition()[0],6.0)
  it = tmp.MZEnd(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION((Iterator MZBegin(Iterator begin, CoordinateType mz, Iterator end)))
  MSSpectrum<> tmp;
  MSSpectrum<>::PeakType rdp;
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

  MSSpectrum<>::Iterator it;

  it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
END_SECTION

START_SECTION((ConstIterator MZBegin(ConstIterator begin, CoordinateType mz, ConstIterator end) const))
  MSSpectrum<> tmp;
  MSSpectrum<>::PeakType rdp;
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

  MSSpectrum<>::ConstIterator it;

  it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
END_SECTION

START_SECTION((Iterator MZEnd(Iterator begin, CoordinateType mz, Iterator end)))
  MSSpectrum<> tmp;
  MSSpectrum<>::PeakType rdp;
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

  MSSpectrum<>::Iterator it;

  it = tmp.MZEnd(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZEnd(tmp.begin(), 5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],6.0)
  it = tmp.MZEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
END_SECTION

START_SECTION((ConstIterator MZEnd(ConstIterator begin, CoordinateType mz, ConstIterator end) const))
  MSSpectrum<> tmp;
  MSSpectrum<>::PeakType rdp;
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

  MSSpectrum<>::ConstIterator it;

  it = tmp.MZEnd(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZEnd(tmp.begin(), 5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],6.0)
  it = tmp.MZEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
END_SECTION

START_SECTION((ConstIterator MZEnd(CoordinateType mz) const))
  MSSpectrum<> tmp;
  MSSpectrum<>::PeakType rdp;
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

  MSSpectrum<>::ConstIterator it;

  it = tmp.MZBegin(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZBegin(5.0);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZBegin(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION((ConstIterator MZBegin(CoordinateType mz) const))
  MSSpectrum<> tmp;
  MSSpectrum<>::PeakType rdp;
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

  MSSpectrum<>::ConstIterator it;

  it = tmp.MZEnd(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZEnd(5.0);
  TEST_EQUAL(it->getPosition()[0],6.0)
  it = tmp.MZEnd(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION((Size findNearest(CoordinateType mz) const))
  MSSpectrum<> tmp;
  Peak1D p;
  p.setIntensity(29.0f); p.setMZ(412.321); tmp.push_back(p); //0
  p.setIntensity(60.0f); p.setMZ(412.824); tmp.push_back(p); //1
  p.setIntensity(34.0f); p.setMZ(413.8); tmp.push_back(p); //2
  p.setIntensity(29.0f); p.setMZ(414.301); tmp.push_back(p); //3
  p.setIntensity(37.0f); p.setMZ(415.287); tmp.push_back(p); //4
  p.setIntensity(31.0f); p.setMZ(416.293); tmp.push_back(p); //5
  p.setIntensity(31.0f); p.setMZ(418.232); tmp.push_back(p); //6
  p.setIntensity(31.0f); p.setMZ(419.113); tmp.push_back(p); //7
  p.setIntensity(201.0f); p.setMZ(420.13); tmp.push_back(p); //8
  p.setIntensity(56.0f); p.setMZ(423.269); tmp.push_back(p); //9
  p.setIntensity(34.0f); p.setMZ(426.292); tmp.push_back(p); //10
  p.setIntensity(82.0f); p.setMZ(427.28); tmp.push_back(p); //11
  p.setIntensity(87.0f); p.setMZ(428.322); tmp.push_back(p); //12
  p.setIntensity(30.0f); p.setMZ(430.269); tmp.push_back(p); //13
  p.setIntensity(29.0f); p.setMZ(431.246); tmp.push_back(p); //14
  p.setIntensity(42.0f); p.setMZ(432.289); tmp.push_back(p); //15
  p.setIntensity(32.0f); p.setMZ(436.161); tmp.push_back(p); //16
  p.setIntensity(54.0f); p.setMZ(437.219); tmp.push_back(p); //17
  p.setIntensity(40.0f); p.setMZ(439.186); tmp.push_back(p); //18
  p.setIntensity(40); p.setMZ(440.27); tmp.push_back(p); //19
  p.setIntensity(23.0f); p.setMZ(441.224); tmp.push_back(p); //20

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
  MSSpectrum<> tmp2;
  TEST_PRECONDITION_VIOLATED(tmp2.findNearest(427.3));
END_SECTION

START_SECTION((Size findNearest(CoordinateType mz, CoordinateType tolerance) const))
  MSSpectrum<> tmp;
  Peak1D p;
  p.setIntensity(29.0f); p.setMZ(412.321); tmp.push_back(p); //0
  p.setIntensity(60.0f); p.setMZ(412.824); tmp.push_back(p); //1
  p.setIntensity(34.0f); p.setMZ(413.8); tmp.push_back(p); //2
  p.setIntensity(29.0f); p.setMZ(414.301); tmp.push_back(p); //3
  p.setIntensity(37.0f); p.setMZ(415.287); tmp.push_back(p); //4
  p.setIntensity(31.0f); p.setMZ(416.293); tmp.push_back(p); //5
  p.setIntensity(31.0f); p.setMZ(418.232); tmp.push_back(p); //6
  p.setIntensity(31.0f); p.setMZ(419.113); tmp.push_back(p); //7
  p.setIntensity(201.0f); p.setMZ(420.13); tmp.push_back(p); //8
  p.setIntensity(56.0f); p.setMZ(423.269); tmp.push_back(p); //9
  p.setIntensity(34.0f); p.setMZ(426.292); tmp.push_back(p); //10
  p.setIntensity(82.0f); p.setMZ(427.28); tmp.push_back(p); //11
  p.setIntensity(87.0f); p.setMZ(428.322); tmp.push_back(p); //12
  p.setIntensity(30.0f); p.setMZ(430.269); tmp.push_back(p); //13
  p.setIntensity(29.0f); p.setMZ(431.246); tmp.push_back(p); //14
  p.setIntensity(42.0f); p.setMZ(432.289); tmp.push_back(p); //15
  p.setIntensity(32.0f); p.setMZ(436.161); tmp.push_back(p); //16
  p.setIntensity(54.0f); p.setMZ(437.219); tmp.push_back(p); //17
  p.setIntensity(40.0f); p.setMZ(439.186); tmp.push_back(p); //18
  p.setIntensity(40); p.setMZ(440.27); tmp.push_back(p); //19
  p.setIntensity(23.0f); p.setMZ(441.224); tmp.push_back(p); //20

  //test outside mass range
  TEST_EQUAL(tmp.findNearest(400.0, 1.0), -1);
  TEST_EQUAL(tmp.findNearest(500.0, 1.0), -1);

  //test mass range borders
  TEST_EQUAL(tmp.findNearest(412.4, 0.01), -1);
  TEST_EQUAL(tmp.findNearest(412.4, 0.1), 0);
  TEST_EQUAL(tmp.findNearest(441.3, 0.01),-1);
  TEST_EQUAL(tmp.findNearest(441.3, 0.1), 20);

  //test inside scan
  TEST_EQUAL(tmp.findNearest(426.29, 0.1), 10);
  TEST_EQUAL(tmp.findNearest(426.3, 0.1), 10);
  TEST_EQUAL(tmp.findNearest(427.2, 0.1), 11);
  TEST_EQUAL(tmp.findNearest(427.3, 0.1), 11);
  TEST_EQUAL(tmp.findNearest(427.3, 0.001), -1);

  //empty spectrum
  MSSpectrum<> tmp2;
  TEST_EQUAL(tmp2.findNearest(427.3, 1.0, 1.0), -1);
END_SECTION
START_SECTION((Size findNearest(CoordinateType mz, CoordinateType left_tolerance, CoordinateType right_tolerance) const))
  MSSpectrum<> tmp;
  Peak1D p;
  p.setIntensity(29.0f); p.setMZ(412.321); tmp.push_back(p); //0
  p.setIntensity(60.0f); p.setMZ(412.824); tmp.push_back(p); //1
  p.setIntensity(34.0f); p.setMZ(413.8); tmp.push_back(p); //2
  p.setIntensity(29.0f); p.setMZ(414.301); tmp.push_back(p); //3
  p.setIntensity(37.0f); p.setMZ(415.287); tmp.push_back(p); //4
  p.setIntensity(31.0f); p.setMZ(416.293); tmp.push_back(p); //5
  p.setIntensity(31.0f); p.setMZ(418.232); tmp.push_back(p); //6
  p.setIntensity(31.0f); p.setMZ(419.113); tmp.push_back(p); //7
  p.setIntensity(201.0f); p.setMZ(420.13); tmp.push_back(p); //8
  p.setIntensity(56.0f); p.setMZ(423.269); tmp.push_back(p); //9
  p.setIntensity(34.0f); p.setMZ(426.292); tmp.push_back(p); //10
  p.setIntensity(82.0f); p.setMZ(427.28); tmp.push_back(p); //11
  p.setIntensity(87.0f); p.setMZ(428.322); tmp.push_back(p); //12
  p.setIntensity(30.0f); p.setMZ(430.269); tmp.push_back(p); //13
  p.setIntensity(29.0f); p.setMZ(431.246); tmp.push_back(p); //14
  p.setIntensity(42.0f); p.setMZ(432.289); tmp.push_back(p); //15
  p.setIntensity(32.0f); p.setMZ(436.161); tmp.push_back(p); //16
  p.setIntensity(54.0f); p.setMZ(437.219); tmp.push_back(p); //17
  p.setIntensity(40.0f); p.setMZ(439.186); tmp.push_back(p); //18
  p.setIntensity(40); p.setMZ(440.27); tmp.push_back(p); //19
  p.setIntensity(23.0f); p.setMZ(441.224); tmp.push_back(p); //20

  //test outside mass range
  TEST_EQUAL(tmp.findNearest(400.0, 1.0, 1.0), -1);
  TEST_EQUAL(tmp.findNearest(500.0, 1.0, 1.0), -1);

  //test mass range borders
  TEST_EQUAL(tmp.findNearest(412.4, 0.01, 0.01), -1);
  TEST_EQUAL(tmp.findNearest(412.4, 0.1, 0.1), 0);
  TEST_EQUAL(tmp.findNearest(441.3, 0.01, 0.01),-1);
  TEST_EQUAL(tmp.findNearest(441.3, 0.1, 0.1), 20);

  //test inside scan
  TEST_EQUAL(tmp.findNearest(426.29, 0.1, 0.1), 10);
  TEST_EQUAL(tmp.findNearest(426.3, 0.1, 0.1), 10);
  TEST_EQUAL(tmp.findNearest(427.2, 0.1, 0.1), 11);
  TEST_EQUAL(tmp.findNearest(427.3, 0.1, 0.1), 11);
  TEST_EQUAL(tmp.findNearest(427.3, 0.001, 0.001), -1);

  TEST_EQUAL(tmp.findNearest(427.3, 0.1, 0.001), 11);
  TEST_EQUAL(tmp.findNearest(427.3, 0.001, 1.01), -1);
  TEST_EQUAL(tmp.findNearest(427.3, 0.001, 1.1), 12); 

  //empty spectrum
  MSSpectrum<> tmp2;
  TEST_EQUAL(tmp2.findNearest(427.3, 1.0, 1.0), -1);
END_SECTION

START_SECTION(void clear(bool clear_meta_data))
  MSSpectrum<> edit;
  edit.getInstrumentSettings().getScanWindows().resize(1);
  edit.resize(1);
  edit.setMetaValue("label",String("bla"));
  edit.setRT(5);
  edit.setDriftTime(6);
  edit.setMSLevel(5);
  edit.getFloatDataArrays().resize(5);
  edit.getIntegerDataArrays().resize(5);
  edit.getStringDataArrays().resize(5);

  edit.clear(false);
  TEST_EQUAL(edit.size(),0)
  TEST_EQUAL(edit==MSSpectrum<>(),false)

  edit.clear(true);
  TEST_EQUAL(edit==MSSpectrum<>(),true)
END_SECTION

START_SECTION(([MSSpectrum::RTLess] bool operator()(const MSSpectrum &a, const MSSpectrum &b) const))
  vector< MSSpectrum<> > v;

  MSSpectrum<> sp1;
  sp1.setRT(3.0f);
  v.push_back(sp1);

  MSSpectrum<> sp2;
  sp2.setRT(2.0f);
  v.push_back(sp2);

  MSSpectrum<> sp3;
  sp3.setRT(1.0f);
  v.push_back(sp3);

  std::sort(v.begin(),v.end(), MSSpectrum<>::RTLess());

  TEST_REAL_SIMILAR(v[0].getRT(), 1.0);
  TEST_REAL_SIMILAR(v[1].getRT(), 2.0);
  TEST_REAL_SIMILAR(v[2].getRT(), 3.0);

  ///
  MSSpectrum<> s1;
  s1.setRT(0.451);

  MSSpectrum<> s2;
  s2.setRT(0.5);

  TEST_EQUAL(MSSpectrum<>::RTLess()(s1,s2), true);
  TEST_EQUAL(MSSpectrum<>::RTLess()(s2,s1), false);
  TEST_EQUAL(MSSpectrum<>::RTLess()(s2,s2), false);
END_SECTION

START_SECTION(([EXTRA] std::ostream& operator << (std::ostream& os, const MSSpectrum<PeakT>& spec)))
{
  MSSpectrum<> spec;
  Peak1D p;
  p.setIntensity(29.0f); p.setMZ(412.321); spec.push_back(p); //0
  p.setIntensity(60.0f); p.setMZ(412.824); spec.push_back(p); //1
  p.setIntensity(34.0f); p.setMZ(413.8);   spec.push_back(p); //2
  p.setIntensity(29.0f); p.setMZ(414.301); spec.push_back(p); //3
  p.setIntensity(37.0f); p.setMZ(415.287); spec.push_back(p); //4
  p.setIntensity(31.0f); p.setMZ(416.293); spec.push_back(p); //5
  p.setIntensity(31.0f); p.setMZ(418.232); spec.push_back(p); //6
  p.setIntensity(31.0f); p.setMZ(419.113); spec.push_back(p); //7
  p.setIntensity(201.0f); p.setMZ(420.13); spec.push_back(p); //8
  p.setIntensity(56.0f); p.setMZ(423.269); spec.push_back(p); //9
  p.setIntensity(34.0f); p.setMZ(426.292); spec.push_back(p); //10

  spec.getInstrumentSettings().getScanWindows().resize(1);
  spec.setMetaValue("label",5.0);
  spec.setMSLevel(17);
  spec.setRT(7.0);
  spec.setName("bla");

  ostringstream test_stream;
  test_stream << spec;

  TEST_EQUAL(test_stream.str(), "-- MSSPECTRUM BEGIN --\n"
                                "-- SPECTRUMSETTINGS BEGIN --\n"
                                "-- SPECTRUMSETTINGS END --\n"
                                "POS: 412.321 INT: 29\n"
                                "POS: 412.824 INT: 60\n"
                                "POS: 413.8 INT: 34\n"
                                "POS: 414.301 INT: 29\n"
                                "POS: 415.287 INT: 37\n"
                                "POS: 416.293 INT: 31\n"
                                "POS: 418.232 INT: 31\n"
                                "POS: 419.113 INT: 31\n"
                                "POS: 420.13 INT: 201\n"
                                "POS: 423.269 INT: 56\n"
                                "POS: 426.292 INT: 34\n"
                                "-- MSSPECTRUM END --\n")

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



