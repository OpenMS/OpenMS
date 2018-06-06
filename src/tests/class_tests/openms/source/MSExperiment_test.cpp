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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Peak2D.h>

///////////////////////////

START_TEST(MSExperiment, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PeakMap* ptr = nullptr;
PeakMap* nullPointer = nullptr;
START_SECTION((MSExperiment()))
{
	ptr = new PeakMap;
  TEST_NOT_EQUAL(ptr, nullPointer);
}
END_SECTION

START_SECTION(([EXTRA]~MSExperiment()))
{
	delete ptr;
}
END_SECTION

START_SECTION((MSExperiment(const MSExperiment& source)))
{
  PeakMap tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setFirstName("Name");
  tmp.resize(1);

  PeakMap tmp2(tmp);
  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL(tmp2.getContacts()[0].getFirstName(),"Name");
  TEST_EQUAL(tmp2.size(),1);
}
END_SECTION

START_SECTION((MSExperiment& operator= (const MSExperiment& source)))
  PeakMap tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setFirstName("Name");
  tmp.resize(1);
  Peak1D p;
  p.setMZ(5.0);
  tmp[0].push_back(p);
  p.setMZ(10.0);
  tmp[0].push_back(p);
  tmp.updateRanges();

  PeakMap tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL(tmp2.getContacts()[0].getFirstName(),"Name");
  TEST_EQUAL(tmp2.size(),1);
	TEST_REAL_SIMILAR(tmp2.getMinMZ(),5.0);
	TEST_REAL_SIMILAR(tmp2.getMaxMZ(),10.0);

  tmp2 = PeakMap();
  TEST_EQUAL(tmp2.getContacts().size(),0);
  TEST_EQUAL(tmp2.size(),0);
END_SECTION

START_SECTION((bool operator== (const MSExperiment& rhs) const))
  PeakMap edit,empty;

	TEST_EQUAL(edit==empty, true);

	edit.getContacts().resize(1);
  TEST_EQUAL(edit==empty, false);

	edit = empty;
  edit.resize(1);
	TEST_EQUAL(edit==empty, false);
END_SECTION

START_SECTION((bool operator!= (const MSExperiment& rhs) const))
  PeakMap edit,empty;

	TEST_EQUAL(edit!=empty, false);

	edit.getContacts().resize(1);
  TEST_EQUAL(edit!=empty, true);

	edit = empty;
  edit.resize(1);
	TEST_EQUAL(edit!=empty, true);
END_SECTION

START_SECTION((template<class Container> void get2DData(Container& cont) const))
	PeakMap exp;
	PeakMap::SpectrumType spec;
	PeakMap::PeakType peak;

	// first spectrum (MS)
	spec.setRT(11.1);
	spec.setMSLevel(1);
	peak.getPosition()[0] = 5;
	peak.setIntensity(47.11f);
	spec.push_back(peak);
	peak.getPosition()[0] = 10;
	peak.setIntensity(48.11f);
	spec.push_back(peak);
	peak.getPosition()[0] = 15;
	spec.push_back(peak);
	exp.addSpectrum(spec);

	// second spectrum (MS/MS)
	spec.clear(true);
	spec.setRT(11.5);
	spec.setMSLevel(2);
	peak.getPosition()[0] = 6;
	spec.push_back(peak);
	peak.getPosition()[0] = 11;
	spec.push_back(peak);
	exp.addSpectrum(spec);

	// third spectrum (MS)
	spec.clear(true);
	spec.setRT(12.2);
	spec.setMSLevel(1);
	peak.getPosition()[0] = 20;
	spec.push_back(peak);
	peak.getPosition()[0] = 25;
	spec.push_back(peak);
	exp.addSpectrum(spec);

	// forth spectrum (MS/MS)
	spec.clear(true);
	spec.setRT(12.5);
	spec.setMSLevel(2);
	peak.getPosition()[0] = 21;
	spec.push_back(peak);
	peak.getPosition()[0] = 26;
	spec.push_back(peak);
	peak.getPosition()[0] = 31;
	spec.push_back(peak);
	exp.addSpectrum(spec);

	//Convert
	std::vector<Peak2D> a;
	exp.get2DData(a);

	//Tests
	TEST_EQUAL(a.size(),5);
	TEST_REAL_SIMILAR(a[0].getRT(),11.1);
	TEST_REAL_SIMILAR(a[0].getMZ(),5);
	TEST_REAL_SIMILAR(a[0].getIntensity(),47.11);
	TEST_REAL_SIMILAR(a[1].getRT(),11.1);
	TEST_REAL_SIMILAR(a[1].getMZ(),10);
	TEST_REAL_SIMILAR(a[1].getIntensity(),48.11);
	TEST_REAL_SIMILAR(a[2].getRT(),11.1);
	TEST_REAL_SIMILAR(a[2].getMZ(),15);
	TEST_REAL_SIMILAR(a[3].getRT(),12.2);
	TEST_REAL_SIMILAR(a[3].getMZ(),20);
	TEST_REAL_SIMILAR(a[4].getRT(),12.2);
	TEST_REAL_SIMILAR(a[4].getMZ(),25);

	//Convert
	std::vector<Peak2D> list;
	exp.get2DData(list);

	//Tests
	TEST_EQUAL(list.size(),5);
	std::vector<Peak2D>::const_iterator it = list.begin();
	TEST_REAL_SIMILAR(it->getRT(),11.1);
	TEST_REAL_SIMILAR(it->getMZ(),5);
	TEST_REAL_SIMILAR(it->getIntensity(),47.11);
	++it;
	TEST_REAL_SIMILAR(it->getRT(),11.1);
	TEST_REAL_SIMILAR(it->getMZ(),10);
	TEST_REAL_SIMILAR(it->getIntensity(),48.11);
	++it;
	TEST_REAL_SIMILAR(it->getRT(),11.1);
	TEST_REAL_SIMILAR(it->getMZ(),15);
	++it;
	TEST_REAL_SIMILAR(it->getRT(),12.2);
	TEST_REAL_SIMILAR(it->getMZ(),20);
	++it;
	TEST_REAL_SIMILAR(it->getRT(),12.2);
	TEST_REAL_SIMILAR(it->getMZ(),25);
END_SECTION

START_SECTION((template <class Container> void set2DData(const Container& cont, const StringList& store_metadata_names = StringList())))
   NOT_TESTABLE // tested below
END_SECTION

START_SECTION((template <bool add_mass_traces, class Container> void set2DData(const Container& cont, const StringList& store_metadata_names = StringList())))
	PeakMap exp;

	// create sample data
	std::vector<Peak2D> input;

	Peak2D p1(Peak2D::PositionType(2.0, 3.0), 1.0);
	input.push_back(p1);

	Peak2D p2(Peak2D::PositionType(5.0, 6.0), 4.0);
	input.push_back(p2);

	Peak2D p3(Peak2D::PositionType(8.5, 9.5), 7.5);
	input.push_back(p3);

	exp.set2DData(input);

	// retrieve data again and check for changes
	std::vector<Peak2D> output;

	exp.get2DData(output);
	TEST_EQUAL(output==input,true);
  
  //////////////////////////////////////////////////////////////////////////
  // test if meta values are added as floatDataArrays in MSSpectra
  std::vector<RichPeak2D> inputr;

  RichPeak2D pr1(RichPeak2D::PositionType(2.0, 3.0), 1.0);
  pr1.setMetaValue("meta1", 111.1);
  inputr.push_back(pr1);
  RichPeak2D pr2(RichPeak2D::PositionType(5.0, 6.0), 4.0);
  inputr.push_back(pr2);
  RichPeak2D pr3(RichPeak2D::PositionType(8.5, 9.5), 7.5);
  pr3.setMetaValue("meta3", 333.3);
  inputr.push_back(pr3);
  
  // create float data arrays for these two meta values (missing values in data will be set to NaN)
  exp.set2DData(inputr, ListUtils::create<String>("meta1,meta3"));
  TEST_EQUAL(exp.getNrSpectra(), 3);
  // retrieve data again and check for changes
  std::vector<Peak2D> outputr;
  exp.get2DData(outputr);
  TEST_EQUAL(outputr==input, true); // we compare to non-meta output, since floatdata is not converted back to metavalues
  // check for meta data
  TEST_EQUAL(exp[0].getFloatDataArrays().size(), 2);
  TEST_EQUAL(exp[0].getFloatDataArrays()[0][0], 111.1);
  TEST_EQUAL(exp[1].getFloatDataArrays().size(), 2); // present but all NaN
  TEST_EQUAL(exp[2].getFloatDataArrays().size(), 2);
  TEST_EQUAL(exp[2].getFloatDataArrays()[1][0], 333.3);

  ///////////////////////////////////////
  // test adding of mass traces
	FeatureMap fm, fm2, fm_out;
  Feature f1;
	f1.setIntensity(7.5f);
	f1.setRT(8.5);
	f1.setMZ(9.5);
  Feature f2;
	f2.setIntensity(17.5f);
	f2.setRT(18.5);
	f2.setMZ(19.5);
  fm.push_back(f1);
  fm.push_back(f2);  
  fm2 = fm; // copy without meta values (get2DData will not have them)
  fm.back().setMetaValue("num_of_masstraces", 2);
  fm.back().setMetaValue("masstrace_intensity_0", 11.0f);
  fm.back().setMetaValue("masstrace_intensity_1", 12.0f);
  fm.back().setCharge(2);
  exp.set2DData<true>(fm);
	exp.get2DData(fm_out);
  Feature f2_x = f2;
  f2_x.setIntensity(11.0f);
  f2_x.setMZ(f2_x.getMZ() + Constants::C13C12_MASSDIFF_U/2*0);
  fm2.back() = (f2_x);   // replace
  f2_x.setIntensity(12.0f);
  f2_x.setMZ(f2_x.getMZ() + Constants::C13C12_MASSDIFF_U/2*1);
  fm2.push_back(f2_x);  // add +1Th trace
	TEST_EQUAL(fm_out.size(), fm2.size());
  //std::cout << fm_out[1] << "\n\n" << fm2[1] << std::endl;
  //std::cout << fm_out.back() << "\n\n" << fm2.back() << std::endl;
	TEST_EQUAL(fm_out==fm2,true);

	// test precondition (input sorted by RT)
	input.push_back(p1);
	TEST_PRECONDITION_VIOLATED(exp.set2DData(input));

END_SECTION

START_SECTION(([EXTRA] PeakMap()))
	PeakMap tmp;
	tmp.resize(1);
	tmp[0].resize(1);
	tmp[0][0].getPosition()[0] = 47.11;
	TEST_REAL_SIMILAR(tmp[0][0].getPosition()[0],47.11)
END_SECTION

START_SECTION((CoordinateType getMinMZ() const))
	PeakMap tmp;
	TEST_REAL_SIMILAR(tmp.getMinMZ(),numeric_limits<DPosition<2>::CoordinateType>::max())
END_SECTION

START_SECTION((CoordinateType getMaxMZ() const))
	PeakMap tmp;
	TEST_REAL_SIMILAR(tmp.getMaxMZ(),-numeric_limits<DPosition<2>::CoordinateType>::max())
END_SECTION

START_SECTION((CoordinateType getMinRT() const))
	PeakMap tmp;
	TEST_REAL_SIMILAR(tmp.getMinRT(),numeric_limits<DPosition<2>::CoordinateType>::max())
END_SECTION

START_SECTION((CoordinateType getMaxRT() const))
	PeakMap tmp;
	TEST_REAL_SIMILAR(tmp.getMaxRT(),-numeric_limits<DPosition<2>::CoordinateType>::max())
END_SECTION

START_SECTION((const std::vector<UInt>& getMSLevels() const))
	PeakMap tmp;
	TEST_EQUAL(tmp.getMSLevels().size(),0)
END_SECTION

START_SECTION((UInt64 getSize() const ))
	PeakMap tmp;
	TEST_EQUAL(tmp.getSize(),0)
END_SECTION

START_SECTION((const AreaType& getDataRange() const))
	PeakMap tmp;
	TEST_REAL_SIMILAR(tmp.getDataRange().minPosition()[1],numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_SIMILAR(tmp.getDataRange().maxPosition()[1],-numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_SIMILAR(tmp.getDataRange().minPosition()[0],numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_SIMILAR(tmp.getDataRange().maxPosition()[0],-numeric_limits<DPosition<2>::CoordinateType>::max())
END_SECTION

START_SECTION((virtual void updateRanges()))
	PeakMap tmp;
	MSSpectrum s;
	Peak1D p;

	s.setMSLevel(1);
	s.setRT(30.0);
	p.getPosition()[0] = 5.0;
	p.setIntensity(-5.0f);
	s.push_back(p);
	tmp.addSpectrum(s);

	s.clear(true);
	s.setMSLevel(1);
	s.setRT(40.0);
	p.getPosition()[0] = 7.0;
	p.setIntensity(-7.0f);
	s.push_back(p);
	tmp.addSpectrum(s);

	s.clear(true);
	s.setMSLevel(3);
	s.setRT(45.0);
	p.getPosition()[0] = 9.0;
	p.setIntensity(-10.0f);
	s.push_back(p);
	tmp.addSpectrum(s);

	s.clear(true);
	s.setMSLevel(3);
	s.setRT(50.0);
	p.getPosition()[0] = 10.0;
	p.setIntensity(-9.0f);
	s.push_back(p);
	tmp.addSpectrum(s);

	tmp.updateRanges();
	tmp.updateRanges(); //second time to check the initialization

	TEST_REAL_SIMILAR(tmp.getMinMZ(),5.0)
	TEST_REAL_SIMILAR(tmp.getMaxMZ(),10.0)
	TEST_REAL_SIMILAR(tmp.getMinInt(),-10.0)
	TEST_REAL_SIMILAR(tmp.getMaxInt(),-5.0)
	TEST_REAL_SIMILAR(tmp.getMinRT(),30.0)
	TEST_REAL_SIMILAR(tmp.getMaxRT(),50.0)
	TEST_EQUAL(tmp.getMSLevels().size(),2)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getMSLevels()[1],3)
	TEST_EQUAL(tmp.getSize(),4)
	tmp.updateRanges();
	TEST_REAL_SIMILAR(tmp.getMinMZ(),5.0)
	TEST_REAL_SIMILAR(tmp.getMaxMZ(),10.0)
	TEST_REAL_SIMILAR(tmp.getMinInt(),-10.0)
	TEST_REAL_SIMILAR(tmp.getMaxInt(),-5.0)
	TEST_REAL_SIMILAR(tmp.getMinRT(),30.0)
	TEST_REAL_SIMILAR(tmp.getMaxRT(),50.0)

	TEST_REAL_SIMILAR(tmp.getDataRange().minPosition()[1],5.0)
	TEST_REAL_SIMILAR(tmp.getDataRange().maxPosition()[1],10.0)
	TEST_REAL_SIMILAR(tmp.getDataRange().minPosition()[0],30.0)
	TEST_REAL_SIMILAR(tmp.getDataRange().maxPosition()[0],50.0)

	TEST_EQUAL(tmp.getMSLevels().size(),2)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getMSLevels()[1],3)

	TEST_EQUAL(tmp.getSize(),4)

	//Update for MS level 1

	tmp.updateRanges(1);
	tmp.updateRanges(1);
	TEST_REAL_SIMILAR(tmp.getMinMZ(),5.0)
	TEST_REAL_SIMILAR(tmp.getMaxMZ(),7.0)
	TEST_REAL_SIMILAR(tmp.getMinInt(),-7.0)
	TEST_REAL_SIMILAR(tmp.getMaxInt(),-5.0)
	TEST_REAL_SIMILAR(tmp.getMinRT(),30.0)
	TEST_REAL_SIMILAR(tmp.getMaxRT(),40.0)
	TEST_EQUAL(tmp.getMSLevels().size(),1)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getSize(),2)
	tmp.updateRanges(1);
	TEST_REAL_SIMILAR(tmp.getMinMZ(),5.0)
	TEST_REAL_SIMILAR(tmp.getMaxMZ(),7.0)
	TEST_REAL_SIMILAR(tmp.getMinInt(),-7.0)
	TEST_REAL_SIMILAR(tmp.getMaxInt(),-5.0)
	TEST_REAL_SIMILAR(tmp.getMinRT(),30.0)
	TEST_REAL_SIMILAR(tmp.getMaxRT(),40.0)
	TEST_EQUAL(tmp.getMSLevels().size(),1)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getSize(),2)

	//test with only one peak

	PeakMap tmp2;
	MSSpectrum s2;
	Peak1D p2;

	s2.setRT(30.0);
	p2.getPosition()[0] = 5.0;
	p2.setIntensity(-5.0f);
	s2.push_back(p2);
	tmp2.addSpectrum(s2);

	tmp2.updateRanges();
	TEST_REAL_SIMILAR(tmp2.getMinMZ(),5.0)
	TEST_REAL_SIMILAR(tmp2.getMaxMZ(),5.0)
	TEST_REAL_SIMILAR(tmp2.getMinInt(),-5.0)
	TEST_REAL_SIMILAR(tmp2.getMaxInt(),-5.0)
	TEST_REAL_SIMILAR(tmp2.getMinRT(),30.0)
	TEST_REAL_SIMILAR(tmp2.getMaxRT(),30.0)

	tmp2.updateRanges(1);
	TEST_REAL_SIMILAR(tmp2.getMinMZ(),5.0)
	TEST_REAL_SIMILAR(tmp2.getMaxMZ(),5.0)
	TEST_REAL_SIMILAR(tmp2.getMinInt(),-5.0)
	TEST_REAL_SIMILAR(tmp2.getMaxInt(),-5.0)
	TEST_REAL_SIMILAR(tmp2.getMinRT(),30.0)
	TEST_REAL_SIMILAR(tmp2.getMaxRT(),30.0)

END_SECTION

START_SECTION((void updateRanges(Int ms_level)))
	PeakMap tmp;
	MSSpectrum s;
	Peak1D p;

	s.setMSLevel(1);
	s.setRT(30.0);
	p.getPosition()[0] = 5.0;
	p.setIntensity(-5.0f);
	s.push_back(p);
	tmp.addSpectrum(s);

	s.clear(true);
	s.setMSLevel(1);
	s.setRT(40.0);
	p.getPosition()[0] = 7.0;
	p.setIntensity(-7.0f);
	s.push_back(p);
	tmp.addSpectrum(s);

	s.clear(true);
	s.setMSLevel(3);
	s.setRT(45.0);
	p.getPosition()[0] = 9.0;
	p.setIntensity(-10.0f);
	s.push_back(p);
	tmp.addSpectrum(s);

	s.clear(true);
	s.setMSLevel(3);
	s.setRT(50.0);
	p.getPosition()[0] = 10.0;
	p.setIntensity(-9.0f);
	s.push_back(p);
	tmp.addSpectrum(s);

	//Update for MS level 1

	tmp.updateRanges(1);
	tmp.updateRanges(1);
	TEST_REAL_SIMILAR(tmp.getMinMZ(),5.0)
	TEST_REAL_SIMILAR(tmp.getMaxMZ(),7.0)
	TEST_REAL_SIMILAR(tmp.getMinInt(),-7.0)
	TEST_REAL_SIMILAR(tmp.getMaxInt(),-5.0)
	TEST_REAL_SIMILAR(tmp.getMinRT(),30.0)
	TEST_REAL_SIMILAR(tmp.getMaxRT(),40.0)
	TEST_EQUAL(tmp.getMSLevels().size(),1)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getSize(),2)
	tmp.updateRanges(1);
	TEST_REAL_SIMILAR(tmp.getMinMZ(),5.0)
	TEST_REAL_SIMILAR(tmp.getMaxMZ(),7.0)
	TEST_REAL_SIMILAR(tmp.getMinInt(),-7.0)
	TEST_REAL_SIMILAR(tmp.getMaxInt(),-5.0)
	TEST_REAL_SIMILAR(tmp.getMinRT(),30.0)
	TEST_REAL_SIMILAR(tmp.getMaxRT(),40.0)
	TEST_EQUAL(tmp.getMSLevels().size(),1)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getSize(),2)

	//test with only one peak

	PeakMap tmp2;
	MSSpectrum s2;
	Peak1D p2;

	s2.setRT(30.0);
	p2.getPosition()[0] = 5.0;
	p2.setIntensity(-5.0f);
	s2.push_back(p2);
	tmp2.addSpectrum(s2);

	tmp2.updateRanges(1);
	TEST_REAL_SIMILAR(tmp2.getMinMZ(),5.0)
	TEST_REAL_SIMILAR(tmp2.getMaxMZ(),5.0)
	TEST_REAL_SIMILAR(tmp2.getMinInt(),-5.0)
	TEST_REAL_SIMILAR(tmp2.getMaxInt(),-5.0)
	TEST_REAL_SIMILAR(tmp2.getMinRT(),30.0)
	TEST_REAL_SIMILAR(tmp2.getMaxRT(),30.0)

END_SECTION

START_SECTION((ConstAreaIterator areaEndConst() const))
NOT_TESTABLE
END_SECTION

START_SECTION((ConstAreaIterator areaBeginConst(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz) const))
	std::vector< Peak2D> plist;

	Peak2D p1;
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 2.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 3.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 2.0;
	p1.getPosition()[1] = 10.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 2.0;
	p1.getPosition()[1] = 11.0;
	plist.push_back(p1);

	PeakMap exp;
	exp.set2DData(plist);

	PeakMap::ConstAreaIterator it = exp.areaBeginConst(0,15,0,15);

	TEST_EQUAL(it->getPosition()[0],2.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],3.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],10.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],11.0);
	it++;
	TEST_EQUAL(it==exp.areaEndConst(),true);

	TEST_PRECONDITION_VIOLATED(exp.areaBeginConst(15,0,0,15));
	TEST_PRECONDITION_VIOLATED(exp.areaBeginConst(0,15,15,0));
	TEST_PRECONDITION_VIOLATED(exp.areaBeginConst(15,0,15,0));

END_SECTION

START_SECTION((AreaIterator areaEnd()))
NOT_TESTABLE
END_SECTION

START_SECTION((AreaIterator areaBegin(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz)))
	std::vector< Peak2D> plist;

	Peak2D p1;
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 2.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 3.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 2.0;
	p1.getPosition()[1] = 10.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 2.0;
	p1.getPosition()[1] = 11.0;
	plist.push_back(p1);

	PeakMap exp;
	exp.set2DData(plist);

	PeakMap::AreaIterator it = exp.areaBegin(0,15,0,15);

	TEST_EQUAL(it->getPosition()[0],2.0);
	it->getPosition()[0] = 4711.0;
	TEST_EQUAL(it->getPosition()[0],4711.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],3.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],10.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],11.0);
	it++;
	TEST_EQUAL(it==exp.areaEnd(),true);

	TEST_PRECONDITION_VIOLATED(exp.areaBegin(15,0,0,15));
	TEST_PRECONDITION_VIOLATED(exp.areaBegin(0,15,15,0));
	TEST_PRECONDITION_VIOLATED(exp.areaBegin(15,0,15,0));

END_SECTION

START_SECTION((Iterator RTBegin(CoordinateType rt)))
	PeakMap tmp;
	MSSpectrum s;
	Peak1D p;

	s.setRT(30.0);
	tmp.addSpectrum(s);
	s.setRT(40.0);
	tmp.addSpectrum(s);
	s.setRT(45.0);
	tmp.addSpectrum(s);
	s.setRT(50.0);
	tmp.addSpectrum(s);

	PeakMap::Iterator it;

	it = tmp.RTBegin(20.0);
	TEST_REAL_SIMILAR(it->getRT(),30.0)
	it = tmp.RTBegin(30.0);
	TEST_REAL_SIMILAR(it->getRT(),30.0)
	it = tmp.RTBegin(31.0);
	TEST_REAL_SIMILAR(it->getRT(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true)
END_SECTION

START_SECTION((Iterator RTEnd(CoordinateType rt)))
	PeakMap tmp;
	MSSpectrum s;
	Peak1D p;

	s.setRT(30.0);
	tmp.addSpectrum(s);
	s.setRT(40.0);
	tmp.addSpectrum(s);
	s.setRT(45.0);
	tmp.addSpectrum(s);
	s.setRT(50.0);
	tmp.addSpectrum(s);

	PeakMap::Iterator it;

	it = tmp.RTEnd(20.0);
	TEST_REAL_SIMILAR(it->getRT(),30.0)
	it = tmp.RTEnd(30.0);
	TEST_REAL_SIMILAR(it->getRT(),40.0)
	it = tmp.RTEnd(31.0);
	TEST_REAL_SIMILAR(it->getRT(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true)
END_SECTION

START_SECTION((ConstIterator RTBegin(CoordinateType rt) const))
	PeakMap tmp;
	MSSpectrum s;
	Peak1D p;

	s.setRT(30.0);
	tmp.addSpectrum(s);
	s.setRT(40.0);
	tmp.addSpectrum(s);
	s.setRT(45.0);
	tmp.addSpectrum(s);
	s.setRT(50.0);
	tmp.addSpectrum(s);

	PeakMap::Iterator it;

	it = tmp.RTBegin(20.0);
	TEST_REAL_SIMILAR(it->getRT(),30.0)
	it = tmp.RTBegin(30.0);
	TEST_REAL_SIMILAR(it->getRT(),30.0)
	it = tmp.RTBegin(31.0);
	TEST_REAL_SIMILAR(it->getRT(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true)
END_SECTION

START_SECTION((ConstIterator RTEnd(CoordinateType rt) const))
	PeakMap tmp;
	MSSpectrum s;
	Peak1D p;

	s.setRT(30.0);
	tmp.addSpectrum(s);
	s.setRT(40.0);
	tmp.addSpectrum(s);
	s.setRT(45.0);
	tmp.addSpectrum(s);
	s.setRT(50.0);
	tmp.addSpectrum(s);

	PeakMap::Iterator it;

	it = tmp.RTEnd(20.0);
	TEST_REAL_SIMILAR(it->getRT(),30.0)
	it = tmp.RTEnd(30.0);
	TEST_REAL_SIMILAR(it->getRT(),40.0)
	it = tmp.RTEnd(31.0);
	TEST_REAL_SIMILAR(it->getRT(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true)
END_SECTION

START_SECTION((void sortSpectra(bool sort_mz = true)))
	std::vector< Peak2D> plist;

	Peak2D p1;
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 5.0;
	plist.push_back(p1);

	Peak2D p2;
	p2.getPosition()[0] = 1.0;
	p2.getPosition()[1] = 3.0;
	plist.push_back(p2);

	Peak2D p3;
	p3.getPosition()[0] = 2.0;
	p3.getPosition()[1] = 14.0;
	plist.push_back(p3);

	Peak2D p4;
	p4.getPosition()[0] = 2.0;
	p4.getPosition()[1] = 11.0;
	plist.push_back(p4);

	PeakMap exp;
	exp.set2DData(plist);

	exp.sortSpectra(true);

	TEST_REAL_SIMILAR(exp[0][0].getMZ(),3.0);
	TEST_REAL_SIMILAR(exp[0][1].getMZ(),5.0);
	TEST_REAL_SIMILAR(exp[1][0].getMZ(),11.0);
	TEST_REAL_SIMILAR(exp[1][1].getMZ(),14.0);
END_SECTION

START_SECTION(bool isSorted(bool check_mz = true ) const)
	//make test dataset
	PeakMap exp;
	exp.resize(2);
	exp[0].setRT(1.0);
	exp[1].setRT(2.0);
	
	Peak1D p;
	p.setIntensity(1.0);
	p.setMZ(1000.0);
	exp[0].push_back(p);
	exp[1].push_back(p);
	
	p.setIntensity(1.0);
	p.setMZ(1001.0);
	exp[0].push_back(p);
	exp[1].push_back(p);
	
	p.setIntensity(1.0);
	p.setMZ(1002.0);
	exp[0].push_back(p);
	exp[1].push_back(p);
	
	//test with identical RTs
	TEST_EQUAL(exp.isSorted(false),true)
	TEST_EQUAL(exp.isSorted(),true)
	
	//test with acending RTs
	exp[0].setRT(1.0);
	exp[1].setRT(2.0);
	TEST_EQUAL(exp.isSorted(false),true)
	TEST_EQUAL(exp.isSorted(),true)

	//test with a reversed spectrum
	reverse(exp[0].begin(),exp[0].end());
	TEST_EQUAL(exp.isSorted(false),true)
	TEST_EQUAL(exp.isSorted(),false)

	//test with reversed RTs
	reverse(exp.begin(),exp.end());
	TEST_EQUAL(exp.isSorted(false),false)
	TEST_EQUAL(exp.isSorted(),false)
END_SECTION

START_SECTION((void reset()))
	std::vector< Peak2D> plist;

	Peak2D p;
	p.getPosition()[0] = 1.0;
	p.getPosition()[1] = 5.0;
	plist.push_back(p);
	p.getPosition()[0] = 2.0;
	p.getPosition()[1] = 3.0;
	plist.push_back(p);

	PeakMap exp;
	exp.set2DData(plist);
	exp.updateRanges();

	exp.reset();

	TEST_EQUAL(exp==PeakMap(),true);
END_SECTION

START_SECTION((const ExperimentalSettings& getExperimentalSettings() const))
	PeakMap exp;
	exp.setComment("test");
	TEST_EQUAL(exp.getExperimentalSettings().getComment(),"test");
END_SECTION

START_SECTION((ExperimentalSettings& getExperimentalSettings()))
	PeakMap exp;
	exp.getExperimentalSettings().setComment("test");
	TEST_EQUAL(exp.getExperimentalSettings().getComment(),"test");
END_SECTION

START_SECTION((MSExperiment& operator=(const ExperimentalSettings &source)))
	PeakMap exp,exp2;
	exp.getExperimentalSettings().setComment("test");
	exp2 = exp.getExperimentalSettings();
	TEST_EQUAL(exp2.getExperimentalSettings().getComment(),"test");
END_SECTION

START_SECTION((ConstIterator getPrecursorSpectrum(ConstIterator iterator) const))
	PeakMap exp;
	exp.resize(10);
	exp[0].setMSLevel(1);
	exp[1].setMSLevel(2);
	exp[2].setMSLevel(1);
	exp[3].setMSLevel(2);
	exp[4].setMSLevel(2);

	TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin())==exp.end(),true)
	TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+1)==exp.begin(),true)
	TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+2)==exp.end(),true)
	TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+3)==exp.begin()+2,true)
	TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+4)==exp.begin()+2,true)
	TEST_EQUAL(exp.getPrecursorSpectrum(exp.end())==exp.end(),true)

	exp[0].setMSLevel(2);
	exp[1].setMSLevel(1);
	exp[2].setMSLevel(1);
	exp[3].setMSLevel(1);
	exp[4].setMSLevel(1);

	TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin())==exp.end(),true)
	TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+1)==exp.end(),true)
	TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+2)==exp.end(),true)
	TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+3)==exp.end(),true)
	TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+4)==exp.end(),true)
	TEST_EQUAL(exp.getPrecursorSpectrum(exp.end())==exp.end(),true)

END_SECTION

START_SECTION((bool clearMetaDataArrays()))
	PeakMap exp;
	exp.resize(5);
	exp[0].getFloatDataArrays().resize(5);
	exp[0].getIntegerDataArrays().resize(5);
	exp[0].getStringDataArrays().resize(5);
	exp.clearMetaDataArrays();
	TEST_EQUAL(exp[0].getFloatDataArrays().size(),0)
	TEST_EQUAL(exp[0].getIntegerDataArrays().size(),0)
	TEST_EQUAL(exp[0].getStringDataArrays().size(),0)
END_SECTION

START_SECTION((void swap(MSExperiment &from)))
	PeakMap exp1, exp2;
	exp1.setComment("stupid comment");
	exp1.resize(1);
	exp1[0].setMSLevel(2);
	exp1[0].resize(2);
	exp1[0][0].setIntensity(0.5f);
	exp1[0][1].setIntensity(1.7f);
	exp1.updateRanges();
	
	exp1.swap(exp2);
	
	TEST_EQUAL(exp1.getComment(),"")
	TEST_EQUAL(exp1.size(),0)
	TEST_REAL_SIMILAR(exp1.getMinInt(),DRange<1>().minPosition()[0])
	TEST_EQUAL(exp1.getMSLevels().size(),0)
	TEST_EQUAL(exp1.getSize(),0);
	
	TEST_EQUAL(exp2.getComment(),"stupid comment")
	TEST_EQUAL(exp2.size(),1)
	TEST_REAL_SIMILAR(exp2.getMinInt(),0.5)
	TEST_EQUAL(exp2.getMSLevels().size(),1)
	TEST_EQUAL(exp2.getSize(),2);
	
END_SECTION

START_SECTION(void clear(bool clear_meta_data))
  PeakMap edit;
  edit.getSample().setName("bla");
	edit.resize(5);
	edit.updateRanges();
	edit.setMetaValue("label",String("bla"));
	vector<MSChromatogram > tmp;
	tmp.resize(5);
	edit.setChromatograms(tmp);

	edit.clear(false);
	TEST_EQUAL(edit.size(),0)
	TEST_EQUAL(edit==PeakMap(),false)

	edit.clear(true);
	TEST_EQUAL(edit==PeakMap(),true)
END_SECTION

START_SECTION((void sortChromatograms(bool sort_rt=true)))
  PeakMap exp;
  MSChromatogram chrom1, chrom2;
  ChromatogramPeak p1, p2, p3;
  p1.setRT(0.3);
  p1.setIntensity(10.0f);
  p2.setRT(0.2);
  p2.setIntensity(10.2f);
  p3.setRT(0.1);
  p3.setIntensity(10.4f);

	Product prod1;
	prod1.setMZ(100.0);
	chrom1.setProduct(prod1);
  chrom1.push_back(p1);
  chrom1.push_back(p2);

	Product prod2;
	prod2.setMZ(80.0);
	chrom2.setProduct(prod2);
  chrom2.push_back(p2);
  chrom2.push_back(p3);

  vector<MSChromatogram > chroms;
  chroms.push_back(chrom1);
  chroms.push_back(chrom2);
  exp.setChromatograms(chroms);
	TEST_EQUAL(exp.getChromatograms().size(), 2)
	TEST_REAL_SIMILAR(exp.getChromatograms()[0].getMZ(), 100.0)
	TEST_REAL_SIMILAR(exp.getChromatograms()[1].getMZ(), 80.0)

	// first sort without rt
	exp.sortChromatograms(false);
	TEST_REAL_SIMILAR(exp.getChromatograms()[0].getMZ(), 80.0)
	TEST_REAL_SIMILAR(exp.getChromatograms()[1].getMZ(), 100.0)

	TEST_REAL_SIMILAR(exp.getChromatograms()[1][0].getRT(), 0.3)
	TEST_REAL_SIMILAR(exp.getChromatograms()[1][1].getRT(), 0.2)

	// now also sort rt
	exp.sortChromatograms();

	TEST_REAL_SIMILAR(exp.getChromatograms()[0].getMZ(), 80.0)
	TEST_REAL_SIMILAR(exp.getChromatograms()[1].getMZ(), 100.0)

	TEST_REAL_SIMILAR(exp.getChromatograms()[1][0].getRT(), 0.2)
	TEST_REAL_SIMILAR(exp.getChromatograms()[1][1].getRT(), 0.3)

END_SECTION

START_SECTION((void setChromatograms(const std::vector< MSChromatogram > &chromatograms)))
	PeakMap exp;
	MSChromatogram chrom1, chrom2;
	ChromatogramPeak p1, p2, p3;
	p1.setRT(0.1);
	p1.setIntensity(10.0f);
	p2.setRT(0.2);
	p2.setIntensity(10.2f);
	p3.setRT(0.3);
	p3.setIntensity(10.4f);
	chrom1.push_back(p1);
	chrom1.push_back(p2);
	chrom2.push_back(p2);
	chrom2.push_back(p3);
	vector<MSChromatogram > chroms;
	chroms.push_back(chrom1);
	chroms.push_back(chrom2);
	exp.setChromatograms(chroms);
	TEST_EQUAL(exp.getChromatograms().size(), 2)
	TEST_EQUAL(exp.getChromatograms()[0] == chrom1, true)
	TEST_EQUAL(exp.getChromatograms()[1] == chrom2, true)
END_SECTION

START_SECTION((void addChromatogram(const MSChromatogram &chromatogram)))
  PeakMap exp;
  MSChromatogram chrom1, chrom2;
  ChromatogramPeak p1, p2, p3;
  p1.setRT(0.1);
  p1.setIntensity(10.0f);
  p2.setRT(0.2);
  p2.setIntensity(10.2f);
  p3.setRT(0.3);
  p3.setIntensity(10.4f);
  chrom1.push_back(p1);
  chrom1.push_back(p2);
  chrom2.push_back(p2);
  chrom2.push_back(p3);

	TEST_EQUAL(exp.getChromatograms().size(), 0)
	exp.addChromatogram(chrom1);
	TEST_EQUAL(exp.getChromatograms().size(), 1)
	TEST_EQUAL(exp.getChromatograms()[0] == chrom1, true)
	exp.addChromatogram(chrom2);
	TEST_EQUAL(exp.getChromatograms().size(), 2)
	TEST_EQUAL(exp.getChromatograms()[0] == chrom1, true)	
	TEST_EQUAL(exp.getChromatograms()[1] == chrom2, true)	
END_SECTION

START_SECTION((const std::vector<MSChromatogram >& getChromatograms() const))
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION((std::vector<MSChromatogram >& getChromatograms()))
  PeakMap exp;
  vector<MSChromatogram > chromatograms(2);
  exp.getChromatograms().swap(chromatograms);
  TEST_EQUAL(exp.getChromatograms().size(), 2);
  TEST_EQUAL(chromatograms.size(), 0);
  exp.getChromatograms().swap(chromatograms);
  TEST_EQUAL(exp.getChromatograms().size(), 0);
  TEST_EQUAL(chromatograms.size(), 2);
END_SECTION

START_SECTION((const MSChromatogram getTIC() const))
  PeakMap tmp;
  tmp.resize(2);
  Peak1D p;
  p.setMZ(5.0);
  p.setIntensity(3);
  tmp[0].push_back(p);
  p.setMZ(10.0);
  p.setIntensity(5);
  tmp[0].push_back(p);
  p.setMZ(5.0);
  p.setIntensity(2);
  tmp[1].push_back(p);
  tmp.updateRanges();
  MSChromatogram chrom = tmp.getTIC();
  TEST_EQUAL(chrom.size(), 2);
  TEST_EQUAL(chrom[0].getIntensity(), 8);
  TEST_EQUAL(chrom[1].getIntensity(), 2);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

