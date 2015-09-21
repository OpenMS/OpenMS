// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>

#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <string>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(FeatureMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


FeatureMap* pl_ptr = 0;
FeatureMap* nullPointer = 0;
START_SECTION((FeatureMap()))
	pl_ptr = new FeatureMap();
  TEST_NOT_EQUAL(pl_ptr, nullPointer)

	TEST_EQUAL(pl_ptr->getMin(), FeatureMap::PositionType::maxPositive())
	TEST_EQUAL(pl_ptr->getMax(), FeatureMap::PositionType::minNegative())
	TEST_REAL_SIMILAR(pl_ptr->getMinInt(), numeric_limits<double>::max())
	TEST_REAL_SIMILAR(pl_ptr->getMaxInt(), -numeric_limits<double>::max())
END_SECTION

START_SECTION((virtual ~FeatureMap()))
	delete pl_ptr;
END_SECTION

std::vector<PeptideIdentification> ids(1);
PeptideHit hit;
hit.setSequence(AASequence::fromString("ABCDE"));
ids[0].setHits(std::vector<PeptideHit>(1, hit));

Feature feature1;
feature1.getPosition()[0] = 2.0;
feature1.getPosition()[1] = 3.0;
feature1.setIntensity(1.0f);
feature1.setPeptideIdentifications(ids);  // single hit

Feature feature2;
feature2.getPosition()[0] = 0.0;
feature2.getPosition()[1] = 2.5;
feature2.setIntensity(0.5f);
ids.resize(2);
ids[1].setHits(std::vector<PeptideHit>(1, hit)); // same as first hit
feature2.setPeptideIdentifications(ids);

Feature feature3;
feature3.getPosition()[0] = 10.5;
feature3.getPosition()[1] = 0.0;
feature3.setIntensity(0.01f);
hit.setSequence(AASequence::fromString("KRGH"));
ids[1].setHits(std::vector<PeptideHit>(1, hit)); // different to first hit
feature3.setPeptideIdentifications(ids);

//feature with convex hulls
Feature feature4;
feature4.getPosition()[0] = 5.25;
feature4.getPosition()[1] = 1.5;
feature4.setIntensity(0.5f);
std::vector< ConvexHull2D > hulls(1);
hulls[0].addPoint(DPosition<2>(-1.0,2.0));
hulls[0].addPoint(DPosition<2>(4.0,1.2));
hulls[0].addPoint(DPosition<2>(5.0,3.123));
feature4.setConvexHulls(hulls);

START_SECTION((const std::vector<ProteinIdentification>& getProteinIdentifications() const))
	FeatureMap tmp;
	TEST_EQUAL(tmp.getProteinIdentifications().size(),0)
END_SECTION

START_SECTION((std::vector<ProteinIdentification>& getProteinIdentifications()))
	FeatureMap tmp;
	tmp.getProteinIdentifications().resize(1);
	TEST_EQUAL(tmp.getProteinIdentifications().size(),1)
END_SECTION

START_SECTION((void setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications)))
	FeatureMap tmp;
	tmp.setProteinIdentifications(std::vector<ProteinIdentification>(2));
	TEST_EQUAL(tmp.getProteinIdentifications().size(),2)
END_SECTION

START_SECTION((const std::vector<PeptideIdentification>& getUnassignedPeptideIdentifications() const))
	FeatureMap tmp;
	TEST_EQUAL(tmp.getUnassignedPeptideIdentifications().size(),0)
END_SECTION

START_SECTION((std::vector<PeptideIdentification>& getUnassignedPeptideIdentifications()))
	FeatureMap tmp;
	tmp.getUnassignedPeptideIdentifications().resize(1);
	TEST_EQUAL(tmp.getUnassignedPeptideIdentifications().size(),1)
END_SECTION

START_SECTION((void setUnassignedPeptideIdentifications(const std::vector<PeptideIdentification>& unassigned_peptide_identifications)))
	FeatureMap tmp;
	tmp.setUnassignedPeptideIdentifications(std::vector<PeptideIdentification>(2));
	TEST_EQUAL(tmp.getUnassignedPeptideIdentifications().size(),2)
END_SECTION

START_SECTION((const std::vector<DataProcessing>& getDataProcessing() const))
  FeatureMap tmp;
  TEST_EQUAL(tmp.getDataProcessing().size(),0);
END_SECTION

START_SECTION((std::vector<DataProcessing>& getDataProcessing()))
  FeatureMap tmp;
  tmp.getDataProcessing().resize(1);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION((void setDataProcessing(const std::vector< DataProcessing > &processing_method)))
  FeatureMap tmp;
  std::vector<DataProcessing> dummy;
  dummy.resize(1);
  tmp.setDataProcessing(dummy);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION((void updateRanges()))
	//test without convex hulls
  FeatureMap s;
  s.push_back(feature1);
  s.push_back(feature2);
  s.push_back(feature3);

  s.updateRanges();
  s.updateRanges(); //second time to check the initialization

  TEST_REAL_SIMILAR(s.getMaxInt(),1.0)
  TEST_REAL_SIMILAR(s.getMinInt(),0.01)
  TEST_REAL_SIMILAR(s.getMax()[0],10.5)
  TEST_REAL_SIMILAR(s.getMax()[1],3.0)
  TEST_REAL_SIMILAR(s.getMin()[0],0.0)
  TEST_REAL_SIMILAR(s.getMin()[1],0.0)

  //test with convex hull
  s.push_back(feature4);
  s.updateRanges();
  TEST_REAL_SIMILAR(s.getMaxInt(),1.0)
  TEST_REAL_SIMILAR(s.getMinInt(),0.01)
  TEST_REAL_SIMILAR(s.getMax()[0],10.5)
  TEST_REAL_SIMILAR(s.getMax()[1],3.123)
  TEST_REAL_SIMILAR(s.getMin()[0],-1.0)
  TEST_REAL_SIMILAR(s.getMin()[1],0.0)

END_SECTION

START_SECTION((FeatureMap(const FeatureMap &source)))
  FeatureMap map1;
  map1.setMetaValue("meta",String("value"));
  map1.push_back(feature1);
  map1.push_back(feature2);
  map1.push_back(feature3);
  map1.updateRanges();
  map1.setIdentifier("lsid");;
  map1.getDataProcessing().resize(1);
  map1.getProteinIdentifications().resize(1);
  map1.getUnassignedPeptideIdentifications().resize(1);

  FeatureMap map2(map1);

  TEST_EQUAL(map2.size(),3);
  TEST_EQUAL(map2.getMetaValue("meta").toString(),"value")
  TEST_REAL_SIMILAR(map2.getMaxInt(),1.0)
  TEST_STRING_EQUAL(map2.getIdentifier(),"lsid")
  TEST_EQUAL(map2.getDataProcessing().size(),1)
  TEST_EQUAL(map2.getProteinIdentifications().size(),1);
  TEST_EQUAL(map2.getUnassignedPeptideIdentifications().size(),1);
END_SECTION

START_SECTION((FeatureMap& operator = (const FeatureMap& rhs)))
	FeatureMap map1;
  map1.setMetaValue("meta",String("value"));
	map1.push_back(feature1);
	map1.push_back(feature2);
	map1.push_back(feature3);
	map1.updateRanges();
	map1.setIdentifier("lsid");
	map1.getDataProcessing().resize(1);
	map1.getProteinIdentifications().resize(1);
	map1.getUnassignedPeptideIdentifications().resize(1);

	//assignment
	FeatureMap map2;
	map2 = map1;

	TEST_EQUAL(map2.size(),3);
  TEST_EQUAL(map2.getMetaValue("meta").toString(),"value")
  TEST_REAL_SIMILAR(map2.getMaxInt(),1.0)
  TEST_STRING_EQUAL(map2.getIdentifier(),"lsid")
  TEST_EQUAL(map2.getDataProcessing().size(),1)
	TEST_EQUAL(map2.getProteinIdentifications().size(),1);
	TEST_EQUAL(map2.getUnassignedPeptideIdentifications().size(),1);

  //assignment of empty object
   map2 = FeatureMap();

	TEST_EQUAL(map2.size(),0);
	TEST_REAL_SIMILAR(map2.getMinInt(), numeric_limits<double>::max())
	TEST_REAL_SIMILAR(map2.getMaxInt(), -numeric_limits<double>::max())
  TEST_STRING_EQUAL(map2.getIdentifier(),"")
  TEST_EQUAL(map2.getDataProcessing().size(),0)
	TEST_EQUAL(map2.getProteinIdentifications().size(),0);
	TEST_EQUAL(map2.getUnassignedPeptideIdentifications().size(),0);
END_SECTION

START_SECTION((bool operator == (const FeatureMap& rhs) const))
	FeatureMap empty,edit;

	TEST_EQUAL(empty==edit, true);

	edit.setIdentifier("lsid");;
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.push_back(feature1);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.getDataProcessing().resize(1);
	TEST_EQUAL(empty==edit, false);

  edit = empty;
	edit.getProteinIdentifications().resize(1);
  TEST_EQUAL(edit==empty, false);

	edit = empty;
	edit.getUnassignedPeptideIdentifications().resize(10);
	TEST_EQUAL(empty==edit, false);


	edit = empty;
	edit.push_back(feature1);
	edit.push_back(feature2);
	edit.updateRanges();
	edit.clear(false);
	TEST_EQUAL(empty==edit, false);
END_SECTION

START_SECTION((bool operator != (const FeatureMap& rhs) const))
	FeatureMap empty,edit;

	TEST_EQUAL(empty!=edit, false);

	edit.setIdentifier("lsid");;
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.push_back(feature1);
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.getDataProcessing().resize(1);
	TEST_EQUAL(empty!=edit, true);

  edit = empty;
	edit.getProteinIdentifications().resize(10);
  TEST_EQUAL(edit!=empty, true);

	edit = empty;
	edit.getUnassignedPeptideIdentifications().resize(10);
	TEST_EQUAL(empty!=edit, true);


	edit = empty;
	edit.push_back(feature1);
	edit.push_back(feature2);
	edit.updateRanges();
	edit.clear(false);
	TEST_EQUAL(empty!=edit, true);
END_SECTION

START_SECTION((FeatureMap operator + (const FeatureMap& rhs) const))
	// just some basic testing... most is done in operator +=()
	FeatureMap m1, m2, m3;

	TEST_EQUAL(m1+m2, m3);

	Feature f1;
	f1.setMZ(100.12);
	m1.push_back(f1);
	m3 = m1;
	TEST_EQUAL(m1+m2, m3);

END_SECTION

START_SECTION((FeatureMap& operator+= (const FeatureMap& rhs)))
	FeatureMap m1, m2, m3;

	// adding empty maps has no effect:
	m1+=m2;
	TEST_EQUAL(m1, m3);

	// with content:
	Feature f1;
	f1.setMZ(100.12);
	m1.push_back(f1);
	m3 = m1;
	m1+=m2;
	TEST_EQUAL(m1, m3);

	// test basic classes
	m1.setIdentifier ("123");
	m1.getDataProcessing().resize(1);
	m1.getProteinIdentifications().resize(1);
	m1.getUnassignedPeptideIdentifications().resize(1);
	m1.ensureUniqueId();

	m2.setIdentifier ("321");
	m2.getDataProcessing().resize(2);
	m2.getProteinIdentifications().resize(2);
	m2.getUnassignedPeptideIdentifications().resize(2);
	m2.push_back(Feature());
	m2.push_back(Feature());


	m1+=m2;
	TEST_EQUAL(m1.getIdentifier(), "");
	TEST_EQUAL(UniqueIdInterface::isValid(m1.getUniqueId()), false);
	TEST_EQUAL(m1.getDataProcessing().size(), 3);
	TEST_EQUAL(m1.getProteinIdentifications().size(),3);
	TEST_EQUAL(m1.getUnassignedPeptideIdentifications().size(),3);
	TEST_EQUAL(m1.size(),3);


END_SECTION

START_SECTION((void sortByIntensity(bool reverse=false)))

	FeatureMap to_be_sorted;

	Feature f1;
	f1.setIntensity(10.0f);
	to_be_sorted.push_back(f1);

	Feature f2;
	f2.setIntensity(5.0f);
	to_be_sorted.push_back(f2);

	Feature f3;
	f3.setIntensity(3.0f);
	to_be_sorted.push_back(f3);

	to_be_sorted.sortByIntensity();

	TEST_EQUAL(to_be_sorted[0].getIntensity(),3);
	TEST_EQUAL(to_be_sorted[1].getIntensity(),5);
	TEST_EQUAL(to_be_sorted[2].getIntensity(),10);


	to_be_sorted.sortByIntensity(true);

	TEST_EQUAL(to_be_sorted[0].getIntensity(),10);
	TEST_EQUAL(to_be_sorted[1].getIntensity(),5);
	TEST_EQUAL(to_be_sorted[2].getIntensity(),3);

END_SECTION

START_SECTION((void sortByPosition()))

	FeatureMap to_be_sorted;

	Feature f1;
	f1.getPosition()[0] = 10;
	to_be_sorted.push_back(f1);

	Feature f2;
	f2.getPosition()[0] = 5;
	to_be_sorted.push_back(f2);

	Feature f3;
	f3.getPosition()[0] = 3;
	to_be_sorted.push_back(f3);

	to_be_sorted.sortByPosition();

	TEST_EQUAL(to_be_sorted[0].getPosition()[0],3);
	TEST_EQUAL(to_be_sorted[1].getPosition()[0],5);
	TEST_EQUAL(to_be_sorted[2].getPosition()[0],10);

END_SECTION

START_SECTION((void sortByMZ()))

	FeatureMap to_be_sorted;

	Feature f1;
	f1.getPosition()[0] = 10;
	f1.getPosition()[1] = 25;
	to_be_sorted.push_back(f1);

	Feature f2;
	f2.getPosition()[0] = 5;
	f2.getPosition()[1] = 15;
	to_be_sorted.push_back(f2);

	Feature f3;
	f3.getPosition()[0] = 3;
	f3.getPosition()[1] = 10;
	to_be_sorted.push_back(f3);

	to_be_sorted.sortByMZ();

	TEST_EQUAL(to_be_sorted[0].getPosition()[1],10);
	TEST_EQUAL(to_be_sorted[1].getPosition()[1],15);
	TEST_EQUAL(to_be_sorted[2].getPosition()[1],25);

END_SECTION

START_SECTION((void sortByRT()))

	FeatureMap to_be_sorted;

	Feature f1;
	f1.getPosition()[0] = 10;
	f1.getPosition()[1] = 25;
	to_be_sorted.push_back(f1);

	Feature f2;
	f2.getPosition()[0] = 5;
	f2.getPosition()[1] = 15;
	to_be_sorted.push_back(f2);

	Feature f3;
	f3.getPosition()[0] = 3;
	f3.getPosition()[1] = 10;
	to_be_sorted.push_back(f3);

	to_be_sorted.sortByRT();

	TEST_EQUAL(to_be_sorted[0].getPosition()[0],3);
	TEST_EQUAL(to_be_sorted[1].getPosition()[0],5);
	TEST_EQUAL(to_be_sorted[2].getPosition()[0],10);

END_SECTION

START_SECTION((void swap(FeatureMap& from)))
{
  FeatureMap map1, map2;
	map1.setIdentifier("stupid comment");
	map1.push_back(feature1);
	map1.push_back(feature2);
	map1.updateRanges();
	map1.getDataProcessing().resize(1);
	map1.getProteinIdentifications().resize(1);
	map1.getUnassignedPeptideIdentifications().resize(1);

	map1.swap(map2);

	TEST_EQUAL(map1.getIdentifier(),"")
	TEST_EQUAL(map1.size(),0)
	TEST_REAL_SIMILAR(map1.getMinInt(),DRange<1>().minPosition()[0])
  TEST_EQUAL(map1.getDataProcessing().size(),0)
	TEST_EQUAL(map1.getProteinIdentifications().size(),0);
	TEST_EQUAL(map1.getUnassignedPeptideIdentifications().size(),0);

	TEST_EQUAL(map2.getIdentifier(),"stupid comment")
	TEST_EQUAL(map2.size(),2)
	TEST_REAL_SIMILAR(map2.getMinInt(),0.5)
  TEST_EQUAL(map2.getDataProcessing().size(),1)
	TEST_EQUAL(map2.getProteinIdentifications().size(),1);
	TEST_EQUAL(map2.getUnassignedPeptideIdentifications().size(),1);
}
END_SECTION

START_SECTION((void swapFeaturesOnly(FeatureMap& from)))
{
  FeatureMap map1, map2;
	map1.setIdentifier("stupid comment");
	map1.push_back(feature1);
	map1.push_back(feature2);
	map1.updateRanges();
	map1.getDataProcessing().resize(1);
	map1.getProteinIdentifications().resize(1);
	map1.getUnassignedPeptideIdentifications().resize(1);

	map1.swapFeaturesOnly(map2);

	TEST_EQUAL(map1.getIdentifier(),"stupid comment")
	TEST_EQUAL(map1.size(),0)
	TEST_REAL_SIMILAR(map1.getMinInt(),DRange<1>().minPosition()[0])
  TEST_EQUAL(map1.getDataProcessing().size(),1)
	TEST_EQUAL(map1.getProteinIdentifications().size(),1);
	TEST_EQUAL(map1.getUnassignedPeptideIdentifications().size(),1);

	TEST_EQUAL(map2.getIdentifier(),"")
	TEST_EQUAL(map2.size(),2)
	TEST_REAL_SIMILAR(map2.getMinInt(),0.5)
  TEST_EQUAL(map2.getDataProcessing().size(),0)
	TEST_EQUAL(map2.getProteinIdentifications().size(),0);
	TEST_EQUAL(map2.getUnassignedPeptideIdentifications().size(),0);
}
END_SECTION

START_SECTION((void sortByOverallQuality(bool reverse=false)))

	FeatureMap to_be_sorted;

	Feature f1;
	f1.getPosition()[0] = 1;
	f1.getPosition()[1] = 1;
	f1.setOverallQuality(10);
	to_be_sorted.push_back(f1);

	Feature f2;
	f2.getPosition()[0] = 2;
	f2.getPosition()[1] = 2;
	f2.setOverallQuality(30);
	to_be_sorted.push_back(f2);

	Feature f3;
	f3.getPosition()[0] = 3;
	f3.getPosition()[1] = 3;
	f3.setOverallQuality(20);
	to_be_sorted.push_back(f3);

	to_be_sorted.sortByOverallQuality();

	TEST_EQUAL(to_be_sorted[0].getPosition()[0],1);
	TEST_EQUAL(to_be_sorted[1].getPosition()[0],3);
	TEST_EQUAL(to_be_sorted[2].getPosition()[0],2);

	TEST_EQUAL(to_be_sorted[0].getOverallQuality(),10);
	TEST_EQUAL(to_be_sorted[1].getOverallQuality(),20);
	TEST_EQUAL(to_be_sorted[2].getOverallQuality(),30);

	to_be_sorted.sortByOverallQuality(true);

	TEST_EQUAL(to_be_sorted[0].getPosition()[0],2);
	TEST_EQUAL(to_be_sorted[1].getPosition()[0],3);
	TEST_EQUAL(to_be_sorted[2].getPosition()[0],1);

	TEST_EQUAL(to_be_sorted[0].getOverallQuality(),30);
	TEST_EQUAL(to_be_sorted[1].getOverallQuality(),20);
	TEST_EQUAL(to_be_sorted[2].getOverallQuality(),10);

END_SECTION

START_SECTION((void clear(bool clear_meta_data=true)))
  FeatureMap map1;
	map1.setIdentifier("stupid comment");
	map1.push_back(feature1);
	map1.push_back(feature2);
	map1.updateRanges();
	map1.getDataProcessing().resize(1);
	map1.getProteinIdentifications().resize(1);
	map1.getUnassignedPeptideIdentifications().resize(1);

	map1.clear(false);
	TEST_EQUAL(map1.size(),0)
	TEST_EQUAL(map1==FeatureMap(),false)

	map1.clear(true);
	TEST_EQUAL(map1==FeatureMap(),true)
END_SECTION


START_SECTION(([EXTRA] void uniqueIdToIndex()))
{
	  FeatureMap fm;
	  Feature f;
	  f.setMZ(23.9);
	  std::vector< std::pair < Size, UInt64 > > pairs;
	  const Size num_features = 4;
	  for ( Size i = 0; i < num_features; ++i )
	  {
	    f.setRT(i*100);
	    f.setUniqueId();
	    pairs.push_back(make_pair(i,f.getUniqueId()));
      fm.push_back(f);
	  }
    for ( Size i = 0; i < num_features; ++i )
    {
      TEST_EQUAL(fm.uniqueIdToIndex(pairs[i].second),pairs[i].first);
    }
    STATUS("shuffling ...");
    std::random_shuffle(pairs.begin(),pairs.end());
    std::random_shuffle(fm.begin(),fm.end());
    for ( Size i = 0; i < num_features; ++i )
    {
      STATUS("pairs[i]:  " << pairs[i].first << ", " << pairs[i].second )
      TEST_EQUAL(fm.uniqueIdToIndex(fm[pairs[i].first].getUniqueId()),pairs[i].first);
      TEST_EQUAL(fm[fm.uniqueIdToIndex(pairs[i].second)].getUniqueId(),pairs[i].second);
    }

    f.setRT(98765421);
    f.setUniqueId();
    pairs.push_back(make_pair(987654321,f.getUniqueId()));

    TEST_EQUAL(fm.uniqueIdToIndex(pairs.back().second),Size(-1));
    fm.push_back(f);
    TEST_EQUAL(fm.uniqueIdToIndex(pairs.back().second),fm.size()-1);

    fm.push_back(Feature());
    fm.push_back(f);
    fm.push_back(Feature());
    fm.push_back(Feature());
    STATUS("fm: " << fm);
    fm.erase(fm.begin()+1);
    fm.erase(fm.begin()+2);
    STATUS("fm: " << fm);
    TEST_EXCEPTION_WITH_MESSAGE(Exception::Postcondition,fm.updateUniqueIdToIndex(),"Duplicate valid unique ids detected!   RandomAccessContainer has size()==7, num_valid_unique_id==4, uniqueid_to_index_.size()==3");
}
END_SECTION

START_SECTION((template < typename Type > Size applyMemberFunction(Size(Type::*member_function)())))
{
  FeatureMap fm;
  fm.push_back(Feature());
  fm.push_back(Feature());
  fm.back().getSubordinates().push_back(Feature());

  TEST_EQUAL(fm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),4);
  fm.setUniqueId();
  TEST_EQUAL(fm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),3);
  fm.applyMemberFunction(&UniqueIdInterface::setUniqueId);
  TEST_EQUAL(fm.applyMemberFunction(&UniqueIdInterface::hasValidUniqueId),4);
  TEST_EQUAL(fm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),0);
  fm.begin()->clearUniqueId();
  TEST_EQUAL(fm.applyMemberFunction(&UniqueIdInterface::hasValidUniqueId),3);
  TEST_EQUAL(fm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),1);
}
END_SECTION

START_SECTION((template < typename Type > Size applyMemberFunction(Size(Type::*member_function)() const ) const ))
{
  FeatureMap fm;
  FeatureMap const & fmc(fm);
  fm.push_back(Feature());
  fm.push_back(Feature());
  fm.back().getSubordinates().push_back(Feature());

  TEST_EQUAL(fmc.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),4);
  fm.setUniqueId();
  TEST_EQUAL(fmc.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),3);
  fm.applyMemberFunction(&UniqueIdInterface::setUniqueId);
  TEST_EQUAL(fmc.applyMemberFunction(&UniqueIdInterface::hasValidUniqueId),4);
  TEST_EQUAL(fm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),0);
  fm.begin()->clearUniqueId();
  TEST_EQUAL(fmc.applyMemberFunction(&UniqueIdInterface::hasValidUniqueId),3);
  TEST_EQUAL(fmc.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),1);
}
END_SECTION

START_SECTION((  AnnotationStatistics getAnnotationStatistics() const ))
  FeatureMap fm;

  AnnotationStatistics stats, res;
  stats = fm.getAnnotationStatistics();
  TEST_EQUAL(stats == res, true)

  fm.push_back(feature1); // single hit
  stats = fm.getAnnotationStatistics();
  ++res.states[BaseFeature::FEATURE_ID_SINGLE];
  std::cout << res;
  TEST_EQUAL(stats == res, true)

  fm.push_back(feature4); // single hit + no hit
  stats = fm.getAnnotationStatistics();
  ++res.states[BaseFeature::FEATURE_ID_NONE];
  std::cout << res;
  TEST_EQUAL(stats == res, true)

  fm.push_back(feature4); // single hit + 2x no hit
  stats = fm.getAnnotationStatistics();
  ++res.states[BaseFeature::FEATURE_ID_NONE];
  std::cout << res;
  TEST_EQUAL(stats == res, true)

  fm.push_back(feature2); // single hit + 2x no hit + multi-hit (same)
  stats = fm.getAnnotationStatistics();
  ++res.states[BaseFeature::FEATURE_ID_MULTIPLE_SAME];
  std::cout << res;
  TEST_EQUAL(stats == res, true)

  fm.push_back(feature3); // single hit + 2x no hit + multi-hit (same) + multi (divergent)
  stats = fm.getAnnotationStatistics();
  ++res.states[BaseFeature::FEATURE_ID_MULTIPLE_DIVERGENT];
  std::cout << res;
  std::cout << stats;
  TEST_EQUAL(stats == res, true)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
