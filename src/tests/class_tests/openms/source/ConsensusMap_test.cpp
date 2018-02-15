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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
///////////////////////////

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/DataProcessing.h>

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusMap* ptr = nullptr;
ConsensusMap* nullPointer = nullptr;
START_SECTION((ConsensusMap()))
	ptr = new ConsensusMap();
	TEST_NOT_EQUAL(ptr, nullPointer)
	TEST_EQUAL(ptr->isMetaEmpty(),true)
	TEST_REAL_SIMILAR(ptr->getMinInt(), numeric_limits<double>::max())
	TEST_REAL_SIMILAR(ptr->getMaxInt(), -numeric_limits<double>::max())
END_SECTION

START_SECTION((~ConsensusMap()))
	delete ptr;
END_SECTION

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
  ConsensusMap tmp;
  TEST_EQUAL(tmp.getDataProcessing().size(),0);
END_SECTION

START_SECTION((std::vector<DataProcessing>& getDataProcessing()))
  ConsensusMap tmp;
  tmp.getDataProcessing().resize(1);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION((void setDataProcessing(const std::vector< DataProcessing > &processing_method)))
  ConsensusMap tmp;
  std::vector<DataProcessing> dummy;
  dummy.resize(1);
  tmp.setDataProcessing(dummy);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

Feature feature1;
feature1.getPosition()[0] = 2.0;
feature1.getPosition()[1] = 3.0;
feature1.setIntensity(1.0f);

Feature feature2;
feature2.getPosition()[0] = 0.0;
feature2.getPosition()[1] = 2.5;
feature2.setIntensity(0.5f);

Feature feature3;
feature3.getPosition()[0] = 10.5;
feature3.getPosition()[1] = 0.0;
feature3.setIntensity(0.01f);

Feature feature4;
feature4.getPosition()[0] = 5.25;
feature4.getPosition()[1] = 1.5;
feature4.setIntensity(0.5f);

START_SECTION((void updateRanges()))
  ConsensusMap map;
  feature1.setUniqueId(1);
	ConsensusFeature f;
	f.setIntensity(1.0f);
	f.setRT(2.0);
	f.setMZ(3.0);
	f.insert(1,feature1);
	map.push_back(f);

  map.updateRanges();
  TEST_REAL_SIMILAR(map.getMaxInt(),1.0)
  TEST_REAL_SIMILAR(map.getMinInt(),1.0)
  TEST_REAL_SIMILAR(map.getMax()[0],2.0)
  TEST_REAL_SIMILAR(map.getMax()[1],3.0)
  TEST_REAL_SIMILAR(map.getMin()[0],2.0)
  TEST_REAL_SIMILAR(map.getMin()[1],3.0)

  //second time to check the initialization
  map.updateRanges();

  TEST_REAL_SIMILAR(map.getMaxInt(),1.0)
  TEST_REAL_SIMILAR(map.getMinInt(),1.0)
  TEST_REAL_SIMILAR(map.getMax()[0],2.0)
  TEST_REAL_SIMILAR(map.getMax()[1],3.0)
  TEST_REAL_SIMILAR(map.getMin()[0],2.0)
  TEST_REAL_SIMILAR(map.getMin()[1],3.0)

  //two points
  feature2.setUniqueId(2);
	f.insert(1,feature2);
	map.push_back(f);
	map.updateRanges();

  TEST_REAL_SIMILAR(map.getMaxInt(),1.0)
  TEST_REAL_SIMILAR(map.getMinInt(),0.5)
  TEST_REAL_SIMILAR(map.getMax()[0],2.0)
  TEST_REAL_SIMILAR(map.getMax()[1],3.0)
  TEST_REAL_SIMILAR(map.getMin()[0],0.0)
  TEST_REAL_SIMILAR(map.getMin()[1],2.5)

	//four points
  feature3.setUniqueId(3);
	f.insert(1,feature3);
  feature4.setUniqueId(4);
	f.insert(1,feature4);
	map.push_back(f);
	map.updateRanges();

  TEST_REAL_SIMILAR(map.getMaxInt(),1.0)
  TEST_REAL_SIMILAR(map.getMinInt(),0.01)
  TEST_REAL_SIMILAR(map.getMax()[0],10.5)
  TEST_REAL_SIMILAR(map.getMax()[1],3.0)
  TEST_REAL_SIMILAR(map.getMin()[0],0.0)
  TEST_REAL_SIMILAR(map.getMin()[1],0.0)

END_SECTION

START_SECTION((ConsensusMap& operator+=(const ConsensusMap &rhs)))
{
  ConsensusMap m1, m2, m3;
  // adding empty maps has no effect:
  m1+=m2;
  TEST_EQUAL(m1, m3);

  // with content:
  ConsensusFeature f1;
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
  m1.getFileDescriptions()[0].filename = "m1";

  m2.setIdentifier ("321");
  m2.getDataProcessing().resize(2);
  m2.getProteinIdentifications().resize(2);
  m2.getUnassignedPeptideIdentifications().resize(2);
  m2.push_back(ConsensusFeature());
  m2.push_back(ConsensusFeature());
  m2.getFileDescriptions()[1].filename = "m2";

  m1+=m2;
  TEST_EQUAL(m1.getIdentifier(), "");
  TEST_EQUAL(UniqueIdInterface::isValid(m1.getUniqueId()), false);
  TEST_EQUAL(m1.getDataProcessing().size(), 3);
  TEST_EQUAL(m1.getProteinIdentifications().size(),3);
  TEST_EQUAL(m1.getUnassignedPeptideIdentifications().size(),3);
  TEST_EQUAL(m1.size(),3);
  TEST_EQUAL(m1.getFileDescriptions().size(), 2);
}
END_SECTION

START_SECTION((ConsensusMap& operator = (const ConsensusMap& source)))
  ConsensusMap map1;
  map1.setMetaValue("meta",String("value"));
  map1.setIdentifier("lsid");
  map1.getFileDescriptions()[0].filename = "blub";
  map1.getFileDescriptions()[0].size = 47;
  map1.getFileDescriptions()[0].label = "label";
	map1.getFileDescriptions()[0].setMetaValue("meta",String("meta"));
	map1.getDataProcessing().resize(1);
	map1.setExperimentType("itraq");
	map1.getProteinIdentifications().resize(1);
	map1.getUnassignedPeptideIdentifications().resize(1);

  //assignment
  ConsensusMap map2;
  map2 = map1;
  TEST_EQUAL(map2.getIdentifier(),"lsid")
  TEST_EQUAL(map2.getMetaValue("meta").toString(),"value")
  TEST_EQUAL(map2.getFileDescriptions()[0].filename == "blub", true)
  TEST_EQUAL(map2.getFileDescriptions()[0].label == "label", true)
  TEST_EQUAL(map2.getFileDescriptions()[0].size == 47, true)
	TEST_EQUAL(map2.getFileDescriptions()[0].getMetaValue("meta") == "meta", true)
  TEST_EQUAL(map2.getExperimentType(), "itraq")
  TEST_EQUAL(map2.getDataProcessing().size(),1)
	TEST_EQUAL(map2.getProteinIdentifications().size(),1);
	TEST_EQUAL(map2.getUnassignedPeptideIdentifications().size(),1);

  //assignment of empty object
  map2 = ConsensusMap();
  TEST_EQUAL(map2.getIdentifier(),"")
  TEST_EQUAL(map2.getFileDescriptions().size(),0)
  TEST_EQUAL(map2.getExperimentType(),"")
  TEST_EQUAL(map2.getDataProcessing().size(),0)
	TEST_EQUAL(map2.getProteinIdentifications().size(),0);
	TEST_EQUAL(map2.getUnassignedPeptideIdentifications().size(),0);
END_SECTION

START_SECTION((ConsensusMap(const ConsensusMap& source)))
  ConsensusMap map1;
  map1.setMetaValue("meta",String("value"));
  map1.setIdentifier("lsid");
  map1.getFileDescriptions()[0].filename = "blub";
  map1.getFileDescriptions()[0].size = 47;
  map1.getFileDescriptions()[0].label = "label";
	map1.getFileDescriptions()[0].setMetaValue("meta",String("meta"));
	map1.getDataProcessing().resize(1);
	map1.setExperimentType("itraq");
	map1.getProteinIdentifications().resize(1);
	map1.getUnassignedPeptideIdentifications().resize(1);

  ConsensusMap map2(map1);

  TEST_EQUAL(map2.getIdentifier(),"lsid")
  TEST_EQUAL(map2.getMetaValue("meta").toString(),"value")
  TEST_EQUAL(map2.getFileDescriptions()[0].filename == "blub", true)
  TEST_EQUAL(map2.getFileDescriptions()[0].label == "label", true)
  TEST_EQUAL(map2.getFileDescriptions()[0].size == 47, true)
	TEST_EQUAL(map2.getFileDescriptions()[0].getMetaValue("meta") == "meta", true)
  TEST_EQUAL(map2.getExperimentType(),"itraq")
  TEST_EQUAL(map2.getDataProcessing().size(),1)
	TEST_EQUAL(map2.getProteinIdentifications().size(),1);
	TEST_EQUAL(map2.getUnassignedPeptideIdentifications().size(),1);
END_SECTION

START_SECTION((ConsensusMap(Base::size_type n)))
  ConsensusMap cons_map(5);

  TEST_EQUAL(cons_map.size(),5)
END_SECTION

/////

ConsensusMap::FileDescription* fd_ptr = nullptr;
ConsensusMap::FileDescription* fd_nullPointer = nullptr;

START_SECTION(([ConsensusMap::FileDescription] FileDescription()))
fd_ptr = new ConsensusMap::FileDescription();
TEST_NOT_EQUAL(fd_ptr, fd_nullPointer)
END_SECTION

START_SECTION((const FileDescriptions& getFileDescriptions() const))
  ConsensusMap cons_map;

  TEST_EQUAL(cons_map.getFileDescriptions().size(),0)
END_SECTION

START_SECTION((FileDescriptions& getFileDescriptions()))
  ConsensusMap cons_map;

  cons_map.getFileDescriptions()[0].filename = "blub";
  TEST_EQUAL(cons_map.getFileDescriptions()[0].filename == "blub", true)
END_SECTION

START_SECTION((const String& getExperimentType() const))
  ConsensusMap cons_map;
	TEST_EQUAL(cons_map.getExperimentType() == "", true)
END_SECTION

START_SECTION((void setExperimentType(const String& experiment_type)))
  ConsensusMap cons_map;
	cons_map.setExperimentType("itraq");
  TEST_EQUAL(cons_map.getExperimentType(),"itraq")
END_SECTION

START_SECTION((void swap(ConsensusMap& from)))
	ConsensusMap map1, map2;
	ConsensusFeature f;
	f.insert(1,Feature());
	map1.push_back(f);
  map1.getFileDescriptions()[1].filename = "bla";
	map1.getFileDescriptions()[1].size = 5;
	map1.setIdentifier("LSID");
	map1.setExperimentType("itraq");
	map1.getDataProcessing().resize(1);
	map1.getProteinIdentifications().resize(1);
	map1.getUnassignedPeptideIdentifications().resize(1);

	map1.swap(map2);

	TEST_EQUAL(map1.size(),0)
	TEST_EQUAL(map1.getFileDescriptions().size(),0)
	TEST_EQUAL(map1.getIdentifier(),"")
  TEST_EQUAL(map1.getDataProcessing().size(),0)
	TEST_EQUAL(map1.getProteinIdentifications().size(),0);
	TEST_EQUAL(map1.getUnassignedPeptideIdentifications().size(),0);

	TEST_EQUAL(map2.size(),1)
	TEST_EQUAL(map2.getFileDescriptions().size(),1)
	TEST_EQUAL(map2.getIdentifier(),"LSID")
  TEST_EQUAL(map2.getExperimentType(),"itraq")
  TEST_EQUAL(map2.getDataProcessing().size(),1)
	TEST_EQUAL(map2.getProteinIdentifications().size(),1);
	TEST_EQUAL(map2.getUnassignedPeptideIdentifications().size(),1);
END_SECTION

START_SECTION((bool operator == (const ConsensusMap& rhs) const))
	ConsensusMap empty,edit;

	TEST_EQUAL(empty==edit, true);

	edit.setIdentifier("lsid");;
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.push_back(ConsensusFeature(feature1));
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.getDataProcessing().resize(1);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.setMetaValue("bla", 4.1);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.getFileDescriptions()[0].filename = "bla";
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.setExperimentType("bla");
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.getProteinIdentifications().resize(10);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.getUnassignedPeptideIdentifications().resize(10);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.setExperimentType("bla");
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.push_back(ConsensusFeature(feature1));
	edit.push_back(ConsensusFeature(feature2));
	edit.updateRanges();
	edit.clear(false);
	TEST_EQUAL(empty==edit, false);
END_SECTION

START_SECTION((bool operator != (const ConsensusMap& rhs) const))
	ConsensusMap empty,edit;

	TEST_EQUAL(empty!=edit, false);

	edit.setIdentifier("lsid");;
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.push_back(ConsensusFeature(feature1));
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.getDataProcessing().resize(1);
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.setMetaValue("bla", 4.1);
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.getFileDescriptions()[0].filename = "bla";
	TEST_EQUAL(empty!=edit, true)

	edit = empty;
	edit.setExperimentType("bla");
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.getProteinIdentifications().resize(10);
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.getUnassignedPeptideIdentifications().resize(10);
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.push_back(ConsensusFeature(feature1));
	edit.push_back(ConsensusFeature(feature2));
	edit.updateRanges();
	edit.clear(false);
	TEST_EQUAL(empty!=edit, true);
END_SECTION


START_SECTION((void sortByIntensity(bool reverse=false)))
{
  NOT_TESTABLE; // tested within TOPP TextExporter
}
END_SECTION

START_SECTION((void sortByRT()))
{
  NOT_TESTABLE; // tested within TOPP TextExporter
}
END_SECTION

START_SECTION((void sortByMZ()))
{
  NOT_TESTABLE; // tested within TOPP TextExporter
}
END_SECTION

START_SECTION((void sortByPosition()))
{
  NOT_TESTABLE; // tested within TOPP TextExporter
}
END_SECTION

START_SECTION((void sortByQuality(bool reverse=false)))
{
  NOT_TESTABLE; // tested within TOPP TextExporter
}
END_SECTION

START_SECTION((void sortBySize()))
{
  NOT_TESTABLE; // tested within TOPP TextExporter
}
END_SECTION

START_SECTION((void sortByMaps()))
{
  NOT_TESTABLE; // tested within TOPP TextExporter
}
END_SECTION

START_SECTION((void sortPeptideIdentificationsByMapIndex()))
{
  NOT_TESTABLE; // tested within TOPP IDMapper
}
END_SECTION

START_SECTION((void clear(bool clear_meta_data = true)))
{
  ConsensusMap map1;
	ConsensusFeature f;
	f.insert(1,Feature());
	map1.push_back(f);
  map1.getFileDescriptions()[1].filename = "bla";
	map1.getFileDescriptions()[1].size = 5;
	map1.setIdentifier("LSID");
	map1.setExperimentType("itraq");
	map1.getDataProcessing().resize(1);
	map1.getProteinIdentifications().resize(1);
	map1.getUnassignedPeptideIdentifications().resize(1);
	
	map1.clear(false);
	TEST_EQUAL(map1.size(),0)
	TEST_EQUAL(map1==ConsensusMap(),false)

	map1.clear(true);
	TEST_EQUAL(map1==ConsensusMap(),true)
}
END_SECTION

START_SECTION((template < typename Type > Size applyMemberFunction(Size(Type::*member_function)())))
{
  ConsensusMap cm;
  cm.push_back(ConsensusFeature());
  cm.push_back(ConsensusFeature());
  cm.push_back(ConsensusFeature());

  TEST_EQUAL(cm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),4);
  cm.setUniqueId();
  TEST_EQUAL(cm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),3);
  cm.applyMemberFunction(&UniqueIdInterface::setUniqueId);
  TEST_EQUAL(cm.applyMemberFunction(&UniqueIdInterface::hasValidUniqueId),4);
  TEST_EQUAL(cm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),0);
  cm.begin()->clearUniqueId();
  TEST_EQUAL(cm.applyMemberFunction(&UniqueIdInterface::hasValidUniqueId),3);
  TEST_EQUAL(cm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),1);
}
END_SECTION

START_SECTION((template < typename Type > Size applyMemberFunction(Size(Type::*member_function)() const ) const ))
{
  ConsensusMap cm;
  ConsensusMap const & cmc(cm);
  cm.push_back(ConsensusFeature());
  cm.push_back(ConsensusFeature());
  cm.push_back(ConsensusFeature());

  TEST_EQUAL(cmc.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),4);
  cm.setUniqueId();
  TEST_EQUAL(cmc.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),3);
  cm.applyMemberFunction(&UniqueIdInterface::setUniqueId);
  TEST_EQUAL(cmc.applyMemberFunction(&UniqueIdInterface::hasValidUniqueId),4);
  TEST_EQUAL(cm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),0);
  cm.begin()->clearUniqueId();
  TEST_EQUAL(cmc.applyMemberFunction(&UniqueIdInterface::hasValidUniqueId),3);
  TEST_EQUAL(cmc.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId),1);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



