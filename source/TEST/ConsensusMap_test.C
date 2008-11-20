// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusMap* ptr = 0;
START_SECTION((ConsensusMap()))
	ptr = new ConsensusMap();
	TEST_NOT_EQUAL(ptr, 0)
	TEST_EQUAL(ptr->isMetaEmpty(),true)
	TEST_REAL_SIMILAR(ptr->getMinInt(), numeric_limits<DoubleReal>::max())
	TEST_REAL_SIMILAR(ptr->getMaxInt(), -numeric_limits<DoubleReal>::max())
END_SECTION

START_SECTION((~ConsensusMap()))
	delete ptr;
END_SECTION

START_SECTION((const String& getIdentifier() const))
	ConsensusMap tmp;
	TEST_EQUAL(tmp.getIdentifier(), "");
END_SECTION

START_SECTION((void setIdentifier(const String& identifier)))
	ConsensusMap tmp;
	tmp.setIdentifier("bla");
	TEST_EQUAL(tmp.getIdentifier(), "bla");
END_SECTION

START_SECTION(const DataProcessing& getDataProcessing() const)
  ConsensusMap tmp;
  TEST_EQUAL(tmp.getDataProcessing().size(),0);
END_SECTION

START_SECTION(DataProcessing& getDataProcessing())
  ConsensusMap tmp;
  tmp.getDataProcessing().resize(1);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION(void setDataProcessing(const DataProcessing& data_processing))
  ConsensusMap tmp;
  std::vector<DataProcessing> dummy;
  dummy.resize(1);
  tmp.setDataProcessing(dummy);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

Feature feature1;
feature1.getPosition()[0] = 2.0;
feature1.getPosition()[1] = 3.0;
feature1.setIntensity(1.0);

Feature feature2;
feature2.getPosition()[0] = 0.0;
feature2.getPosition()[1] = 2.5;
feature2.setIntensity(0.5);

Feature feature3;
feature3.getPosition()[0] = 10.5;
feature3.getPosition()[1] = 0.0;
feature3.setIntensity(0.01);

Feature feature4;
feature4.getPosition()[0] = 5.25;
feature4.getPosition()[1] = 1.5;
feature4.setIntensity(0.5);

START_SECTION(void updateRanges())
  ConsensusMap map;
	ConsensusFeature f;
	f.setIntensity(1.0);
	f.setRT(2.0);
	f.setMZ(3.0);
	f.insert(1,1,feature1);
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
	f.insert(1,2,feature2);
	map.push_back(f);
	map.updateRanges();
	
  TEST_REAL_SIMILAR(map.getMaxInt(),1.0)
  TEST_REAL_SIMILAR(map.getMinInt(),0.5)
  TEST_REAL_SIMILAR(map.getMax()[0],2.0)
  TEST_REAL_SIMILAR(map.getMax()[1],3.0)
  TEST_REAL_SIMILAR(map.getMin()[0],0.0)
  TEST_REAL_SIMILAR(map.getMin()[1],2.5)
  
	//four points
	f.insert(1,3,feature3);
	f.insert(1,4,feature4);
	map.push_back(f);
	map.updateRanges();
	
  TEST_REAL_SIMILAR(map.getMaxInt(),1.0)
  TEST_REAL_SIMILAR(map.getMinInt(),0.01)
  TEST_REAL_SIMILAR(map.getMax()[0],10.5)
  TEST_REAL_SIMILAR(map.getMax()[1],3.0)
  TEST_REAL_SIMILAR(map.getMin()[0],0.0)
  TEST_REAL_SIMILAR(map.getMin()[1],0.0)
  	
END_SECTION

START_SECTION((ConsensusMap& operator = (const ConsensusMap& source)))
  ConsensusMap cons_map;
  cons_map.setMetaValue("meta",String("value"));
  cons_map.setIdentifier("lsid");
  cons_map.getFileDescriptions()[0].filename = "blub";
  cons_map.getFileDescriptions()[0].size = 47;
  cons_map.getFileDescriptions()[0].label = "label";
	cons_map.getFileDescriptions()[0].setMetaValue("meta",String("meta"));
	cons_map.getDataProcessing().resize(1);
	cons_map.setExperimentType("itraq");
  vector< FeatureMap<>* > map_vector(4);
  
  ConsensusMap cons_map_copy;
  cons_map_copy = cons_map;

  TEST_EQUAL(cons_map_copy.getIdentifier(),"lsid")
  TEST_EQUAL(cons_map_copy.getMetaValue("meta").toString(),"value")
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].filename == "blub", true)
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].label == "label", true)
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].size == 47, true)
	TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].getMetaValue("meta") == "meta", true)
  TEST_EQUAL(cons_map_copy.getExperimentType() == "itraq", true)
  TEST_EQUAL(cons_map_copy.getDataProcessing().size(),1)
END_SECTION

START_SECTION((ConsensusMap(const ConsensusMap& source)))
  ConsensusMap cons_map;
  cons_map.setMetaValue("meta",String("value"));
  cons_map.setIdentifier("lsid");
  cons_map.getFileDescriptions()[0].filename = "blub";
  cons_map.getFileDescriptions()[0].size = 47;
  cons_map.getFileDescriptions()[0].label = "label";
	cons_map.getFileDescriptions()[0].setMetaValue("meta",String("meta"));
	cons_map.getDataProcessing().resize(1);
	cons_map.setExperimentType("itraq");
  vector< FeatureMap<>* > map_vector(4);
  
  ConsensusMap cons_map_copy(cons_map);

  TEST_EQUAL(cons_map_copy.getIdentifier(),"lsid")
  TEST_EQUAL(cons_map_copy.getMetaValue("meta").toString(),"value")
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].filename == "blub", true)
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].label == "label", true)
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].size == 47, true)
	TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].getMetaValue("meta") == "meta", true)
  TEST_EQUAL(cons_map_copy.getExperimentType() == "itraq", true)
  TEST_EQUAL(cons_map_copy.getDataProcessing().size(),1)
END_SECTION

START_SECTION((ConsensusMap(Base::size_type n)))
  ConsensusMap cons_map(5);
  
  TEST_REAL_SIMILAR(cons_map.size(),5)
END_SECTION

START_SECTION((const FileDescriptions& getFileDescriptions() const ))
  ConsensusMap cons_map;
  
  TEST_REAL_SIMILAR(cons_map.getFileDescriptions().size(),0)
END_SECTION

START_SECTION((FileDescriptions& getFileDescriptions()))
  ConsensusMap cons_map;
	
  cons_map.getFileDescriptions()[0].filename = "blub";
  TEST_EQUAL(cons_map.getFileDescriptions()[0].filename == "blub", true)
END_SECTION		
		
START_SECTION((const String& getExperimentType() const ))
  ConsensusMap cons_map;
	TEST_EQUAL(cons_map.getExperimentType() == "", true)	
END_SECTION		

START_SECTION((void setExperimentType(const String& experiment_type) ))
  ConsensusMap cons_map;
	cons_map.setExperimentType("itraq");
  TEST_EQUAL(cons_map.getExperimentType() == "itraq", true)
END_SECTION		
		
START_SECTION((bool isValid(String& error_message) const))
	String error_message;
	ConsensusMap cm;
	//empty map
	TEST_EQUAL(cm.isValid(error_message),true)
	//one, valid feature
	ConsensusFeature f;
	f.insert(1,1,Feature());
	cm.push_back(f);
  cm.getFileDescriptions()[1].filename = "bla";
	cm.getFileDescriptions()[1].size = 5;
	TEST_EQUAL(cm.isValid(error_message),true)
	TEST_EQUAL(error_message,"")
	//two valid features
	f.insert(1,4,Feature());
	cm.push_back(f);
	TEST_EQUAL(cm.isValid(error_message),true)
	TEST_EQUAL(error_message,"")
	//two valid and one invalid feature (map index)
	f.insert(2,1,Feature());
	cm.push_back(f);
	TEST_EQUAL(cm.isValid(error_message),false)
	TEST_NOT_EQUAL(error_message,"")
	//one invalid feature (element index)
	cm.clear();
	ConsensusFeature f2;
	f2.insert(2,1,Feature());
	cm.push_back(f2);
	
END_SECTION

START_SECTION(void swap(ConsensusMap& from))
	ConsensusMap map1, map2;	
	ConsensusFeature f;
	f.insert(1,1,Feature());
	map1.push_back(f);
  map1.getFileDescriptions()[1].filename = "bla";
	map1.getFileDescriptions()[1].size = 5;
	map1.setIdentifier("LSID");
	map1.setExperimentType("itraq");
	map1.getDataProcessing().resize(1);
	
	map1.swap(map2);

	TEST_EQUAL(map1.size(),0)
	TEST_EQUAL(map1.getFileDescriptions().size(),0)
	TEST_EQUAL(map1.getIdentifier(),"")	
  TEST_EQUAL(map1.getDataProcessing().size(),0)
	
	TEST_EQUAL(map2.size(),1)
	TEST_EQUAL(map2.getFileDescriptions().size(),1)
	TEST_EQUAL(map2.getIdentifier(),"LSID")	
  TEST_EQUAL(map2.getExperimentType() == "itraq", true)
  TEST_EQUAL(map2.getDataProcessing().size(),1)
END_SECTION

START_SECTION(bool operator == (const ConsensusMap& rhs) const)
	ConsensusMap empty,edit;
	
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
	edit.setMetaValue("bla", 4.1);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.getFileDescriptions()[0].filename = "bla";
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.setExperimentType("bla");
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.push_back(feature1);
	edit.push_back(feature2);
	edit.updateRanges();
	edit.clear();
	TEST_EQUAL(empty==edit, false);
END_SECTION

START_SECTION(bool operator != (const ConsensusMap& rhs) const)
	ConsensusMap empty,edit;
	
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
	edit.setMetaValue("bla", 4.1);
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.getFileDescriptions()[0].filename = "bla";
	TEST_EQUAL(empty!=edit, true)

	edit = empty;
	edit.setExperimentType("bla");
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.push_back(feature1);
	edit.push_back(feature2);
	edit.updateRanges();
	edit.clear();
	TEST_EQUAL(empty!=edit, true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



