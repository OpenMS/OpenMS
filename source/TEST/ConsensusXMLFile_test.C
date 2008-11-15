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
// $Maintainer: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;


DRange<1> makeRange(float a, float b)
{
  DPosition<1> pa(a), pb(b);
  return DRange<1>(pa, pb);
}

START_TEST(ConsensusXMLFile, "$Id$")

#if 0

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusXMLFile* ptr = 0;
CHECK((ConsensusXMLFile()))
	ptr = new ConsensusXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~ConsensusXMLFile()))
	delete ptr;
RESULT

PRECISION(0.01)

CHECK((void load(const String &filename, ConsensusMap &map)))
  ConsensusMap cons_map;
  ConsensusXMLFile cons_file;
  cons_file.load("data/ConsensusXMLFile.consensusXML", cons_map);
	TEST_EQUAL(cons_map.getExperimentType()=="label-free",true)
	TEST_EQUAL(cons_map.getMetaValue("name1")==DataValue("value1"),true)
	TEST_EQUAL(cons_map.getMetaValue("name2")==DataValue(2),true)
  TEST_EQUAL(cons_map.getIdentifier(),"lsid");
  TEST_EQUAL(cons_map.getFileDescriptions()[0].filename == "data/MapAlignmentFeatureMap1.xml", true)
  TEST_EQUAL(cons_map.getFileDescriptions()[0].label,"label")
  TEST_EQUAL(cons_map.getFileDescriptions()[0].size, 144)
	TEST_EQUAL(cons_map.getFileDescriptions()[0].getMetaValue("name3")==DataValue("value3"),true)
	TEST_EQUAL(cons_map.getFileDescriptions()[0].getMetaValue("name4")==DataValue(4),true)
  TEST_STRING_EQUAL(cons_map.getFileDescriptions()[1].filename,"data/MapAlignmentFeatureMap2.xml")
  TEST_EQUAL(cons_map.getFileDescriptions()[1].label,"")
  TEST_EQUAL(cons_map.getFileDescriptions()[1].size, 0)
	TEST_EQUAL(cons_map.getFileDescriptions()[1].getMetaValue("name5")==DataValue("value5"),true)
	TEST_EQUAL(cons_map.getFileDescriptions()[1].getMetaValue("name6")==DataValue(6.0),true)

  ConsensusFeature cons_feature = cons_map[0];
  TEST_REAL_EQUAL(cons_feature.getRT(),1273.27)  
  TEST_REAL_EQUAL(cons_feature.getMZ(),904.47)
  TEST_REAL_EQUAL(cons_feature.getIntensity(),3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[0],1273.27)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[0],1273.27)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().min()[0],3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().max()[0],3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getQuality(),1.1)
  TEST_EQUAL(cons_feature.getMetaValue("peptide_id")==DataValue("RefSeq:NC_1234"),true)
  ConsensusFeature::HandleSetType::const_iterator it = cons_feature.begin();
  TEST_REAL_EQUAL(it->getIntensity(),3.12539e+07)
  
  cons_feature = cons_map[5];
  TEST_REAL_EQUAL(cons_feature.getRT(),1194.82)  
  TEST_REAL_EQUAL(cons_feature.getMZ(),777.101)
  TEST_REAL_EQUAL(cons_feature.getIntensity(),1.78215e+07)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[0],1194.82)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[0],1194.82)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[1],777.101)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[1],777.101)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().min()[0],1.78215e+07)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().max()[0],1.78215e+07)
  TEST_REAL_EQUAL(cons_feature.getQuality(),0.0)
  it = cons_feature.begin();
  TEST_REAL_EQUAL(it->getIntensity(),1.78215e+07)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),1.78215e+07)


	//PeakFileOptions tests
	
	cons_file.getOptions().setRTRange(makeRange(815, 818));
	cons_file.load("data/ConsensusXMLFile3.consensusXML",cons_map);
	TEST_EQUAL(cons_map.size(),1)
	TEST_REAL_EQUAL(cons_map[0].getRT(), 817.266)

	cons_file.getOptions() = PeakFileOptions();
	cons_file.getOptions().setMZRange(makeRange(944, 945));
	cons_file.load("data/ConsensusXMLFile3.consensusXML",cons_map);
	TEST_EQUAL(cons_map.size(),1)
	TEST_REAL_EQUAL(cons_map[0].getMZ(), 944.96)

	cons_file.getOptions() = PeakFileOptions();
	cons_file.getOptions().setIntensityRange(makeRange(15000,24000));
	cons_file.load("data/ConsensusXMLFile3.consensusXML",cons_map);
	TEST_EQUAL(cons_map.size(),1)
	TEST_REAL_EQUAL(cons_map[0].getIntensity(),23000.238)

RESULT

CHECK((void store(const String &filename, const ConsensusMap &map)))
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  
  ConsensusMap cons_map, cons_map2;
  ConsensusXMLFile cons_file;
  
  cons_file.load("data/ConsensusXMLFile2.consensusXML",cons_map);  
  cons_file.store(tmp_filename,cons_map);
  cons_file.load(tmp_filename,cons_map2);  
  TEST_EQUAL(cons_map.size(),cons_map2.size())
  TEST_EQUAL(cons_map[0]==cons_map2[0],true)
  TEST_EQUAL(cons_map[1]==cons_map2[1],true)
RESULT

CHECK([EXTRA] (bool isValid(const String& filename)))
  ConsensusXMLFile cons_file;
  TEST_EQUAL(cons_file.isValid("data/ConsensusXMLFile.consensusXML"),true);
  TEST_EQUAL(cons_file.isValid("data/ConsensusXMLFile2.consensusXML"),true);
  TEST_EQUAL(cons_file.isValid("data/ConsensusXMLFile3.consensusXML"),true);
RESULT

CHECK([EXTRA] ProteinIdentification PeptideIdentification I/O)
{
	std::vector<ProteinHit> protein_hits;
	ProteinHit protein_hit;

	protein_hit.setAccession("urn:lsid:ach0du1schreck2wo3iss4er5denn");
	protein_hit.setMetaValue("dadada",String("dududu"));
	protein_hit.setRank(768);
	protein_hit.setScore(70.3);
	protein_hit.setSequence("ABCDEFG");

	protein_hits.push_back(protein_hit);

	protein_hit.setAccession("urn:lsid:rumpelstielzchen");
	protein_hit.setMetaValue("dadada",String("doppeltsogut"));
	protein_hit.setRank(543);
	protein_hit.setScore(140.6);
	protein_hit.setSequence("HIJKLMN");

	protein_hits.push_back(protein_hit);

	ProteinIdentification protein_identification;

	ConsensusMap cons_map;
	
	{
		ProteinIdentification hits;
		hits.setDateTime(DateTime::now());
		hits.setSignificanceThreshold(56.7643);
		hits.insertHit(protein_hits[0]);
		hits.insertHit(protein_hits[1]);
		hits.setIdentifier("id");
		hits.setScoreType("score_type");
		hits.setHigherScoreBetter(true);
		hits.setSearchEngine("MaxKotzt");
		hits.setSearchEngineVersion("2.1");
		ProteinIdentification::SearchParameters param;
		param.db = "RefSeq";
		hits.setSearchParameters(param);
		hits.setMetaValue("pi",3.14159);

		cons_map.getProteinIdentifications().push_back(hits);
	}

	ConsensusMap::FileDescription file_description;
	file_description.label = "dummy_label";
	file_description.size = 19;
	file_description.filename = "dummy/filename";

	cons_map.getFileDescriptions()[222] = file_description;
	
	Peak2D peak2d;
	peak2d.setMZ(736.53445);
	peak2d.setRT(66324.47324);
	peak2d.setIntensity(123123.321);
	FeatureHandle feat_handle(222,9,peak2d);

	ConsensusFeature cons_feat;
	cons_feat.insert(feat_handle);

	float peptide_significance_threshold = 42.3;
	std::vector<PeptideHit> peptide_hits;

	PeptideHit peptide_hit;
	peptide_hit.setProteinAccessions(std::vector<String>(1,"urn:lsid:rumpelstielzchen"));
	peptide_hit.setScore(4324.433);
	peptide_hit.setSequence("HAL");
	peptide_hit.setCharge(23);
	peptide_hit.setAABefore('X');
	peptide_hit.setAAAfter('Y');

	peptide_hits.push_back(peptide_hit);

	{
		PeptideIdentification hits;

		hits.setHits(peptide_hits);

		hits.setHigherScoreBetter(false);
		hits.setIdentifier("id");
		hits.setMetaValue("label",17);
		hits.setScoreType("score_type");
		hits.setSignificanceThreshold(peptide_significance_threshold);

		cons_feat.getPeptideIdentifications().push_back(hits);
	}

	cons_map.push_back(cons_feat);

	ConsensusXMLFile cons_file;

	String tmp_filename;
  NEW_TMP_FILE(tmp_filename);

	cons_file.store(tmp_filename,cons_map);
	
	ConsensusMap cons_map_reloaded;
	cons_file.load(tmp_filename,cons_map_reloaded);

	String tmp_filename_reloaded_then_stored;
	NEW_TMP_FILE(tmp_filename_reloaded_then_stored);

	cons_file.store(tmp_filename_reloaded_then_stored,cons_map_reloaded);

	TEST_FILE(tmp_filename.c_str(),tmp_filename_reloaded_then_stored.c_str());

	// STATUS("Will intentionally fail on next line so we can have a look at the tmp files");
	// TEST_EQUAL(0,1);
}
RESULT
#endif

CHECK([EXTRA] load/store)
{
	ConsensusXMLFile cm_file;
	ConsensusMap cm;
	cm_file.load("data/ItraqQuantifier.consensusXML",cm);
    
	String cm_file_out;// = "data/ItraqQuantifier.consensusXML";
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm);
    
	FuzzyStringComparator fsc;
	fsc.setAcceptableAbsolute(0.01);
	TEST_EQUAL(fsc.compare_files(cm_file_out,"data/ItraqQuantifier.consensusXML"), true);
}
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



 
