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
///////////////////////////

using namespace OpenMS;
using namespace std;


DRange<1> makeRange(float a, float b)
{
  DPosition<1> pa(a), pb(b);
  return DRange<1>(pa, pb);
}

START_TEST(ConsensusXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusXMLFile* ptr = 0;
START_SECTION((ConsensusXMLFile()))
	ptr = new ConsensusXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~ConsensusXMLFile()))
	delete ptr;
END_SECTION

TOLERANCE_ABSOLUTE(0.01)

START_SECTION((void load(const String &filename, ConsensusMap &map)))
  ConsensusMap map;
  ConsensusXMLFile file;
  file.load("data/ConsensusXMLFile_1.consensusXML", map);
	//meta data
  TEST_EQUAL(map.getIdentifier(),"lsid")
	TEST_EQUAL(map.getExperimentType()=="label-free",true)
	TEST_EQUAL(map.getMetaValue("name1")==DataValue("value1"),true)
	TEST_EQUAL(map.getMetaValue("name2")==DataValue(2),true)
	//file descriptions
  TEST_EQUAL(map.getFileDescriptions()[0].filename == "data/MapAlignmentFeatureMap1.xml", true)
  TEST_EQUAL(map.getFileDescriptions()[0].label,"label")
  TEST_EQUAL(map.getFileDescriptions()[0].size, 144)
	TEST_EQUAL(map.getFileDescriptions()[0].getMetaValue("name3")==DataValue("value3"),true)
	TEST_EQUAL(map.getFileDescriptions()[0].getMetaValue("name4")==DataValue(4),true)
  TEST_STRING_EQUAL(map.getFileDescriptions()[1].filename,"data/MapAlignmentFeatureMap2.xml")
  TEST_EQUAL(map.getFileDescriptions()[1].label,"")
  TEST_EQUAL(map.getFileDescriptions()[1].size, 0)
	TEST_EQUAL(map.getFileDescriptions()[1].getMetaValue("name5")==DataValue("value5"),true)
	TEST_EQUAL(map.getFileDescriptions()[1].getMetaValue("name6")==DataValue(6.0),true)
	//data processing
	TEST_EQUAL(map.getDataProcessing().size(),2)
	TEST_STRING_EQUAL(map.getDataProcessing()[0].getSoftware().getName(),"Software1")
	TEST_STRING_EQUAL(map.getDataProcessing()[0].getSoftware().getVersion(),"0.91a")
	TEST_EQUAL(map.getDataProcessing()[0].getProcessingActions().size(),1)
	TEST_EQUAL(map.getDataProcessing()[0].getProcessingActions().count(DataProcessing::DEISOTOPING),1)
	TEST_STRING_EQUAL(map.getDataProcessing()[0].getMetaValue("name"),"dataProcessing")
	TEST_STRING_EQUAL(map.getDataProcessing()[1].getSoftware().getName(),"Software2")
	TEST_STRING_EQUAL(map.getDataProcessing()[1].getSoftware().getVersion(),"0.92a")
	TEST_EQUAL(map.getDataProcessing()[1].getProcessingActions().size(),2)
	TEST_EQUAL(map.getDataProcessing()[1].getProcessingActions().count(DataProcessing::SMOOTHING),1)
	TEST_EQUAL(map.getDataProcessing()[1].getProcessingActions().count(DataProcessing::BASELINE_REDUCTION),1)

	//features
  ConsensusFeature cons_feature = map[0];
  TEST_REAL_SIMILAR(cons_feature.getRT(),1273.27)  
  TEST_REAL_SIMILAR(cons_feature.getMZ(),904.47)
  TEST_REAL_SIMILAR(cons_feature.getIntensity(),3.12539e+07)
  TEST_REAL_SIMILAR(cons_feature.getPositionRange().min()[0],1273.27)
  TEST_REAL_SIMILAR(cons_feature.getPositionRange().max()[0],1273.27)
  TEST_REAL_SIMILAR(cons_feature.getPositionRange().min()[1],904.47)
  TEST_REAL_SIMILAR(cons_feature.getPositionRange().max()[1],904.47)
  TEST_REAL_SIMILAR(cons_feature.getIntensityRange().min()[0],3.12539e+07)
  TEST_REAL_SIMILAR(cons_feature.getIntensityRange().max()[0],3.12539e+07)
  TEST_REAL_SIMILAR(cons_feature.getQuality(),1.1)
  TEST_EQUAL(cons_feature.getMetaValue("peptide_id")==DataValue("RefSeq:NC_1234"),true)
  ConsensusFeature::HandleSetType::const_iterator it = cons_feature.begin();
  TEST_REAL_SIMILAR(it->getIntensity(),3.12539e+07)
  
  cons_feature = map[5];
  TEST_REAL_SIMILAR(cons_feature.getRT(),1194.82)  
  TEST_REAL_SIMILAR(cons_feature.getMZ(),777.101)
  TEST_REAL_SIMILAR(cons_feature.getIntensity(),1.78215e+07)
  TEST_REAL_SIMILAR(cons_feature.getPositionRange().min()[0],1194.82)
  TEST_REAL_SIMILAR(cons_feature.getPositionRange().max()[0],1194.82)
  TEST_REAL_SIMILAR(cons_feature.getPositionRange().min()[1],777.101)
  TEST_REAL_SIMILAR(cons_feature.getPositionRange().max()[1],777.101)
  TEST_REAL_SIMILAR(cons_feature.getIntensityRange().min()[0],1.78215e+07)
  TEST_REAL_SIMILAR(cons_feature.getIntensityRange().max()[0],1.78215e+07)
  TEST_REAL_SIMILAR(cons_feature.getQuality(),0.0)
  it = cons_feature.begin();
  TEST_REAL_SIMILAR(it->getIntensity(),1.78215e+07)
  ++it;
  TEST_REAL_SIMILAR(it->getIntensity(),1.78215e+07)


	//PeakFileOptions tests
	
	file.getOptions().setRTRange(makeRange(815, 818));
	file.load("data/ConsensusXMLFile_2_options.consensusXML",map);
	TEST_EQUAL(map.size(),1)
	TEST_REAL_SIMILAR(map[0].getRT(), 817.266)

	file.getOptions() = PeakFileOptions();
	file.getOptions().setMZRange(makeRange(944, 945));
	file.load("data/ConsensusXMLFile_2_options.consensusXML",map);
	TEST_EQUAL(map.size(),1)
	TEST_REAL_SIMILAR(map[0].getMZ(), 944.96)

	file.getOptions() = PeakFileOptions();
	file.getOptions().setIntensityRange(makeRange(15000,24000));
	file.load("data/ConsensusXMLFile_2_options.consensusXML",map);
	TEST_EQUAL(map.size(),1)
	TEST_REAL_SIMILAR(map[0].getIntensity(),23000.238)

END_SECTION

START_SECTION((void store(const String &filename, const ConsensusMap &map)))
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  
  ConsensusMap map, map2;
  ConsensusXMLFile f;
  
  f.load("data/ConsensusXMLFile_1.consensusXML",map);  
  f.store(tmp_filename, map);
  f.load(tmp_filename, map2);
  TEST_EQUAL(map==map2, true)
END_SECTION

START_SECTION([EXTRA] (bool isValid(const String& filename)))
  ConsensusXMLFile f;
  TEST_EQUAL(f.isValid("data/ConsensusXMLFile_1.consensusXML"),true);
  TEST_EQUAL(f.isValid("data/ConsensusXMLFile_2_options.consensusXML"),true);
	
  //test if written empty file is invalid, so this is not tested :)
	
	//test if written full file is valid
	ConsensusMap m;
	String tmp_filename;
	NEW_TMP_FILE(tmp_filename);
	f.load("data/ConsensusXMLFile_1.consensusXML",m);
	f.store(tmp_filename, m);	
  TEST_EQUAL(f.isValid(tmp_filename),true);
END_SECTION

START_SECTION([EXTRA] peptide and protein identification I/O)
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

	ProteinIdentification::SearchParameters search_param;
	search_param.db="DB";
	search_param.db_version="1.2.3";
	search_param.taxonomy="wolpertinger";
	search_param.charges="+1,+2,-3";
	search_param.mass_type=ProteinIdentification::MONOISOTOPIC;
	search_param.fixed_modifications.push_back("bli");
	search_param.variable_modifications.push_back("bla");
	search_param.variable_modifications.push_back("bluff");
	search_param.enzyme=ProteinIdentification::TRYPSIN;
	search_param.missed_cleavages=2;
	search_param.peak_mass_tolerance=0.3;
	search_param.precursor_tolerance=0.1;
	
	ProteinIdentification protein_identification;
	protein_identification.setSearchParameters(search_param);

	ConsensusMap map;
	
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

		map.getProteinIdentifications().push_back(hits);
	}

	ConsensusMap::FileDescription file_description;
	file_description.label = "dummy_label";
	file_description.size = 19;
	file_description.filename = "dummy/filename";

	map.getFileDescriptions()[222] = file_description;
	
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

	map.push_back(cons_feat);

	ConsensusXMLFile file;

	String tmp_filename;
  NEW_TMP_FILE(tmp_filename);

	file.store(tmp_filename,map);
	
	ConsensusMap map_reloaded;
	file.load(tmp_filename,map_reloaded);

	String tmp_filename_reloaded_then_stored;
	NEW_TMP_FILE(tmp_filename_reloaded_then_stored);

	file.store(tmp_filename_reloaded_then_stored,map_reloaded);

	TEST_FILE_EQUAL(tmp_filename.c_str(),tmp_filename_reloaded_then_stored.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



 
