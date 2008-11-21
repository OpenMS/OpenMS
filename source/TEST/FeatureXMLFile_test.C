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
// $Maintainer: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
///////////////////////////

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(float a, float b)
{
  DPosition<1> pa(a), pb(b);
  return DRange<1>(pa, pb);
}

///////////////////////////

START_TEST(FeatureXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureXMLFile* ptr = 0;
START_SECTION((FeatureXMLFile()))
	ptr = new FeatureXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~FeatureXMLFile()))
	delete ptr;
END_SECTION
 
START_SECTION((void load(String filename, FeatureMap<>& feature_map)))
	TOLERANCE_ABSOLUTE(0.01)
	
	FeatureMap<> e;
	FeatureXMLFile dfmap_file;
	
	//test exception
	TEST_EXCEPTION( Exception::FileNotFound , dfmap_file.load("dummy/dummy.MzData",e) )
	
	// real test
	dfmap_file.load("data/FeatureXMLFile_1.featureXML",e);
	TEST_EQUAL(e.getIdentifier(),"lsid");
	TEST_EQUAL(e.size(),2)
	TEST_REAL_SIMILAR(e[0].getRT(), 25)
	TEST_REAL_SIMILAR(e[0].getMZ(), 0)
	TEST_REAL_SIMILAR(e[0].getIntensity(), 300)
	TEST_EQUAL(e[0].getMetaValue("stringparametername"),"stringparametervalue")
	TEST_EQUAL((UInt)e[0].getMetaValue("intparametername"),4)
	TEST_REAL_SIMILAR((DoubleReal)e[0].getMetaValue("floatparametername"),4.551)
	TEST_REAL_SIMILAR(e[1].getRT(), 0)
	TEST_REAL_SIMILAR(e[1].getMZ(), 35)
	TEST_REAL_SIMILAR(e[1].getIntensity(), 500)
	//data processing
	TEST_EQUAL(e.getDataProcessing().size(),2)
	TEST_STRING_EQUAL(e.getDataProcessing()[0].getSoftware().getName(),"Software1")
	TEST_STRING_EQUAL(e.getDataProcessing()[0].getSoftware().getVersion(),"0.91a")
	TEST_EQUAL(e.getDataProcessing()[0].getProcessingActions().size(),1)
	TEST_EQUAL(e.getDataProcessing()[0].getProcessingActions().count(DataProcessing::DEISOTOPING),1)
	TEST_STRING_EQUAL(e.getDataProcessing()[0].getMetaValue("name"),"dataProcessing")
	TEST_STRING_EQUAL(e.getDataProcessing()[1].getSoftware().getName(),"Software2")
	TEST_STRING_EQUAL(e.getDataProcessing()[1].getSoftware().getVersion(),"0.92a")
	TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().size(),2)
	TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::SMOOTHING),1)
	TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::BASELINE_REDUCTION),1)

	//test of old file with mzData description (version 1.2)
	//here only the downward-compatibility of the new parser is tested
	//no exception should be thrown
	dfmap_file.load("data/FeatureXMLFile_3_old.featureXML",e);
	TEST_EQUAL(e.size(),1)

	//PeakFileOptions tests
	dfmap_file.getOptions().setRTRange(makeRange(0, 10));
	dfmap_file.load("data/FeatureXMLFile_1.featureXML",e);
	TEST_EQUAL(e.size(),1)
	TEST_REAL_SIMILAR(e[0].getRT(), 0)
	TEST_REAL_SIMILAR(e[0].getMZ(), 35)
	TEST_REAL_SIMILAR(e[0].getIntensity(), 500)

	dfmap_file.getOptions() = PeakFileOptions();
	dfmap_file.getOptions().setMZRange(makeRange(10, 50));
	dfmap_file.load("data/FeatureXMLFile_1.featureXML",e);
	TEST_EQUAL(e.size(),1)
	TEST_REAL_SIMILAR(e[0].getRT(), 0)
	TEST_REAL_SIMILAR(e[0].getMZ(), 35)
	TEST_REAL_SIMILAR(e[0].getIntensity(), 500)

	dfmap_file.getOptions() = PeakFileOptions();
	dfmap_file.getOptions().setIntensityRange(makeRange(400, 600));
	dfmap_file.load("data/FeatureXMLFile_1.featureXML",e);
	TEST_EQUAL(e.size(),1)
	TEST_REAL_SIMILAR(e[0].getRT(), 0)
	TEST_REAL_SIMILAR(e[0].getMZ(), 35)
	TEST_REAL_SIMILAR(e[0].getIntensity(), 500)
	

END_SECTION

START_SECTION((void store(String filename, const FeatureMap<>& feature_map) const))
  
  std::string tmp_filename;
  FeatureMap<> e;
  FeatureXMLFile f;
  
  NEW_TMP_FILE(tmp_filename);
  f.load("data/FeatureXMLFile_1.featureXML",e);
  f.store(tmp_filename,e);
  TEST_FILE_EQUAL(tmp_filename.c_str(),"data/FeatureXMLFile_1.featureXML");

END_SECTION

START_SECTION( PeakFileOptions& getOptions() )
	FeatureXMLFile f;
  FeatureMap<> e;
	f.getOptions().setRTRange(makeRange(1.5, 4.5));
	f.load("data/FeatureXMLFile_2_options.featureXML",e);
	TEST_EQUAL(e.size(), 5)

	f.getOptions().setMZRange(makeRange(1025.0, 2000.0));
	f.load("data/FeatureXMLFile_2_options.featureXML",e);
	TEST_EQUAL(e.size(), 3)

	f.getOptions().setIntensityRange(makeRange(290.0, 310.0));
	f.load("data/FeatureXMLFile_2_options.featureXML",e);
	TEST_EQUAL(e.size(), 1)
	
	f.getOptions().setMetadataOnly(true);
	f.load("data/FeatureXMLFile_2_options.featureXML",e);
	TEST_EQUAL(e.getIdentifier(), "lsid2")
	TEST_EQUAL(e.size(), 0)
END_SECTION

START_SECTION([EXTRA] static bool isValid(const String& filename))
  FeatureXMLFile f;
	TEST_EQUAL(f.isValid("data/FeatureXMLFile_1.featureXML"),true);	
	TEST_EQUAL(f.isValid("data/FeatureXMLFile_2_options.featureXML"),true);	

	FeatureMap<> e;
	String filename;
	
  //test if empty file is valid
	NEW_TMP_FILE(filename)
	f.store(filename,e);	
  TEST_EQUAL(f.isValid(filename),true);	
	
	//test if full file is valid
	NEW_TMP_FILE(filename);
	f.load("data/FeatureXMLFile_1.featureXML",e);
	f.store(filename, e);	
  TEST_EQUAL(f.isValid(filename),true);
END_SECTION

START_SECTION( const PeakFileOptions& getOptions() const )
	FeatureXMLFile f;
 	FeatureMap<> e;
	f.getOptions().setRTRange(makeRange(1.5, 4.5));
	f.getOptions().setIntensityRange(makeRange(290.0, 310.0));
	
	const PeakFileOptions pfo = f.getOptions();
	
	TEST_EQUAL(pfo.getRTRange(),makeRange(1.5, 4.5))
	TEST_EQUAL(pfo.getIntensityRange(),makeRange(290.0, 310.0))
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
	
	FeatureMap<> map;
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

	Feature element;
	element.setMZ(736.53445);
	element.setRT(66324.47324);
	element.setIntensity(123123.321);
	
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

		element.getPeptideIdentifications().push_back(hits);
	}

	map.push_back(element);

	FeatureXMLFile file;

	String tmp_filename;
  NEW_TMP_FILE(tmp_filename);

	file.store(tmp_filename,map);
	
	FeatureMap<> map_reloaded;
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
 
