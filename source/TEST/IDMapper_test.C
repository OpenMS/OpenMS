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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

///////////////////////////

START_TEST(IDMapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


IDMapper* ptr = 0;
START_SECTION((IDMapper()))
	ptr = new IDMapper();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~IDMapper()))
	delete ptr;
END_SECTION

START_SECTION((DoubleReal getRTDelta() const))
	IDMapper mapper;
	TEST_REAL_SIMILAR(mapper.getRTDelta(),5.0)
END_SECTION

START_SECTION((void setRTDelta(DoubleReal rt_delta)))
	IDMapper mapper;
	mapper.setRTDelta(1.5);
	TEST_REAL_SIMILAR(mapper.getRTDelta(),1.5)
END_SECTION

START_SECTION((DoubleReal getMZDelta() const))
	IDMapper mapper;
	TEST_REAL_SIMILAR(mapper.getMZDelta(),0.01)
END_SECTION

START_SECTION((void setMZDelta(DoubleReal mz_delta)))
	IDMapper mapper;
	mapper.setMZDelta(1.05);
	TEST_REAL_SIMILAR(mapper.getMZDelta(),1.05)
END_SECTION


START_SECTION((template <typename PeakType> void annotate(MSExperiment< PeakType > &map, const std::vector< PeptideIdentification > &ids, const std::vector< ProteinIdentification > &protein_ids)))
	//load id
	vector<PeptideIdentification> identifications; 
	vector<ProteinIdentification> protein_identifications;
	String document_id;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_1.idXML"), protein_identifications, identifications, document_id);
	
	TEST_EQUAL(identifications.size(),3)
	TEST_EQUAL(identifications[0].getHits().size(), 2)
	TEST_EQUAL(identifications[1].getHits().size(), 1)
	TEST_EQUAL(identifications[2].getHits().size(), 2)
	TEST_EQUAL(protein_identifications.size(),1)
	TEST_EQUAL(protein_identifications[0].getHits().size(), 2)

	//create experiment
	MSExperiment<> experiment;
	MSSpectrum<> spectrum;
	Precursor precursor;
	precursor.setMZ(0);
	spectrum.setRT(60);
	experiment.push_back(spectrum);							
	experiment[0].getPrecursors().push_back(precursor);
	precursor.setMZ(20);
	spectrum.setRT(181);
	experiment.push_back(spectrum);							
	experiment[1].getPrecursors().push_back(precursor);
	precursor.setMZ(11);
	spectrum.setRT(120.0001);
	experiment.push_back(spectrum);							
	experiment[2].getPrecursors().push_back(precursor);
	
	//map
	IDMapper mapper;
	mapper.setRTDelta(0.5);
	mapper.setMZDelta(0.05);
	mapper.annotate(experiment, identifications, protein_identifications);
	
	//test
	TEST_EQUAL(experiment.getProteinIdentifications().size(), 1)
	TEST_EQUAL(experiment.getProteinIdentifications()[0].getHits().size(),2)
	TEST_EQUAL(experiment.getProteinIdentifications()[0].getHits()[0].getAccession(),"ABCDE")
	TEST_EQUAL(experiment.getProteinIdentifications()[0].getHits()[1].getAccession(),"FGHIJ")
	//scan 1
	TEST_EQUAL(experiment[0].getPeptideIdentifications().size(), 1)
	TEST_EQUAL(experiment[0].getPeptideIdentifications()[0].getHits().size(), 2)
	TEST_EQUAL(experiment[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(experiment[0].getPeptideIdentifications()[0].getHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	//scan 2
	TEST_EQUAL(experiment[1].getPeptideIdentifications().size(), 0)
	//scan 3
	TEST_EQUAL(experiment[2].getPeptideIdentifications().size(), 1)	
	TEST_EQUAL(experiment[2].getPeptideIdentifications()[0].getHits().size(), 1)
	TEST_EQUAL(experiment[2].getPeptideIdentifications()[0].getHits()[0].getSequence(), "HSKLSAK")
END_SECTION



START_SECTION((template <typename FeatureType> void annotate(FeatureMap<FeatureType> &map, const std::vector< PeptideIdentification > &ids, const std::vector< ProteinIdentification >& protein_ids, bool use_centroids=false)))
	//load id data
	vector<PeptideIdentification> identifications; 
	vector<ProteinIdentification> protein_identifications;
	String document_id;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_2.idXML"), protein_identifications, identifications, document_id);
	
	//--------------------------------------------------------------------------------------
	//TEST MAPPING TO CONVEX HULLS
	FeatureMap<> fm;
	FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_2.featureXML"), fm);
	
	IDMapper mapper;
	mapper.setRTDelta(0.0);
	mapper.setMZDelta(0.0);
	mapper.annotate(fm,identifications,protein_identifications);
	
	//test protein ids
	TEST_EQUAL(fm.getProteinIdentifications().size(),1)
	TEST_EQUAL(fm.getProteinIdentifications()[0].getHits().size(),2)
	TEST_EQUAL(fm.getProteinIdentifications()[0].getHits()[0].getAccession(),"ABCDE")
	TEST_EQUAL(fm.getProteinIdentifications()[0].getHits()[1].getAccession(),"FGHIJ")
	
	//test peptide ids
	TEST_EQUAL(fm[0].getPeptideIdentifications().size(),7)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[1].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[2].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[3].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[4].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[5].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[6].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getHits()[0].getSequence(),"A")
	TEST_EQUAL(fm[0].getPeptideIdentifications()[1].getHits()[0].getSequence(),"K")
	TEST_EQUAL(fm[0].getPeptideIdentifications()[2].getHits()[0].getSequence(),"C")
	TEST_EQUAL(fm[0].getPeptideIdentifications()[3].getHits()[0].getSequence(),"D")
	TEST_EQUAL(fm[0].getPeptideIdentifications()[4].getHits()[0].getSequence(),"E")
	TEST_EQUAL(fm[0].getPeptideIdentifications()[5].getHits()[0].getSequence(),"F")
	TEST_EQUAL(fm[0].getPeptideIdentifications()[6].getHits()[0].getSequence(),"I")
	
	//test unassigned peptide ids
	TEST_EQUAL(fm.getUnassignedPeptideIdentifications().size(),3)
	TEST_EQUAL(fm.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence(),"G")
	TEST_EQUAL(fm.getUnassignedPeptideIdentifications()[1].getHits()[0].getSequence(),"H")
	TEST_EQUAL(fm.getUnassignedPeptideIdentifications()[2].getHits()[0].getSequence(),"L")
		
	//--------------------------------------------------------------------------------------
	//TEST MAPPING TO CENTROIDS
	FeatureMap<> fm2;
	FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_2.featureXML"), fm2);
	mapper.setRTDelta(4.0);
	mapper.setMZDelta(1.5);
	mapper.annotate(fm2,identifications,protein_identifications, true);
	
	//test protein ids
	TEST_EQUAL(fm2.getProteinIdentifications().size(),1)
	TEST_EQUAL(fm2.getProteinIdentifications()[0].getHits().size(),2)
	TEST_EQUAL(fm2.getProteinIdentifications()[0].getHits()[0].getAccession(),"ABCDE")
	TEST_EQUAL(fm2.getProteinIdentifications()[0].getHits()[1].getAccession(),"FGHIJ")
	
	//test peptide ids
	TEST_EQUAL(fm2[0].getPeptideIdentifications().size(),2)
	TEST_EQUAL(fm2[0].getPeptideIdentifications()[0].getHits().size(),1)
	TEST_EQUAL(fm2[0].getPeptideIdentifications()[1].getHits().size(),1)
	TEST_EQUAL(fm2[0].getPeptideIdentifications()[0].getHits()[0].getSequence(),"A")
	TEST_EQUAL(fm2[0].getPeptideIdentifications()[1].getHits()[0].getSequence(),"K")

	//test unassigned peptide ids
	TEST_EQUAL(fm2.getUnassignedPeptideIdentifications().size(),8)
	TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence(),"C")
	TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[1].getHits()[0].getSequence(),"D")
	TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[2].getHits()[0].getSequence(),"E")
	TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[3].getHits()[0].getSequence(),"F")
	TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[4].getHits()[0].getSequence(),"G")
	TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[5].getHits()[0].getSequence(),"H")
	TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[6].getHits()[0].getSequence(),"I")
	TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[7].getHits()[0].getSequence(),"L")
END_SECTION


START_SECTION((void annotate(ConsensusMap& map, const std::vector<PeptideIdentification>& ids, const std::vector<ProteinIdentification>& protein_ids, bool measure_from_subelements=false)))
{
	IDMapper mapper;
	TOLERANCE_ABSOLUTE(0.01);
		
	std::vector<ProteinIdentification> protein_ids;
	std::vector<PeptideIdentification> peptide_ids;
	String document_id;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_3.idXML"), protein_ids, peptide_ids, document_id);
	
  ConsensusXMLFile cons_file;
  
	{
	  std::string tmp_filename;
	  NEW_TMP_FILE(tmp_filename);
		ConsensusMap cons_map;
		cons_file.load(OPENMS_GET_TEST_DATA_PATH("IDMapper_3.consensusXML"), cons_map);
		mapper.annotate(cons_map, peptide_ids, protein_ids);
		cons_file.store(tmp_filename,cons_map);
		WHITELIST("<?xml-stylesheet");
		TEST_FILE_SIMILAR(tmp_filename,OPENMS_GET_TEST_DATA_PATH("IDMapper_3_out1.consensusXML"));
	}

	{
	  std::string tmp_filename;
	  NEW_TMP_FILE(tmp_filename);
		ConsensusMap cons_map;
		cons_file.load(OPENMS_GET_TEST_DATA_PATH("IDMapper_3.consensusXML"), cons_map);
		mapper.annotate(cons_map, peptide_ids, protein_ids, true);
		cons_file.store(tmp_filename,cons_map);
		WHITELIST("<?xml-stylesheet");
		TEST_FILE_SIMILAR(tmp_filename,OPENMS_GET_TEST_DATA_PATH("IDMapper_3_out2.consensusXML"));
	}
	
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
