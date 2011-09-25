// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Frederic Lehnert, Fabian Kriegel $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/SIMULATION/DigestSimulation.h>
///////////////////////////
#include <OpenMS/SIMULATION/LABELING/ICPLLabeler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

void createTestFeatureMapSimVector_(FeatureMapSimVector& feature_maps, bool add3rd)
{
	feature_maps.clear();

	FeatureMapSim fm1,fm2,fm3;
	ProteinHit prothit1,prothit2,prothit3,prothit4,prothit5,prothit6,prothit7,prothit8,prothit9,prothit10, prothit11, prothit12;

	// create first map
	prothit1.setSequence("AAAAAAAKAAAAA"); // 2 Fragmente AAAAAAAK und AAAAA und kommt in allen Channels vor
	prothit1.setMetaValue("description", "test sequence 1");
	prothit1.setAccession("ACC1");
	prothit1.setMetaValue("intensity", 200.0);

	prothit2.setSequence("CNARCNCNCN"); // 2 Fragmente CNAR und CNCNCN und kommt in allen Channels vor
	prothit2.setMetaValue("description", "test sequence 2");
	prothit2.setAccession("ACC2");
	prothit2.setMetaValue("intensity", 80.0);

	prothit3.setSequence("CNHAADDAAAAA"); // ungelabelt, einzelnes Fragment
	prothit3.setMetaValue("description", "test sequence 3");
	prothit3.setAccession("ACC3");
	prothit3.setMetaValue("intensity", 100.0);

	prothit4.setSequence("VNAAAAAARVNCNCNAAAA"); // Ergebnis: CNAAAAAAR(Label Medium_R) , CNCNCNAAAA (einmal kommt in allen Channels vor)
	prothit4.setMetaValue("description", "test sequence 4");
	prothit4.setAccession("ACC5");
	prothit4.setMetaValue("intensity", 115.0);

	ProteinIdentification protIdent1;
	protIdent1.insertHit(prothit1);
	protIdent1.insertHit(prothit2);
	protIdent1.insertHit(prothit3);
	protIdent1.insertHit(prothit4);
	vector<ProteinIdentification> protIdents_vec1;
	protIdents_vec1.push_back(protIdent1);
	fm1.setProteinIdentifications(protIdents_vec1);

	// create second map
	prothit5.setSequence("AAAAAAAKAAAAA"); // Ergbenis: AAAAAAAK(Label Medium_K) , AAAAA ( einmal kommt in allen Channels vor)
	prothit5.setMetaValue("description", "test sequence 5");
	prothit5.setAccession("ACC4");
	prothit5.setMetaValue("intensity", 50.0);

	prothit6.setSequence("CNARCNCNCN"); // Ergebnis: CNAR(Label Medium_R) , CNCNCN (einmal kommt in allen Channels vor)
	prothit6.setMetaValue("description", "test sequence 6");
	prothit6.setAccession("ACC5");
	prothit6.setMetaValue("intensity", 100.0);

	prothit7.setSequence("LDRCEL"); // Ergbenis : LDR(label Medium_R) , CEL (einmal kommt in channel 2 und 3 vor)
	prothit7.setMetaValue("description", "test sequence 7");
	prothit7.setAccession("ACC6");
	prothit7.setMetaValue("intensity", 120.0);

	prothit8.setSequence("VNAAAAAARVNCNCNAAAA"); // Ergebnis: CNAAAAAAR(Label Medium_R) , CNCNCNAAAA (einmal kommt in allen Channels vor)
	prothit8.setMetaValue("description", "test sequence 8");
	prothit8.setAccession("ACC5");
	prothit8.setMetaValue("intensity", 110.0);

	ProteinIdentification protIdent2;
	protIdent2.insertHit(prothit5);
	protIdent2.insertHit(prothit6);
	protIdent2.insertHit(prothit7);
	protIdent2.insertHit(prothit8);
	vector<ProteinIdentification> protIdents_vec2;
	protIdents_vec2.push_back(protIdent2);
	fm2.setProteinIdentifications(protIdents_vec2);

	feature_maps.push_back(fm1);
	feature_maps.push_back(fm2);

	if (add3rd)
	{
		prothit9.setSequence("AAAAAAAKAAAAA"); // Ergebnis : AAAAAAAK(Label Heavy_K) , AAAAA ( einmal kommt in allen Channels vor )
		prothit9.setMetaValue("description", "test sequence 9");
		prothit9.setAccession("ACC7");
		prothit9.setMetaValue("intensity", 30.0);

		prothit10.setSequence("CNARCNCNCN"); // Ergebnis: CNAR(Label Heavy_R) , CNCNCN (einmal kommt in allen Channels vor)
		prothit10.setMetaValue("description", "test sequence 10");
		prothit10.setAccession("ACC8");
		prothit10.setMetaValue("intensity", 130.0);

		prothit11.setSequence("LDRCEL"); //Ergebnis: LDR(label Heavy_R) , CEL (einmal kommt in channel 2 und 3 vor)
		prothit11.setMetaValue("description", "test sequence 11");
		prothit11.setAccession("ACC9");
		prothit11.setMetaValue("intensity", 70.0);

		prothit12.setSequence("YCYCY"); //Ergebnis: YCYCY kommt nur in diesem Channel vor
		prothit12.setMetaValue("description", "test sequence 12");
		prothit12.setAccession("ACC10");
		prothit12.setMetaValue("intensity", 80.0);

		ProteinIdentification protIdent3;
		protIdent3.insertHit(prothit9);
		protIdent3.insertHit(prothit10);
		protIdent3.insertHit(prothit11);
		protIdent3.insertHit(prothit12);
		vector<ProteinIdentification> protIdents_vec3;
		protIdents_vec3.push_back(protIdent3);
		fm3.setProteinIdentifications(protIdents_vec3);
		feature_maps.push_back(fm3);
	}
}


void digestFeaturesMapSimVector_(FeatureMapSimVector& feature_maps)
{
	// digest here
	DigestSimulation digest_sim;
	Param p;
	p.setValue("model", "naive");
	p.setValue("model_naive:missed_cleavages", 0);
	digest_sim.setParameters(p);
	std::cout << digest_sim.getParameters() << std::endl;
	for(FeatureMapSimVector::iterator iter = feature_maps.begin() ; iter != feature_maps.end() ; ++iter)
	{
		digest_sim.digest((*iter));
	}
}

START_TEST(ICPLLabeler, "$Id: ICPLLabeler_test.C 7837 2011-05-14 11:41:44Z flehnert $")
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ICPLLabeler* ptr = 0;
START_SECTION(ICPLLabeler())
{
	ptr = new ICPLLabeler();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~ICPLLabeler())
{
	delete ptr;
}
END_SECTION

START_SECTION((void setUpHook(FeatureMapSimVector &)))
{
	std::cerr << "*********************************setUpHook*************************************" ;
	ICPLLabeler labeler;

	FeatureMapSim fm1,fm2,fm3,fm4;
	FeatureMapSimVector fm_vec;

	fm_vec.push_back(fm1);
	TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, labeler.setUpHook(fm_vec),"We currently support only 2- or 3-channel ICPL")
	fm_vec.push_back(fm2);
	labeler.setUpHook(fm_vec);
	fm_vec.push_back(fm3);
	labeler.setUpHook(fm_vec);
	fm_vec.push_back(fm4);
	TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, labeler.setUpHook(fm_vec),"We currently support only 2- or 3-channel ICPL")
}
END_SECTION

START_SECTION((void postDigestHook(FeatureMapSimVector &)))
{
	std::cerr <<"*********************************postDigestHook****************************************";
	FeatureMapSimVector feature_maps;


//*******************************Peptide-Labeling*******************************************
	Param p1;
  p1.setValue("label_proteins","true");		//set to true for protein_labeling


	//************2 Channel Protein Labeler********************
	createTestFeatureMapSimVector_(feature_maps, false);

	ICPLLabeler protein_labeler;
	protein_labeler.setParameters(p1);

	protein_labeler.setUpHook(feature_maps);		// Labeling
	digestFeaturesMapSimVector_(feature_maps);		// Digest
	// maps are digested by now
	protein_labeler.postDigestHook(feature_maps);	// Merge

	//Überprüft, ob das Ergebnis des Labelings und des Verdaus das erwartete Ergebnis ist.
	TEST_EQUAL(feature_maps.size(), 1)
	ABORT_IF(feature_maps.size() != 1)

	TEST_EQUAL(feature_maps[0].size(), 12)
	ABORT_IF(feature_maps[0].size() != 12)

	TEST_EQUAL(feature_maps[0][0].getIntensity(), 50)
	TEST_EQUAL(feature_maps[0][0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))AAAAAAAK")

	TEST_EQUAL(feature_maps[0][1].getIntensity(), 100)
	TEST_EQUAL(feature_maps[0][1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))CNAR")

	TEST_EQUAL(feature_maps[0][2].getIntensity(), 120)
	TEST_EQUAL(feature_maps[0][2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))LDR")

	TEST_EQUAL(feature_maps[0][3].getIntensity(), 110)
	TEST_EQUAL(feature_maps[0][3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))VNAAAAAAR")

	TEST_EQUAL(feature_maps[0][4].getIntensity(), 250)
	TEST_EQUAL(feature_maps[0][4].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAA")

	TEST_EQUAL(feature_maps[0][5].getIntensity(), 120)
	TEST_EQUAL(feature_maps[0][5].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CEL")

	TEST_EQUAL(feature_maps[0][6].getIntensity(), 180)
	TEST_EQUAL(feature_maps[0][6].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNCNCN")

	TEST_EQUAL(feature_maps[0][7].getIntensity(), 225)
	TEST_EQUAL(feature_maps[0][7].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "VNCNCNAAAA")

	TEST_EQUAL(feature_maps[0][8].getIntensity(), 200)
	TEST_EQUAL(feature_maps[0][8].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)AAAAAAAK")

	TEST_EQUAL(feature_maps[0][9].getIntensity(), 80)
	TEST_EQUAL(feature_maps[0][9].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)CNAR")

	TEST_EQUAL(feature_maps[0][10].getIntensity(), 100)
	TEST_EQUAL(feature_maps[0][10].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)CNHAADDAAAAA")

	TEST_EQUAL(feature_maps[0][11].getIntensity(), 115)
	TEST_EQUAL(feature_maps[0][11].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)VNAAAAAAR")

	
	//************3 Channel Protein Labeler********************
	createTestFeatureMapSimVector_(feature_maps, true);

	ICPLLabeler three_channel_protein_labeler;
	three_channel_protein_labeler.setParameters(p1);

	three_channel_protein_labeler.setUpHook(feature_maps);		//.Labeling
	digestFeaturesMapSimVector_(feature_maps);					//.Digest
	// maps are digested by now
	three_channel_protein_labeler.postDigestHook(feature_maps);	//.Merge

	//Überprüft, ob das Ergebnis des Labelings und des Verdaus das erwartete Ergebnis ist.
	TEST_EQUAL(feature_maps.size(), 1)
	ABORT_IF(feature_maps.size() != 1)

	TEST_EQUAL(feature_maps[0].size(), 16)
	ABORT_IF(feature_maps[0].size() != 16)

	TEST_EQUAL(feature_maps[0][0].getIntensity(), 30)
	TEST_EQUAL(feature_maps[0][0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:13C(6))AAAAAAAK")

  TEST_EQUAL(feature_maps[0][7].getIntensity(), 50)
  TEST_EQUAL(feature_maps[0][7].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))AAAAAAAK")

  TEST_EQUAL(feature_maps[0][1].getIntensity(), 130)
  TEST_EQUAL(feature_maps[0][1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:13C(6))CNAR")

  TEST_EQUAL(feature_maps[0][8].getIntensity(), 100)
  TEST_EQUAL(feature_maps[0][8].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))CNAR")

  TEST_EQUAL(feature_maps[0][2].getIntensity(), 70)
  TEST_EQUAL(feature_maps[0][2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:13C(6))LDR")

  TEST_EQUAL(feature_maps[0][9].getIntensity(), 120)
  TEST_EQUAL(feature_maps[0][9].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))LDR")

  TEST_EQUAL(feature_maps[0][6].getIntensity(), 310)
  TEST_EQUAL(feature_maps[0][6].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CNCNCN")

  TEST_EQUAL(feature_maps[0][4].getIntensity(), 280)
  TEST_EQUAL(feature_maps[0][4].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAA")

  TEST_EQUAL(feature_maps[0][5].getIntensity(), 190)
  TEST_EQUAL(feature_maps[0][5].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "CEL")

  TEST_EQUAL(feature_maps[0][3].getIntensity(), 80)
  TEST_EQUAL(feature_maps[0][3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:13C(6))YCYCY")

	TEST_EQUAL(feature_maps[0][10].getIntensity(), 110)
	TEST_EQUAL(feature_maps[0][10].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))VNAAAAAAR")

	TEST_EQUAL(feature_maps[0][11].getIntensity(), 225)
	TEST_EQUAL(feature_maps[0][11].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "VNCNCNAAAA")

	TEST_EQUAL(feature_maps[0][12].getIntensity(), 200)
	TEST_EQUAL(feature_maps[0][12].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)AAAAAAAK")

	TEST_EQUAL(feature_maps[0][13].getIntensity(), 80)
	TEST_EQUAL(feature_maps[0][13].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)CNAR")

	TEST_EQUAL(feature_maps[0][14].getIntensity(), 100)
	TEST_EQUAL(feature_maps[0][14].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)CNHAADDAAAAA")

	TEST_EQUAL(feature_maps[0][15].getIntensity(), 115)
	TEST_EQUAL(feature_maps[0][15].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)VNAAAAAAR")


//*******************************Peptide-Labeling*******************************************

	Param p2;
  p2.setValue("label_proteins","false");		//set to false for peptide_labeling


	//************2 Channel Peptide Labeler
	createTestFeatureMapSimVector_(feature_maps, false);

	ICPLLabeler peptide_labeler;
	peptide_labeler.setParameters(p2);
  peptide_labeler.setUpHook(feature_maps);
	digestFeaturesMapSimVector_(feature_maps);						// Digest
	peptide_labeler.postDigestHook(feature_maps);					// Labeling & Merge

	//Überprüft, ob das Ergebnis des Labelings und des Verdaus das erwartete Ergebnis ist.
	TEST_EQUAL(feature_maps.size(), 1)
	ABORT_IF(feature_maps.size() != 1)

	TEST_EQUAL(feature_maps[0].size(), 15)
	ABORT_IF(feature_maps[0].size() != 15)

	TEST_EQUAL(feature_maps[0][0].getIntensity(), 50)
	TEST_EQUAL(feature_maps[0][0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))AAAAA")

  TEST_EQUAL(feature_maps[0][1].getIntensity(), 50)
  TEST_EQUAL(feature_maps[0][1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))AAAAAAAK")

  TEST_EQUAL(feature_maps[0][2].getIntensity(), 120)
  TEST_EQUAL(feature_maps[0][2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))CEL")

	TEST_EQUAL(feature_maps[0][3].getIntensity(), 100)
	TEST_EQUAL(feature_maps[0][3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))CNAR")

	TEST_EQUAL(feature_maps[0][4].getIntensity(), 100)
	TEST_EQUAL(feature_maps[0][4].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))CNCNCN")

	TEST_EQUAL(feature_maps[0][5].getIntensity(), 120)
	TEST_EQUAL(feature_maps[0][5].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))LDR")

	TEST_EQUAL(feature_maps[0][6].getIntensity(), 110)
	TEST_EQUAL(feature_maps[0][6].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))VNAAAAAAR")

	TEST_EQUAL(feature_maps[0][7].getIntensity(), 110)
	TEST_EQUAL(feature_maps[0][7].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))VNCNCNAAAA")

	TEST_EQUAL(feature_maps[0][8].getIntensity(), 200)
	TEST_EQUAL(feature_maps[0][8].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)AAAAA")

	TEST_EQUAL(feature_maps[0][9].getIntensity(), 200)
	TEST_EQUAL(feature_maps[0][9].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)AAAAAAAK")

	TEST_EQUAL(feature_maps[0][10].getIntensity(), 80)
	TEST_EQUAL(feature_maps[0][10].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)CNAR")

	TEST_EQUAL(feature_maps[0][11].getIntensity(), 80)
	TEST_EQUAL(feature_maps[0][11].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)CNCNCN")

	TEST_EQUAL(feature_maps[0][12].getIntensity(), 100)
	TEST_EQUAL(feature_maps[0][12].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)CNHAADDAAAAA")

	TEST_EQUAL(feature_maps[0][13].getIntensity(), 115)
	TEST_EQUAL(feature_maps[0][13].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)VNAAAAAAR")

	TEST_EQUAL(feature_maps[0][14].getIntensity(), 115)
	TEST_EQUAL(feature_maps[0][14].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)VNCNCNAAAA")


	//************3 Channel Peptide Labeler
	createTestFeatureMapSimVector_(feature_maps, true);

	ICPLLabeler three_channel_peptide_labeler;
	three_channel_peptide_labeler.setParameters(p2);

	digestFeaturesMapSimVector_(feature_maps);						// Digest
	three_channel_peptide_labeler.postDigestHook(feature_maps);		// Labeling & Merge

	//Überprüft, ob das Ergebnis des Labelings und des Verdaus das erwartete Ergebnis ist.
	TEST_EQUAL(feature_maps.size(), 1)
	ABORT_IF(feature_maps.size() != 1)

	TEST_EQUAL(feature_maps[0].size(), 22)
	ABORT_IF(feature_maps[0].size() != 22)

	TEST_EQUAL(feature_maps[0][0].getIntensity(), 30)
	TEST_EQUAL(feature_maps[0][0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:13C(6))AAAAA")

  TEST_EQUAL(feature_maps[0][7].getIntensity(), 50)
  TEST_EQUAL(feature_maps[0][7].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))AAAAA")

  TEST_EQUAL(feature_maps[0][1].getIntensity(), 30)
  TEST_EQUAL(feature_maps[0][1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:13C(6))AAAAAAAK")

  TEST_EQUAL(feature_maps[0][8].getIntensity(), 50)
  TEST_EQUAL(feature_maps[0][8].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))AAAAAAAK")

  TEST_EQUAL(feature_maps[0][2].getIntensity(), 70)
  TEST_EQUAL(feature_maps[0][2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:13C(6))CEL")

  TEST_EQUAL(feature_maps[0][9].getIntensity(), 120)
  TEST_EQUAL(feature_maps[0][9].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))CEL")

  TEST_EQUAL(feature_maps[0][3].getIntensity(), 130)
  TEST_EQUAL(feature_maps[0][3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:13C(6))CNAR")

  TEST_EQUAL(feature_maps[0][10].getIntensity(), 100)
  TEST_EQUAL(feature_maps[0][10].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))CNAR")

  TEST_EQUAL(feature_maps[0][4].getIntensity(), 130)
  TEST_EQUAL(feature_maps[0][4].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:13C(6))CNCNCN")

  TEST_EQUAL(feature_maps[0][11].getIntensity(), 100)
  TEST_EQUAL(feature_maps[0][11].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))CNCNCN")

  TEST_EQUAL(feature_maps[0][5].getIntensity(), 70)
  TEST_EQUAL(feature_maps[0][5].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:13C(6))LDR")

  TEST_EQUAL(feature_maps[0][12].getIntensity(), 120)
  TEST_EQUAL(feature_maps[0][12].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))LDR")

  TEST_EQUAL(feature_maps[0][6].getIntensity(), 80)
  TEST_EQUAL(feature_maps[0][6].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:13C(6))YCYCY")

	TEST_EQUAL(feature_maps[0][13].getIntensity(), 110)
	TEST_EQUAL(feature_maps[0][13].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))VNAAAAAAR")

	TEST_EQUAL(feature_maps[0][14].getIntensity(), 110)
	TEST_EQUAL(feature_maps[0][14].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL:2H(4))VNCNCNAAAA")

	TEST_EQUAL(feature_maps[0][15].getIntensity(), 200)
	TEST_EQUAL(feature_maps[0][15].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)AAAAA")

	TEST_EQUAL(feature_maps[0][16].getIntensity(), 200)
	TEST_EQUAL(feature_maps[0][16].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)AAAAAAAK")

	TEST_EQUAL(feature_maps[0][17].getIntensity(), 80)
	TEST_EQUAL(feature_maps[0][17].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)CNAR")

	TEST_EQUAL(feature_maps[0][18].getIntensity(), 80)
	TEST_EQUAL(feature_maps[0][18].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)CNCNCN")

	TEST_EQUAL(feature_maps[0][19].getIntensity(), 100)
	TEST_EQUAL(feature_maps[0][19].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)CNHAADDAAAAA")

	TEST_EQUAL(feature_maps[0][20].getIntensity(), 115)
	TEST_EQUAL(feature_maps[0][20].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)VNAAAAAAR")

	TEST_EQUAL(feature_maps[0][21].getIntensity(), 115)
	TEST_EQUAL(feature_maps[0][21].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "(ICPL)VNCNCNAAAA")

}
END_SECTION

START_SECTION((void postRawMSHook(FeatureMapSimVector &)))
{
	// TODO
}
END_SECTION

// just to call the methods once
ICPLLabeler dummyLabeler;
FeatureMapSimVector empty;

START_SECTION((void preCheck(Param &param) const ))
{
	Param p;
	dummyLabeler.preCheck(p);

	// preCheck has no content
	NOT_TESTABLE
}
END_SECTION


START_SECTION((void postRTHook(FeatureMapSimVector &)))
{
	// we do not modify the map in this step
	dummyLabeler.postRTHook(empty);
	NOT_TESTABLE
}
END_SECTION

START_SECTION((void postDetectabilityHook(FeatureMapSimVector &)))
{
	// we do not modify the map in this step
	dummyLabeler.postDetectabilityHook(empty);
	NOT_TESTABLE
}
END_SECTION

START_SECTION((void postIonizationHook(FeatureMapSimVector &)))
{
	// we do not modify the map in this step
	dummyLabeler.postIonizationHook(empty);
	NOT_TESTABLE
}
END_SECTION

MSSimExperiment exp;
START_SECTION((void postRawTandemMSHook(FeatureMapSimVector &, MSSimExperiment &)))
{
	// we do not modify the map in this step
	dummyLabeler.postRawTandemMSHook(empty,exp);
	NOT_TESTABLE
}
END_SECTION

START_SECTION((static BaseLabeler* create()))
{
	BaseLabeler* labeler = ICPLLabeler::create();
	TEST_NOT_EQUAL(labeler, 0)
	delete labeler;
}
END_SECTION

START_SECTION((static const String getProductName()))
{
	TEST_EQUAL(ICPLLabeler::getProductName(), "ICPL")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
