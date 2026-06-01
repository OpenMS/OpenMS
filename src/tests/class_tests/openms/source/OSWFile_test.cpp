// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/OSWFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;


/*
Creating an OWS test database from pyProphet osw-db files (which provided the SCORE_MS2 table, which is missing in OpenSwathWorkflow outputs)

Using DBBrowser for SQLite:
1) open a full blown OSW file
2) export the DB schema to a new test DB file (without the data)
3) To fill the empty test DB, we start at the highest (=Protein) level and pick all dependent rows from downstream tables using the SQL-commands below.
The resulting data can be copied "as SQL-commands" and inserted into the test DB. This has to be done table by table.
Sometimes, only the first columns must be selected using the mouse, because the query contains more columns from joins with other tables.
Note: not all tables are populated. Only the ones we need at the moment..

SELECT * FROM PROTEIN WHERE ID IN (2,3866)
--> insert data into PROTEIN table

select * from PEPTIDE
INNER JOIN (SELECT PEPTIDE_ID as PID, PROTEIN_ID as POD FROM PEPTIDE_PROTEIN_MAPPING) AS PEP ON PEP.PID=PEPTIDE.ID WHERE PEP.POD IN (2,3866)
--> insert data into PEPTIDE table

SELECT * FROM PEPTIDE_PROTEIN_MAPPING WHERE PROTEIN_ID IN (2,3866)
--> insert data into PEPTIDE_PROTEIN_MAPPING table


select PEPTIDE_ID, PRECURSOR_ID FROM PRECURSOR_PEPTIDE_MAPPING
INNER JOIN (SELECT PEPTIDE_ID as PID, PROTEIN_ID as POD FROM PEPTIDE_PROTEIN_MAPPING) AS PEP ON PEP.PID=PEPTIDE_ID WHERE PEP.POD IN (2,3866) AND PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID=PEP.PID
--> insert data into PRECURSOR_PEPTIDE_MAPPING table

SELECT * FROM PRECURSOR
INNER JOIN (select PEPTIDE_ID, PRECURSOR_ID FROM PRECURSOR_PEPTIDE_MAPPING
INNER JOIN (SELECT PEPTIDE_ID as PID, PROTEIN_ID as POD FROM PEPTIDE_PROTEIN_MAPPING) AS PEP ON PEP.PID=PEPTIDE_ID WHERE PEP.POD IN (2,3866) AND PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID=PEP.PID) as PC ON PC.PRECURSOR_ID=PRECURSOR.ID
--> insert data into PRECURSOR table

SELECT * from FEATURE
INNER JOIN (SELECT * FROM PRECURSOR
INNER JOIN (select PEPTIDE_ID, PRECURSOR_ID FROM PRECURSOR_PEPTIDE_MAPPING
INNER JOIN (SELECT PEPTIDE_ID as PID, PROTEIN_ID as POD FROM PEPTIDE_PROTEIN_MAPPING) AS PEP ON PEP.PID=PEPTIDE_ID WHERE PEP.POD IN (2,3866) AND PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID=PEP.PID) as PC ON PC.PRECURSOR_ID=PRECURSOR.ID) AS PREC ON PREC.ID=FEATURE.PRECURSOR_ID
--> insert data into FEATURE table

SELECT * from FEATURE_TRANSITION
INNER JOIN (SELECT * from FEATURE
INNER JOIN (SELECT * FROM PRECURSOR
INNER JOIN (select PEPTIDE_ID, PRECURSOR_ID FROM PRECURSOR_PEPTIDE_MAPPING
INNER JOIN (SELECT PEPTIDE_ID as PID, PROTEIN_ID as POD FROM PEPTIDE_PROTEIN_MAPPING) AS PEP ON PEP.PID=PEPTIDE_ID WHERE PEP.POD IN (2,3866) AND PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID=PEP.PID) as PC ON PC.PRECURSOR_ID=PRECURSOR.ID) AS PREC ON PREC.ID=FEATURE.PRECURSOR_ID) AS FEAT ON FEAT.ID=FEATURE_TRANSITION.FEATURE_ID
--> insert data into FEATURE_TRANSITION table

SELECT *  from TRANSITION
INNER JOIN (SELECT DISTINCT TRANSITION_ID from FEATURE_TRANSITION
INNER JOIN (SELECT * from FEATURE
INNER JOIN (SELECT * FROM PRECURSOR
INNER JOIN (select PEPTIDE_ID, PRECURSOR_ID FROM PRECURSOR_PEPTIDE_MAPPING
INNER JOIN (SELECT PEPTIDE_ID as PID, PROTEIN_ID as POD FROM PEPTIDE_PROTEIN_MAPPING) AS PEP ON PEP.PID=PEPTIDE_ID WHERE PEP.POD IN (2,3866) AND PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID=PEP.PID) as PC ON PC.PRECURSOR_ID=PRECURSOR.ID) AS PREC ON PREC.ID=FEATURE.PRECURSOR_ID) AS FEAT ON FEAT.ID=FEATURE_TRANSITION.FEATURE_ID) AS FEATTR ON FEATTR.TRANSITION_ID=TRANSITION.ID
--> insert data into TRANSITION table

SELECT * from SCORE_MS2
INNER JOIN (SELECT * from FEATURE
INNER JOIN (SELECT * FROM PRECURSOR
INNER JOIN (select PEPTIDE_ID, PRECURSOR_ID FROM PRECURSOR_PEPTIDE_MAPPING
INNER JOIN (SELECT PEPTIDE_ID as PID, PROTEIN_ID as POD FROM PEPTIDE_PROTEIN_MAPPING) AS PEP ON PEP.PID=PEPTIDE_ID WHERE PEP.POD IN (2,3866) AND PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID=PEP.PID) as PC ON PC.PRECURSOR_ID=PRECURSOR.ID) AS PREC ON PREC.ID=FEATURE.PRECURSOR_ID)
		AS FEAT ON FEAT.ID=SCORE_MS2.FEATURE_ID
--> insert data into SCORE_MS2 table  (only available after pyProphet ran!)


*/

void checkData(OSWData& res)
{
	const OSWProtein& prot = *res.getProteins().begin();
	TEST_EQUAL(prot.getAccession(), "1/P00167ups|CYB5_HUMAN_UPS");
	const OSWPeptidePrecursor& prec = *prot.getPeptidePrecursors().begin();
	TEST_EQUAL(prec.getCharge(), 2);
	TEST_EQUAL(prec.isDecoy(), false);
	TEST_REAL_SIMILAR(prec.getPCMz(), 1103.4676);
	TEST_EQUAL(prec.getSequence(), "EQAGGDATENFEDVGHSTDAR");
	TEST_EQUAL(prec.getFeatures().size(), 5);
	const std::vector<UInt32> tr{ 236830, 236831, 236832, 236833, 236834 };
	const auto& trd = prec.getFeatures().back().getTransitionIDs();
	TEST_TRUE(trd == tr);
	// check last transition
	const OSWProtein& prot_last = res.getProteins().back();
	TEST_EQUAL(prot_last.getPeptidePrecursors().back().getFeatures().back().getTransitionIDs().back(), 99);

	// all features should have 5 transitions
	for (const auto& prot : res.getProteins())
	{
		for (const auto& pc : prot.getPeptidePrecursors())
		{
			for (const auto& feat : pc.getFeatures())
			{
				TEST_EQUAL(feat.getTransitionIDs().size(), 5);
			}
		}
	}
}

START_TEST(OSWFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(UInt64 getRunID() const)
	OSWFile oswf(OPENMS_GET_TEST_DATA_PATH("OSWFile_test.osw"));
	TEST_EQUAL(oswf.getRunID(), 6996169951924032342);
END_SECTION

START_SECTION(void read(OSWData& swath_result))
	OSWData res;
	OSWFile oswf(OPENMS_GET_TEST_DATA_PATH("OSWFile_test.osw"));
	oswf.read(res);
	TEST_EQUAL(res.getProteins().size(), 2);
	TEST_EQUAL(res.transitionCount(), 140);
	TEST_EQUAL(res.getRunID(), 6996169951924032342);
	checkData(res);
END_SECTION			

START_SECTION(void readMinimal(OSWData & swath_result))
	OSWData res;
	OSWFile oswf(OPENMS_GET_TEST_DATA_PATH("OSWFile_test.osw"));
	oswf.readMinimal(res);
	TEST_EQUAL(res.getProteins().size(), 2);
	TEST_EQUAL(res.transitionCount(), 140);

	TEST_EQUAL(res.getRunID(), 6996169951924032342);

	// make sure proteins are actually empty
	TEST_EQUAL(res.getProteins()[0].getPeptidePrecursors().empty(), true);
	TEST_EQUAL(res.getProteins()[1].getPeptidePrecursors().empty(), true);

	// now fill them...
	for (Size i = 0; i < res.getProteins().size(); ++i)
	{
		oswf.readProtein(res, i);
	}
	checkData(res);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



