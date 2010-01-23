// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <iostream>

using namespace std;
using namespace OpenMS;


/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithmIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithmIdentification* ptr = 0;
START_SECTION((MapAlignmentAlgorithmIdentification()))
	ptr = new MapAlignmentAlgorithmIdentification();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

	
START_SECTION((virtual ~MapAlignmentAlgorithmIdentification()))
	delete ptr;
END_SECTION


START_SECTION((static MapAlignmentAlgorithm* create()))
	TEST_NOT_EQUAL(MapAlignmentAlgorithmIdentification::create(), 0)
END_SECTION


START_SECTION((static String getProductName()))
	TEST_EQUAL(MapAlignmentAlgorithmIdentification::getProductName(), "identification")
END_SECTION


START_SECTION((virtual void alignPeptideIdentifications(std::vector<std::vector<PeptideIdentification> >&, std::vector<TransformationDescription>&)))
{
	vector<vector<PeptideIdentification> > peptides(2);
	vector<ProteinIdentification> proteins;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmIdentification_test_1.idXML"),	proteins, peptides[0]);
 	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmIdentification_test_2.idXML"),	proteins, peptides[1]);
// 	cout << "Input:\n";
// 	for (Size i = 0; i < peptides[0].size(); ++i)
// 	{
// 		cout << left << setw(20)
// 				 << peptides[0][i].getHits()[0].getSequence().toString() + ":"
// 				 << peptides[0][i].getMetaValue("RT") << "\t"
// 				 << peptides[1][i].getMetaValue("RT") << endl;
// 	}
	vector<TransformationDescription> transforms(2);
	MapAlignmentAlgorithm* aligner = Factory<MapAlignmentAlgorithm>::create(
		"identification");

	Param params = aligner->getParameters();
	params.setValue("peptide_score_threshold", 0.0);
	params.setValue("num_breakpoints", 10);
	aligner->setParameters(params);
	aligner->setLogType(ProgressLogger::CMD);
	aligner->alignPeptideIdentifications(peptides, transforms);
// 	cout << "Output (transformed):\n";
// 	for (Size i = 0; i < peptides[0].size(); ++i)
// 	{
// 		cout << left << setw(20)
// 				 << peptides[0][i].getHits()[0].getSequence().toString() + ":"
// 				 << peptides[0][i].getMetaValue("RT") << "\t"
// 				 << peptides[1][i].getMetaValue("RT") << endl;
// 	}
	for (Size i = 0; i < peptides[0].size(); ++i)
	{
		TEST_REAL_SIMILAR(peptides[0][i].getMetaValue("RT"),
											peptides[1][i].getMetaValue("RT"));
	}
	
}
END_SECTION


START_SECTION((virtual void alignPeakMaps(std::vector<MSExperiment<> >&, std::vector<TransformationDescription>&)))
{
	// largely the same as "alignPeptideIdentifications"
  NOT_TESTABLE;
}
END_SECTION


START_SECTION((virtual void alignFeatureMaps(std::vector<FeatureMap<> >&, std::vector<TransformationDescription>&)))
{
	// largely the same as "alignPeptideIdentifications"
  NOT_TESTABLE;
}
END_SECTION


// can't test protected methods...

// START_SECTION((DoubleReal median_(DoubleList&, bool)))
// {
// 	DoubleList values;
// 	TEST_EXCEPTION(Exception::IllegalArgument, median_(values)); // empty list
// 	// -1.0, -0.5, ..., 2.0 scrambled:
// 	values << 0.5 << -1.0 << 0.0 << 1.5 << 1.0 << -0.5 << 2.0;
// 	TEST_EQUAL(median_(values, false), 0.5);
// 	TEST_EQUAL(median_(values, true), 0.5); // should be sorted now
// 	values << 2.5; // even number of values
// 	TEST_EQUAL(median_(values, true), 0.75);
// }
// END_SECTION


// START_SECTION((void computeMedians_(SeqToList&, SeqToValue&, bool)))
// {
// 	map<String, DoubleList> seq_to_list;
// 	map<String, DoubleReal> seq_to_value;
// 	seq_to_list["ABC"] << -1.0 << 2.5 << 0.5 << -3.5;
// 	seq_to_list["DEF"] << 1.5 << -2.5 << -1;
// 	computeMedians_(seq_to_list, seq_to_value, false);
// 	TEST_EQUAL(seq_to_value.size(), 2);
// 	TEST_EQUAL(seq_to_value["ABC"], -0.25);
// 	TEST_EQUAL(seq_to_value["DEF"], -1);
// 	seq_to_value.clear();
// 	computeMedians_(seq_to_list, seq_to_value, true); // should be sorted now
// 	TEST_EQUAL(seq_to_value.size(), 2);
// 	TEST_EQUAL(seq_to_value["ABC"], -0.25);
// 	TEST_EQUAL(seq_to_value["DEF"], -1);
// }
// END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
