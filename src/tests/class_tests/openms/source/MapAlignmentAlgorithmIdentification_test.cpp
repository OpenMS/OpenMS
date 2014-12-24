// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <iostream>

#include <OpenMS/CONCEPT/Factory.h>

using namespace std;
using namespace OpenMS;

/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithmIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithmIdentification* ptr = 0;
MapAlignmentAlgorithmIdentification* nullPointer = 0;
MapAlignmentAlgorithm* base_nullPointer = 0;
START_SECTION((MapAlignmentAlgorithmIdentification()))
	ptr = new MapAlignmentAlgorithmIdentification();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION


START_SECTION((virtual ~MapAlignmentAlgorithmIdentification()))
	delete ptr;
END_SECTION


START_SECTION((static MapAlignmentAlgorithm* create()))
  TEST_NOT_EQUAL(MapAlignmentAlgorithmIdentification::create(), base_nullPointer)
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
	vector<TransformationDescription> transforms(2);
	MapAlignmentAlgorithm* aligner = Factory<MapAlignmentAlgorithm>::create(
		"identification");

	Param params = aligner->getParameters();
	params.setValue("peptide_score_threshold", 0.0);
	aligner->setParameters(params);
	aligner->setLogType(ProgressLogger::CMD);
	aligner->alignPeptideIdentifications(peptides, transforms);
	params.clear();
	aligner->fitModel("interpolated", params, transforms);
  MapAlignmentTransformer::transformPeptideIdentifications(peptides, transforms);
	for (Size i = 0; i < peptides[0].size(); ++i)
	{
		TEST_REAL_SIMILAR(peptides[0][i].getRT(), peptides[1][i].getRT());
	}
}
END_SECTION


START_SECTION((virtual void alignPeakMaps(std::vector<MSExperiment<> >&, std::vector<TransformationDescription>&)))
{
	// largely the same as "alignPeptideIdentifications"
  NOT_TESTABLE;
}
END_SECTION


START_SECTION((virtual void alignFeatureMaps(std::vector<FeatureMap >&, std::vector<TransformationDescription>&)))
{
	// largely the same as "alignPeptideIdentifications"
  NOT_TESTABLE;
}
END_SECTION


START_SECTION((virtual void alignConsensusMaps(std::vector<ConsensusMap>&, std::vector<TransformationDescription>&)))
{
	// largely the same as "alignPeptideIdentifications"
  NOT_TESTABLE;
}
END_SECTION


START_SECTION((virtual void setReference(Size reference_index=0, const String& reference_file="")))
{
	MapAlignmentAlgorithm* aligner = Factory<MapAlignmentAlgorithm>::create(
		"identification");
	aligner->setReference(1); // nothing happens
	TEST_EXCEPTION(Exception::FileNotFound,
								 aligner->setReference(0, "not-a-real-file.idXML"));
}
END_SECTION


// can't test protected methods...

// START_SECTION((double median_(DoubleList&, bool)))
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
// 	map<String, double> seq_to_value;
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
