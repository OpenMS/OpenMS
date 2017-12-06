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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

#include <iostream>

using namespace std;
using namespace OpenMS;

/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithmIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithmIdentification* ptr = nullptr;
MapAlignmentAlgorithmIdentification* nullPointer = nullptr;
START_SECTION((MapAlignmentAlgorithmIdentification()))
	ptr = new MapAlignmentAlgorithmIdentification();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION


START_SECTION((virtual ~MapAlignmentAlgorithmIdentification()))
	delete ptr;
END_SECTION

vector<vector<PeptideIdentification> > peptides(2);
vector<ProteinIdentification> proteins;
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmIdentification_test_1.idXML"),	proteins, peptides[0]);
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmIdentification_test_2.idXML"),	proteins, peptides[1]);

MapAlignmentAlgorithmIdentification aligner;
aligner.setLogType(ProgressLogger::CMD);
Param params = aligner.getParameters();
params.setValue("peptide_score_threshold", 0.0);
aligner.setParameters(params);
vector<double> reference_rts; // needed later

START_SECTION((template <typename DataType> void align(std::vector<DataType>& data, std::vector<TransformationDescription>& transformations, Int reference_index = -1)))
{
  // alignment without reference:
	vector<TransformationDescription> transforms;
	aligner.align(peptides, transforms);

  TEST_EQUAL(transforms.size(), 2);
  TEST_EQUAL(transforms[0].getDataPoints().size(), 10);
  TEST_EQUAL(transforms[1].getDataPoints().size(), 10);

  reference_rts.reserve(10);
  for (Size i = 0; i < transforms[0].getDataPoints().size(); ++i)
  {
    // both RT transforms should map to a common RT scale:
    TEST_REAL_SIMILAR(transforms[0].getDataPoints()[i].second,
                      transforms[1].getDataPoints()[i].second);
    reference_rts.push_back(transforms[0].getDataPoints()[i].first);
  }

  // alignment with internal reference:
  transforms.clear();
  aligner.align(peptides, transforms, 0);

  TEST_EQUAL(transforms.size(), 2);
  TEST_EQUAL(transforms[0].getModelType(), "identity");
  TEST_EQUAL(transforms[1].getDataPoints().size(), 10);

  for (Size i = 0; i < transforms[1].getDataPoints().size(); ++i)
  {
    // RT transform should map to RT scale of the reference:
    TEST_REAL_SIMILAR(transforms[1].getDataPoints()[i].second,
                      reference_rts[i]);
  }

  // algorithm works the same way for other input data types -> no extra tests
}
END_SECTION


START_SECTION((template <typename DataType> void setReference(DataType& data)))
{
  // alignment with external reference:
  aligner.setReference(peptides[0]);
  peptides.erase(peptides.begin());

  vector<TransformationDescription> transforms;
  aligner.align(peptides, transforms);

  TEST_EQUAL(transforms.size(), 1);
  TEST_EQUAL(transforms[0].getDataPoints().size(), 10);

  for (Size i = 0; i < transforms[0].getDataPoints().size(); ++i)
  {
    // RT transform should map to RT scale of the reference:
    TEST_REAL_SIMILAR(transforms[0].getDataPoints()[i].second,
                      reference_rts[i]);
  }
}
END_SECTION


// can't test protected methods...

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
