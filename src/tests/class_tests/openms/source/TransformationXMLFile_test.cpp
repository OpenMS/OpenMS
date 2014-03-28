// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/TransformationXMLFile.h>

///////////////////////////

START_TEST(FASTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TransformationXMLFile* ptr = 0;
TransformationXMLFile* nullPointer = 0;

START_SECTION((TransformationXMLFile()))
{
	ptr = new TransformationXMLFile();
  TEST_NOT_EQUAL(ptr,nullPointer);
}
END_SECTION

START_SECTION([EXTRA] static bool isValid(const String& filename))
{
	TransformationXMLFile f;
	TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_1.trafoXML"), std::cerr),true);
	TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_2.trafoXML"), std::cerr),true);
//	TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_3.trafoXML")),false);
//	TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_4.trafoXML")),true);
}
END_SECTION


START_SECTION(void load(const String& filename, TransformationDescription& transformation))
{
	TransformationDescription trafo;
	TransformationXMLFile trafo_xml;
	Param params;

	trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_1.trafoXML"), trafo);
	TEST_STRING_EQUAL(trafo.getModelType(), "none");
	params = trafo.getModelParameters();
	TEST_EQUAL(params.empty(), true);

	trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_2.trafoXML"),trafo);
	TEST_STRING_EQUAL(trafo.getModelType(), "linear");
	params = trafo.getModelParameters();
	TEST_EQUAL(params.size(), 2);
	TEST_REAL_SIMILAR(params.getValue("slope"), 3.141592653589793238);
	TEST_REAL_SIMILAR(params.getValue("intercept"), 2.718281828459045235);

//	trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_4.trafoXML"), trafo);
//	TEST_STRING_EQUAL(trafo.getModelType(), "interpolated");
//	trafo.getModelParameters(params);
//	TEST_EQUAL(params.getValue("interpolation_type"), "linear");
//	TEST_EQUAL(trafo.getDataPoints().size(), 3);
//	TEST_REAL_SIMILAR(trafo.getDataPoints()[0].first, 1.2);
//	TEST_REAL_SIMILAR(trafo.getDataPoints()[1].first, 2.2);
//	TEST_REAL_SIMILAR(trafo.getDataPoints()[2].first, 3.2);
//	TEST_REAL_SIMILAR(trafo.getDataPoints()[0].second, 5.2);
//	TEST_REAL_SIMILAR(trafo.getDataPoints()[1].second, 6.25);
//	TEST_REAL_SIMILAR(trafo.getDataPoints()[2].second, 7.3);
}
END_SECTION

START_SECTION(void store(String filename, const TransformationDescription& transformation))
{
	TransformationDescription trafo, trafo2;
	TransformationXMLFile trafo_xml;

	String tmp_file_none;
	Param params;
	trafo.fitModel("none", params);
	NEW_TMP_FILE(tmp_file_none);
	trafo_xml.store(tmp_file_none, trafo);
	trafo_xml.load(tmp_file_none, trafo2);
	TEST_STRING_EQUAL(trafo2.getModelType(), "none");
	params = trafo2.getModelParameters();
	TEST_EQUAL(params.empty(), true);
	{
	  DoubleReal pre_image = 234255132.43212;
	  DoubleReal image = trafo.apply(pre_image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
	}
	String tmp_file_linear;
	NEW_TMP_FILE(tmp_file_linear);
	params.setValue("slope", 3.141592653589793238);
	params.setValue("intercept", 2.718281828459045235);
	trafo.fitModel("linear", params);
	trafo_xml.store(tmp_file_linear, trafo);
	trafo_xml.load(tmp_file_linear, trafo2);
	TEST_STRING_EQUAL(trafo.getModelType(), "linear");
	params.clear();
	params = trafo2.getModelParameters();
	TEST_EQUAL(params.size(), 2);
	TEST_REAL_SIMILAR(params.getValue("slope"), 3.141592653589793238);
	TEST_REAL_SIMILAR(params.getValue("intercept"), 2.718281828459045235);

  TEST_EXCEPTION(Exception::IllegalArgument, trafo.fitModel("mumble_pfrwoarpfz"));

	String tmp_file_pairs;
	NEW_TMP_FILE(tmp_file_pairs);
	TransformationDescription::DataPoints pairs;
  pairs.push_back(make_pair(1.2f,5.2f));
  pairs.push_back(make_pair(3.2f,7.3f));
  pairs.push_back(make_pair(2.2f,6.25f));
  pairs.push_back(make_pair(2.2f,3.1f));
  pairs.push_back(make_pair(2.2f,7.25f));
  pairs.push_back(make_pair(3.0f,8.5f));
  pairs.push_back(make_pair(3.1f,4.7f));
  pairs.push_back(make_pair(1.7f,6.0f));
  pairs.push_back(make_pair(2.9f,4.7f));
  pairs.push_back(make_pair(4.2f,5.0f));
  pairs.push_back(make_pair(3.7f,-2.4f));
  trafo.setDataPoints(pairs);
	params.clear();
	trafo.fitModel("interpolated", params);
  trafo_xml.store(tmp_file_pairs,trafo);
  trafo_xml.load(tmp_file_pairs,trafo2);
  TEST_STRING_EQUAL(trafo2.getModelType(), "interpolated");
	params.clear();
	params = trafo2.getModelParameters();
  TEST_EQUAL(trafo2.getDataPoints().size(), 11);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[0].first, 1.2);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[0].second, 5.2);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[10].first, 3.7);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[10].second, -2.4);

}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
