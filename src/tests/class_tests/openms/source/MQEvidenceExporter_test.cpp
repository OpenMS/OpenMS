// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Chris Bielow$
// $Authors: Valentin Noske, Vincent Musch$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/QC/MQEvidenceExporter.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/test_config.h>

///////////////////////////
///////////////////////////

START_TEST(MQEvidence, "$ID$")

using namespace OpenMS;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

File::TempDir dir;
const String path = dir.getPath();

START_SECTION(MQEvidence())
{
  MQEvidence *ptr = nullptr;
  MQEvidence *null_ptr = nullptr;
  ptr = new MQEvidence(path);
  TEST_NOT_EQUAL(ptr, null_ptr);
  delete ptr;
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION((void exportFeatureMap(
                    const OpenMS::FeatureMap &feature_map,
                    const OpenMS::ConsensusMap &cmap)))
{
  {
    MQEvidence evd(path);
    ConsensusMap cmap_one;
    ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_1.consensusXML"), cmap_one);
    ConsensusMap cmap_two;
    ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_2.consensusXML"), cmap_two);
    FeatureMap fmap_one;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_1.featureXML"), fmap_one);
    evd.exportFeatureMap(fmap_one, cmap_two);
    FeatureMap fmap_two;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_2.featureXML"), fmap_two);
    evd.exportFeatureMap(fmap_two, cmap_two);
    FeatureMap fmap_three;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_3.featureXML"), fmap_three);
    evd.exportFeatureMap(fmap_three, cmap_two);
    FeatureMap fmap_four;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_4.featureXML"), fmap_four);
    evd.exportFeatureMap(fmap_four, cmap_one);
    FeatureMap fmap_five;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_5.featureXML"), fmap_five);
    evd.exportFeatureMap(fmap_five, cmap_one);
    FeatureMap fmap_six;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_6.featureXML"), fmap_six);
    evd.exportFeatureMap(fmap_six, cmap_one);
  }
  String filename = path + "/evidence.txt";
  TEST_FILE_SIMILAR(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("MQEvidence_result.txt"));
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST