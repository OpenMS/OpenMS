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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/AbsoluteQuantitationStandards.h>
#include <OpenMS/FORMAT/AbsoluteQuantitationStandardsFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AbsoluteQuantitationStandards, "$Id$")

/////////////////////////////////////////////////////////////

AbsoluteQuantitationStandards* ptr = 0;
AbsoluteQuantitationStandards* null_ptr = 0;
const String calib_conc_path = OPENMS_GET_TEST_DATA_PATH("AbsoluteQuantitationStandardsFile_150516_calibration_concentrations_cut.csv");
const String features_path = OPENMS_GET_TEST_DATA_PATH("170808_Jonathan_yeast_Sacc1_1x.featureXML");

START_SECTION(AbsoluteQuantitationStandards())
{
  ptr = new AbsoluteQuantitationStandards();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~AbsoluteQuantitationStandards())
{
  delete ptr;
}
END_SECTION

ptr = new AbsoluteQuantitationStandards();

START_SECTION(void mapComponentsToConcentrations(
  const std::vector<runConcentration>& run_concentrations,
  const std::vector<FeatureMap>& feature_maps,
  std::map<String, std::vector<featureConcentration>>& components_to_concentrations
) const)
{
  std::vector<AbsoluteQuantitationStandards::runConcentration> runs;
  AbsoluteQuantitationStandardsFile aqsf;
  aqsf.load(calib_conc_path, runs);
  FeatureMap fmap;
  FeatureXMLFile fxmlf;
  fxmlf.load(features_path, fmap);
  fmap.setMetaValue("sample_name", "150516_CM1_Level1");
  std::vector<FeatureMap> feature_maps;
  feature_maps.push_back(fmap);
  TEST_EQUAL(feature_maps.size(), 1)
  std::map<String, AbsoluteQuantitationStandards::featureConcentration> m;
  ptr->mapComponentsToConcentrations(runs, feature_maps, m);
  TEST_NOT_EQUAL(m.size(), 0)
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
