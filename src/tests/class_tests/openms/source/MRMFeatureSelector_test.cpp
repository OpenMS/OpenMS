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
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MRMFeatureQCFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
///////////////////////////

#define TRANSITIONTSVREADER_TESTING 1 

using namespace OpenMS;
using namespace std;

START_TEST(MRMFeatureSelector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMFeatureSelector* ptr = 0;
MRMFeatureSelector* null_ptr = 0;

const String features_path = OPENMS_GET_TEST_DATA_PATH("MRMFeatureSelector_150601_0_BloodProject01_PLT_QC_Broth-1_1.featureXML");
const String target_list_path = OPENMS_GET_TEST_DATA_PATH("MRMFeatureSelector_BloodProject01_SWATH.csv");
// const String components_path = OPENMS_GET_TEST_DATA_PATH("test_pyTOPP_MRMFeatureQCComponents.csv");
// const String components_groups_path = OPENMS_GET_TEST_DATA_PATH("test_pyTOPP_MRMFeatureQCComponentGroups.csv");

FeatureMap feature_map;
FeatureXMLFile feature_file;
feature_file.load(features_path, feature_map);

TransitionTSVFile tsv_reader;
TargetedExperiment targeted_exp;
tsv_reader.convertTSVToTargetedExperiment(target_list_path.c_str(), FileTypes::CSV, targeted_exp);

// MRMFeatureQCFile mrmfqcfile;
// MRMFeatureQC mrmfqc;
// mrmfqcfile.load(components_path, mrmfqc, false); // components file
// mrmfqcfile.load(components_groups_path, mrmfqc, true); // component groups file

START_SECTION(MRMFeatureSelector())
{
  ptr = new MRMFeatureSelector();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MRMFeatureSelector())
{
  delete ptr;
}
END_SECTION

ptr = new MRMFeatureSelector();

START_SECTION(getParameters().getValue("nn_threshold"))
{
  TEST_EQUAL(ptr->getParameters().getValue("nn_threshold"), 4.0)
}
END_SECTION

START_SECTION(setNNThreshold())
{
  TEST_EQUAL(ptr->getNNThreshold(), 4.0)
  ptr->setNNThreshold(5.0);
  TEST_EQUAL(ptr->getNNThreshold(), 5.0)
}
END_SECTION

START_SECTION(getLocalityWeight())
{
  TEST_EQUAL(ptr->getLocalityWeight(), false)
  ptr->setLocalityWeight(true);
  TEST_EQUAL(ptr->getLocalityWeight(), true)
}
END_SECTION

START_SECTION(getSelectTransitionGroup())
{
  TEST_EQUAL(ptr->getSelectTransitionGroup(), true)
  ptr->setSelectTransitionGroup(false);
  TEST_EQUAL(ptr->getSelectTransitionGroup(), false)
}
END_SECTION

START_SECTION(getSegmentWindowLength())
{
  TEST_EQUAL(ptr->getSegmentWindowLength(), 8.0)
  ptr->setSegmentWindowLength(7.0);
  TEST_EQUAL(ptr->getSegmentWindowLength(), 7.0)
}
END_SECTION

START_SECTION(getSegmentStepLength())
{
  TEST_EQUAL(ptr->getSegmentStepLength(), 4.0)
  ptr->setSegmentStepLength(3.0);
  TEST_EQUAL(ptr->getSegmentStepLength(), 3.0)
}
END_SECTION

START_SECTION(getSelectHighestCount())
{
  TEST_EQUAL(ptr->getSelectHighestCount(), false)
  ptr->setSelectHighestCount(true);
  TEST_EQUAL(ptr->getSelectHighestCount(), true)
}
END_SECTION

START_SECTION(getVariableType())
{
  TEST_EQUAL(ptr->getVariableType(), "continuous")
  ptr->setVariableType("integer");
  TEST_EQUAL(ptr->getVariableType(), "integer")
}
END_SECTION

START_SECTION(getOptimalThreshold())
{
  TEST_EQUAL(ptr->getOptimalThreshold(), 0.5)
  ptr->setOptimalThreshold(0.6);
  TEST_EQUAL(ptr->getOptimalThreshold(), 0.6)
}
END_SECTION

START_SECTION(select_MRMFeature_qmip())
{
  ptr->select_MRMFeature_qmip(feature_map, targeted_exp);
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
