// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Chris Bielow
// $Authors: Juliane Schmachtenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/QC/RTAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
///////////////////////////

START_TEST(RTAlignment, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

RTAlignment* ptr = nullptr;
RTAlignment* nullPointer = nullptr;

START_SECTION(RTAlignment())
{
  ptr = new RTAlignment;
  TEST_NOT_EQUAL(ptr, nullPointer);
}
END_SECTION

START_SECTION(~RTAlignment())
{
  delete ptr;
}
END_SECTION

RTAlignment rtA;
START_SECTION(QCBase::Status requires() const override)
{
  TEST_EQUAL(rtA.requires() == (QCBase::Status() | QCBase::Requires::TRAFOALIGN | QCBase::Requires::POSTFDRFEAT), true);
}
END_SECTION

START_SECTION(const String& getName() const override)
{
  TEST_EQUAL(rtA.getName(), "RTAlignment")
}
END_SECTION

START_SECTION((compute(FeatureMap& features, TransformationDescription& trafo)))
{
  //Valid FeatureMap
  FeatureMap fmap;
  PeptideIdentification peptide_ID;
  vector<PeptideIdentification> identifications;
  vector<PeptideIdentification> unassignedIDs;
  Feature feature1, feature2;
  peptide_ID.setRT(0);
  identifications.push_back(peptide_ID);
  peptide_ID.setRT(1);
  identifications.push_back(peptide_ID);
  feature1.setPeptideIdentifications(identifications);
  identifications.clear();
  fmap.push_back(feature1);
  peptide_ID.setRT(10);
  identifications.push_back(peptide_ID);
  peptide_ID.setRT(12);
  identifications.push_back(peptide_ID);
  feature2.setPeptideIdentifications(identifications);
  fmap.push_back(feature2);
  //unassigned PeptideHits
  peptide_ID.setRT(0.5);
  unassignedIDs.push_back(peptide_ID);
  peptide_ID.setRT(2.5);
  unassignedIDs.push_back(peptide_ID);
  fmap.setUnassignedPeptideIdentifications(unassignedIDs);

  //Transformation
  TransformationDescription td;
  td.fitModel("identity", Param());
  td.setDataPoints(vector<pair<double, double> > { { 0.0, 1.0 } , { 0.25, 1.5 }, { 0.5, 2.0 }, { 1.0, 3.0 } });
  td.fitModel("linear");
  RTAlignment rtA;
  rtA.compute(fmap, td);
  //test features
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getMetaValue("rt_align"), 1);
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getMetaValue("rt_raw"), 0);
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[1].getMetaValue("rt_align"), 3);
  TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[1].getMetaValue("rt_raw"), 1);
  TEST_REAL_SIMILAR(fmap[1].getPeptideIdentifications()[0].getMetaValue("rt_align"), 21);
  TEST_REAL_SIMILAR(fmap[1].getPeptideIdentifications()[0].getMetaValue("rt_raw"), 10);
  //test unassigned
  TEST_REAL_SIMILAR(fmap.getUnassignedPeptideIdentifications()[0].getMetaValue("rt_align"), 2);
  TEST_REAL_SIMILAR(fmap.getUnassignedPeptideIdentifications()[0].getMetaValue("rt_raw"), 0.5);
  TEST_REAL_SIMILAR(fmap.getUnassignedPeptideIdentifications()[1].getMetaValue("rt_align"), 6);
  TEST_REAL_SIMILAR(fmap.getUnassignedPeptideIdentifications()[1].getMetaValue("rt_raw"), 2.5);

  //empty FeatureMap
  FeatureMap fmap_empty{};
  rtA.compute(fmap_empty, td);
  //empty feature
  Feature feature_empty{};
  fmap_empty.push_back(feature_empty);
  rtA.compute(fmap_empty, td);
  //empty PeptideIdentifications
  identifications.clear();
  feature1.setPeptideIdentifications(identifications);
  fmap_empty.push_back(feature1);
  rtA.compute(fmap_empty, td);
  //empty PeptideIdentification
  PeptideIdentification peptide_ID_empty{};
  identifications.push_back(peptide_ID_empty);
  feature1.setPeptideIdentifications(identifications);
  fmap_empty.push_back(feature1);
  rtA.compute(fmap_empty, td);

  //data processing: after alignment
  DataProcessing processing_method{};
  std::set<DataProcessing::ProcessingAction> dp{ DataProcessing::ProcessingAction::ALIGNMENT };
  processing_method.setProcessingActions(dp);
  fmap.setDataProcessing({ processing_method });
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, rtA.compute(fmap, td), "Metric RTAlignment received a featureXML AFTER map alignment, but needs a featureXML BEFORE map alignment!");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

