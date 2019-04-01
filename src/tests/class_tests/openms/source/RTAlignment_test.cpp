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
// $Maintainer: Timo Sachsenberg$
// $Authors: Chris Bielow, Juliane Schmachtenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/QC/RTAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/METADATA/DataProcessing.h>

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

START_SECTION(QCBase::Status requires() const override)
{
		RTAlignment rtA;
		TEST_EQUAL(rtA.requires() == (QCBase::Status() | QCBase::Requires::TRAFOALIGN | QCBase::Requires::PREALIGNFEAT), true);
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
		//data processing
		DataProcessing processing_method{};
		std::set<DataProcessing::ProcessingAction> dp{ DataProcessing::ProcessingAction::ALIGNMENT };
		processing_method.setProcessingActions(dp);
		fmap.setDataProcessing({ processing_method });

		//Transformation


}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

