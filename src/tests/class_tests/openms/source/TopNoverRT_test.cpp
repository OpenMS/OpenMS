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
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/QC/TopNoverRT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/METADATA/DataProcessing.h>

///////////////////////////

START_TEST(TopNoverRT, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TopNoverRT* ptr = nullptr;
TopNoverRT* nullPointer = nullptr;

START_SECTION(TopNoverRT())
{
	ptr = new TopNoverRT;
	TEST_NOT_EQUAL(ptr, nullPointer);
}
END_SECTION

START_SECTION(~TopNoverRT())
{
	delete ptr;
}
END_SECTION

START_SECTION(QCBase::Status requires() const override)
{
	TopNoverRT top;
	TEST_EQUAL(top.requires() == (QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT), true);
}
END_SECTION

START_SECTION(compute(const MSExperiment& exp, FeatureMap& features))
{
	//Valid FeatureMap
	FeatureMap fmap;
	PeptideIdentification peptide_ID;
	vector<PeptideIdentification> identifications;
	vector<PeptideIdentification> unassignedIDs;
	Feature feature1;
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
	feature1.setPeptideIdentifications(identifications);
	fmap.push_back(feature1);
	//unassigned PeptideHits
	peptide_ID.setRT(1.5);
	unassignedIDs.push_back(peptide_ID);
	peptide_ID.setRT(2.5);
	unassignedIDs.push_back(peptide_ID);
	fmap.setUnassignedPeptideIdentifications(unassignedIDs);
	
	//MSExperiment
	PeakMap exp;
	MSSpectrum spec;
	std::vector< MSSpectrum> spectra;
	spec.setMSLevel(2);
	spec.setRT(0);
	spectra.push_back(spec);
	spec.setMSLevel(1);
	spec.setRT(0.5);
	spectra.push_back(spec);
	spec.setMSLevel(2);
	spec.setRT(1);
	spectra.push_back(spec);
	spec.setRT(1.5);
	spectra.push_back(spec);
	spec.setRT(2.5);
	spectra.push_back(spec);
	spec.setMSLevel(1);
	spec.setRT(9);
	spectra.push_back(spec);
	spec.setMSLevel(2);
	spec.setRT(10);
	spectra.push_back(spec);
	spec.setRT(12);
	spectra.push_back(spec);
	//not identified
	spec.setRT(20);
	spectra.push_back(spec);
	exp.setSpectra(spectra);
	
	TopNoverRT top;
	top.compute(exp, fmap);

	//test features
	TEST_EQUAL(fmap[0].getPeptideIdentifications()[0].getMetaValue("ScanEventNumber"), 1);
	TEST_EQUAL(fmap[0].getPeptideIdentifications()[0].getMetaValue("identified"), '+');
	TEST_EQUAL(fmap[0].getPeptideIdentifications()[1].getMetaValue("ScanEventNumber"), 1);
	TEST_EQUAL(fmap[1].getPeptideIdentifications()[0].getMetaValue("ScanEventNumber"), 1);
	TEST_EQUAL(fmap[1].getPeptideIdentifications()[1].getMetaValue("ScanEventNumber"), 2);
	//test unassigned
	TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[0].getMetaValue("ScanEventNumber"), 2);
	TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[0].getMetaValue("identified"), '+');
	TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[1].getMetaValue("ScanEventNumber"), 3);
	TEST_REAL_SIMILAR(fmap.getUnassignedPeptideIdentifications()[2].getRT(), 20);
	TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[2].getMetaValue("ScanEventNumber"), 3);
	TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[2].getMetaValue("identified"), '-');

	//empty FeatureMap
	FeatureMap fmap_empty{};
	top.compute(exp, fmap_empty);
	TEST_EQUAL(fmap_empty.getUnassignedPeptideIdentifications().size(), 7);
	//empty feature
	fmap_empty.clear();
	Feature feature_empty{};
	fmap_empty.push_back(feature_empty);
	top.compute(exp, fmap_empty);
	TEST_EQUAL(fmap_empty.getUnassignedPeptideIdentifications().size(), 7);
	//empty PeptideIdentifications
	identifications.clear();
	fmap_empty.clear();
	feature_empty.setPeptideIdentifications(identifications);
	fmap_empty.setUnassignedPeptideIdentifications(identifications);
	fmap_empty.push_back(feature_empty);
	top.compute(exp, fmap_empty);
	TEST_EQUAL(fmap_empty.getUnassignedPeptideIdentifications().size(), 7);
	//empty MSExperiment
	PeakMap exp_empty{};
	TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, top.compute(exp_empty, fmap), "The mzml file / MSExperiment is empty.\n");

	//test exceptions spectrum.getRT() - peptide_ID.getRT() > EPSILON_
	exp.getSpectra()[0].setRT(0.1);
	TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, top.compute(exp, fmap), "PeptideID with RT " + to_string(0.0) + " s does not have a matching MS2 spectrum. Closest RT was " + to_string(0.1) + ", which seems too far off.\n");
	//test exception rt>end()
	exp.getSpectra()[0].setRT(0);
	fmap[1].getPeptideIdentifications()[1].setRT(50);
	TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, top.compute(exp, fmap), "The retention time of the MZML and featureXML file does not match.");
	//test exception if closest RT to PeptideID has MS-Level=1
	exp.getSpectra()[0].setMSLevel(1);
	TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, top.compute(exp, fmap), "The matching retention time of the MZML has the wrong MSLevel");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

