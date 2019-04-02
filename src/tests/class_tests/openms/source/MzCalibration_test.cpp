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
// $Maintainer: Chris Bielow$
// $Authors: Juliane Schmachtenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
///////////////////////////
#include <OpenMS/QC/MzCalibration.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
///////////////////////////

START_TEST(MzCalibration, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
using namespace OpenMS;
using namespace std;

MzCalibration* ptr = nullptr;
MzCalibration* nullPointer = nullptr;
START_SECTION(MzCalibration())
ptr = new MzCalibration();
TEST_NOT_EQUAL(ptr, nullPointer);
END_SECTION

START_SECTION(~MzCalibration())
delete ptr;
END_SECTION

START_SECTION(QCBase::Status requires() const override)
{
		MzCalibration mzCal;
		TEST_EQUAL(mzCal.requires() == (QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT), true);
}
END_SECTION

// PeakMap
PeakMap exp;
MSSpectrum spec;
Precursor pre;
std::vector< MSSpectrum> spectra;
pre.setMetaValue("mz_raw", 5);
spec.setMSLevel(2);
spec.setPrecursors({ pre });
spec.setRT(0);
spectra.push_back(spec);
pre.setMetaValue("mz_raw", 6);
spec.setPrecursors({ pre });
spec.setRT(0.5);
spectra.push_back(spec);
pre.setMetaValue("mz_raw", 7);
spec.setPrecursors({ pre });
spec.setRT(1);
spectra.push_back(spec);
exp.setSpectra(spectra);

//FeatureMap
FeatureMap fmap;
PeptideHit peptide_hit;
std::vector<PeptideHit> peptide_hits;
PeptideIdentification peptide_ID;
vector<PeptideIdentification> identifications;
vector<PeptideIdentification> unassignedIDs;
Feature feature1;
peptide_hit.setSequence(AASequence::fromString("AAAA"));
peptide_hit.setCharge(2);
peptide_hits.push_back(peptide_hit);
peptide_ID.setHits(peptide_hits);
peptide_ID.setMZ(5.5);
peptide_hits.clear();
peptide_ID.setRT(0);
identifications.push_back(peptide_ID);
peptide_hit.setSequence(AASequence::fromString("WWWW"));
peptide_hit.setCharge(3);
peptide_hits.push_back(peptide_hit);
peptide_ID.setHits(peptide_hits);
peptide_hits.clear();
peptide_ID.setRT(1);
identifications.push_back(peptide_ID);
feature1.setPeptideIdentifications(identifications);
fmap.push_back(feature1);
//unassigned PeptideHits
peptide_hit.setSequence(AASequence::fromString("YYYY"));
peptide_hit.setCharge(2);
peptide_hits.push_back(peptide_hit);
peptide_ID.setHits(peptide_hits);
peptide_hits.clear();
peptide_ID.setRT(0.5);
unassignedIDs.push_back(peptide_ID);
fmap.setUnassignedPeptideIdentifications(unassignedIDs);

MzCalibration cal;
//tests compute function
START_SECTION(void compute(FeatureMap& features, const MSExperiment& exp))
{
	cal.compute(fmap, exp);
	for (Feature f:fmap)
	{
		for (auto pepID : f.getPeptideIdentifications())
		{
			ABORT_IF(!pepID.getHits()[0].metaValueExists("mz_raw"));
		}
	}
	//test with valid input
	TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_raw"), 5);
	TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[1].getHits()[0].getMetaValue("mz_raw"), 7);
	//test unassigned
	TEST_REAL_SIMILAR(fmap.getUnassignedPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_raw"), 6);
	
	//test refMZ
	double ref = AASequence::fromString("AAAA").getMonoWeight(OpenMS::Residue::Full, 2) / 2;
	TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_ref"), ref);
	TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[1].getHits()[0].getMetaValue("mz_ref"), AASequence::fromString("WWWW").getMonoWeight(OpenMS::Residue::Full, 3) / 3);
	TEST_REAL_SIMILAR(fmap.getUnassignedPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_ref"), AASequence::fromString("YYYY").getMonoWeight(OpenMS::Residue::Full, 2) / 2);
	
	//test  mz_error
	TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("uncalibrated_mz_error_ppm"),(5-ref)/ref*1000000);
	TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("calibrated_mz_error_ppm"), (5.5-ref)/ref*1000000);
	
	// test exception for empty MSExperiment
	MSExperiment exp_empty{};
	TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, cal.compute(fmap, exp_empty), "The PeakMap is empty.");
	
	// test empty FeatureMap
	FeatureMap fmap_empty{};
	cal.compute(fmap_empty, exp);
	TEST_EQUAL(fmap_empty.isMetaEmpty(), true);
	
	//test feature is empty
	Feature feature_empty{};
	fmap_empty.push_back(feature_empty);
	fmap_empty.push_back(feature1);
	cal.compute(fmap_empty, exp);
	std::cout << "\n Compute empty Feature without error.\n";

	//test empty PeptideIdentification
	fmap_empty.clear();
	PeptideIdentification peptide_ID_empty {};
	identifications.push_back(peptide_ID_empty);
	feature1.setPeptideIdentifications(identifications);
	fmap.push_back(feature1);
	cal.compute(fmap_empty, exp);
	std::cout << "\n Compute empty PeptideIdentification without error.\n";

	//test empty hit
	peptide_ID.setHits(std::vector<PeptideHit>{});
	identifications.clear();
	identifications.push_back(peptide_ID);
	fmap.push_back(feature1);
	cal.compute(fmap_empty, exp);
	std::cout << "\n Compute empty PeptideHit without error.\n";

	//test exception if RT-values of the feature map are greater than the greatest RT-values in MSexpriment
	pre.setMetaValue("mz_raw", 4.4);
	spec.setPrecursors({ pre });
	spec.setRT(4);
	spectra.push_back(spec);
	exp.setSpectra(spectra);
	
	peptide_hit.setCharge(2);
	peptide_hits.push_back(peptide_hit);
	peptide_ID.setHits(peptide_hits);
	peptide_ID.setRT(5);
	identifications.push_back(peptide_ID);
	feature1.setPeptideIdentifications(identifications);
	fmap.clear();
	fmap.push_back(feature1);
	TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, cal.compute(fmap, exp), "The retention time of the MZML and featureXML file does not match.");

	//test exception if RT difference between mz_raw (get from MSexpriment) and mz after calibration (in FeatureMap) is greater than EPSILON_
	pre.setMetaValue("mz_raw", 4.4);
	spec.setPrecursors({ pre });
	spec.setRT(4);
	spectra.push_back(spec);
	pre.setMetaValue("mz_raw", 2.4);
	spec.setPrecursors({ pre });
	spec.setRT(5);
	spectra.push_back(spec);
	exp.setSpectra(spectra);

	peptide_ID.setRT(4.1);
	identifications.clear();
	identifications.push_back(peptide_ID);
	feature1.setPeptideIdentifications(identifications);
	fmap.clear();
	fmap.push_back(feature1);
	TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, cal.compute(fmap, exp), "PeptideID with RT " + to_string(4.1) + " s does not have a matching MS2 spectrum. Closest RT was " + to_string(5.0) + ", which seems too far off.\n");

	//test exception if no calibrtion was performed and the meta_value mz_raw does not exist
	Precursor pre_without_meta;
	pre_without_meta.setMetaValue("dummy", 1);
	spec.setPrecursors({ pre_without_meta });
	spec.setRT(4);
	spectra.clear();
	spectra.push_back(spec);
	exp.setSpectra(spectra);
	
	peptide_ID.setRT(4.01);
	identifications.clear();
	identifications.push_back(peptide_ID);
	feature1.setPeptideIdentifications(identifications);
	fmap.clear();
	fmap.push_back(feature1);
	TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, cal.compute(fmap, exp), "InternalCalibration was not called; MetaValue 'mz_raw' missing from MSExperiment.");
	//test wrong MS-Level
	spec.setMSLevel(1);
	spec.setPrecursors({ pre });
	spec.setRT(40);
	spectra.clear();
	spectra.push_back(spec);
	exp.setSpectra(spectra);
	peptide_ID.setRT(40);
	identifications.clear();
	identifications.push_back(peptide_ID);
	feature1.setPeptideIdentifications(identifications);
	fmap.clear();
	fmap.push_back(feature1);
	TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, cal.compute(fmap, exp), "The matching retention time of the MZML has the wrong MSLevel");
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST