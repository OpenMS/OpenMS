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
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/QC/MzCalibration.h>

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

//

// peptide_ID>1
// [0]: sequence= "AATT" charge=2
// RT fmap <0.5
// fmap &/or exp empty
//keine richtige seq

// Dummy peakmap
PeakMap exp;
MSSpectrum spec;
Precursor pre;
std::vector< MSSpectrum> spectra;

//Dummy FeatureMap
FeatureMap fmap;
PeptideHit peptide_hit;
std::vector<PeptideHit> peptide_hits;
PeptideIdentification peptide_ID;
vector<PeptideIdentification> identifications;
vector<PeptideIdentification> unassignedIDs;
Feature feature1;

// Spectra: [0]: MZ_raw=5.5;6.5 RT=0 MS-Level=1, [1]: MZ_raw=5.55, 6.55 RT=0.5 MS-Level=2, #precursor=2/1, #spectrum=2
pre.setMetaValue("mz_raw", 5.55);
spec.setMSLevel(2);
spec.setPrecursors({ pre });
spec.setRT(0);
spectra.push_back(spec);
pre.setMetaValue("mz_raw", 6.55);
spec.setPrecursors({ pre });
spec.setRT(0.5);
spectra.push_back(spec);
pre.setMetaValue("mz_raw", 7.55);
spec.setPrecursors({ pre });
spec.setRT(1);
spectra.push_back(spec);
exp.setSpectra(spectra);

//FeatureMap:
//Feature:
peptide_hit.setSequence(AASequence::fromString("AAAA"));
peptide_hit.setCharge(2);
peptide_hits.push_back(peptide_hit);
peptide_hit.setSequence(AASequence::fromString("TTTT"));
peptide_hit.setCharge(1);
peptide_hits.push_back(peptide_hit);
peptide_ID.setHits(peptide_hits);
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
		std::cout << "\n1\n";
		TEST_EQUAL(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_raw"), 5.55);
		std::cout << "\n2\n";
		TEST_EQUAL(fmap[0].getPeptideIdentifications()[1].getHits()[0].getMetaValue("mz_raw"), 7.55);
		
		//test unassigned
		TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_raw"), 6.55);
		//test refMZ
		TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_ref"), AASequence::fromString("AAAA").getMonoWeight(OpenMS::Residue::Full, 2) / 2);
		TEST_REAL_SIMILAR(fmap[0].getPeptideIdentifications()[1].getHits()[0].getMetaValue("mz_ref"), AASequence::fromString("WWWW").getMonoWeight(OpenMS::Residue::Full, 3) / 3);
		TEST_REAL_SIMILAR(fmap.getUnassignedPeptideIdentifications()[0].getHits()[0].getMetaValue("mz_ref"), AASequence::fromString("YYYY").getMonoWeight(OpenMS::Residue::Full, 2) / 2);
		
		// empty MSExperiment
		MSExperiment exp_empty{};
		TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, cal.compute(fmap, exp_empty), "The PeakMap is empty.");
		// empty FeatureMap
		FeatureMap fmap_empty{};
		MzCalibration cal2;
		cal2.compute(fmap_empty, exp);
		TEST_EQUAL(fmap_empty.isMetaEmpty(), true);

		//RT feature map > RT msexpriment
		PeakMap exp2;
		pre.setMetaValue("mz_raw", 4.4);
		spec.setMSLevel(2);
		spec.setPrecursors({ pre });
		spec.setRT(4);
		spectra.push_back(spec);
		exp2.setSpectra(spectra);
		//Dummy FeatureMap
		FeatureMap fmap2;
		peptide_hit.setSequence(AASequence::fromString("AAAA"));
		peptide_hit.setCharge(2);
		peptide_hits.push_back(peptide_hit);
		peptide_ID.setHits(peptide_hits);
		peptide_ID.setRT(5);
		identifications.push_back(peptide_ID);
		feature1.setPeptideIdentifications(identifications);
		fmap2.push_back(feature1);
		TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, cal.compute(fmap2, exp2), "The retention time of the MZML and featureXML file does not match.");
		//RT difference between feature map and RT msexpriment is greater than EPSILON (error tolerance)
		PeakMap exp3;
		pre.setMetaValue("mz_raw", 4.4);
		spec.setMSLevel(2);
		spec.setPrecursors({ pre });
		spec.setRT(4);
		spectra.push_back(spec);
		pre.setMetaValue("mz_raw", 2.4);
		spec.setPrecursors({ pre });
		spec.setRT(5);
		spectra.push_back(spec);
		exp3.setSpectra(spectra);
		//Dummy FeatureMap
		FeatureMap fmap3;
		peptide_ID.setRT(4.1);
		identifications.clear();
		identifications.push_back(peptide_ID);
		feature1.setPeptideIdentifications(identifications);
		fmap3.push_back(feature1);
		TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, cal.compute(fmap3, exp3), "The retention time of the MZML and featureXML file does not match.");

}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST