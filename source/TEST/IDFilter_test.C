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
// $Maintainer: Mathias Walzer $
// $Authors: Nico Pfeifer, Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <string>

#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <vector>

///////////////////////////

START_TEST(IDFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

///load input data
std::vector< ProteinIdentification > protein_identifications;
std::vector< PeptideIdentification > identifications;
String document_id;
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test.idXML"), protein_identifications, identifications, document_id);
PeptideIdentification identification = identifications[0];
ProteinIdentification protein_identification = protein_identifications[0];

/// Proteins for search
vector< FASTAFile::FASTAEntry > proteins;
proteins.push_back(FASTAFile::FASTAEntry("Q824A5", "test description 1", "LHASGITVTEIPVTATNFK"));
proteins.push_back(FASTAFile::FASTAEntry("Q872T5", "test description 2", "THPYGHAIVAGIERYPSK"));

IDFilter* ptr = 0;
IDFilter* nullPointer = 0;

START_SECTION((IDFilter()))
	ptr = new IDFilter();
  TEST_NOT_EQUAL(ptr, nullPointer);
END_SECTION

START_SECTION((~IDFilter()))
	delete ptr;
END_SECTION

START_SECTION((void filterIdentificationsByProteins(const ProteinIdentification& identification, const std::vector<FASTAFile::FASTAEntry> &proteins, ProteinIdentification& filtered_identification)))
	ProteinIdentification protein_identification2;

	IDFilter().filterIdentificationsByProteins(protein_identification, proteins, protein_identification2);

	TEST_EQUAL(protein_identification2.getScoreType() , "Mascot")
	TEST_EQUAL(protein_identification2.getHits().size(), 2)
	TEST_EQUAL(protein_identification2.getHits()[0].getAccession(), "Q824A5")
	TEST_EQUAL(protein_identification2.getHits()[1].getAccession(), "Q872T5")
END_SECTION

START_SECTION((void filterIdentificationsByProteins(const PeptideIdentification &identification, const std::vector< FASTAFile::FASTAEntry > &proteins, PeptideIdentification &filtered_identification, bool no_protein_identifiers=false)))
	PeptideIdentification identification2;

	IDFilter().filterIdentificationsByProteins(identification, proteins, identification2);

	TEST_EQUAL(identification2.getScoreType() , "Mascot")
	TEST_EQUAL(identification2.getHits().size(), 2)
	TEST_EQUAL(identification2.getHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(identification2.getHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
END_SECTION

START_SECTION((template <class IdentificationType> void filterIdentificationsByThreshold(const IdentificationType &identification, DoubleReal threshold_fraction, IdentificationType &filtered_identification)))
	PeptideIdentification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

	TEST_EQUAL(identification.getHits().size(), 10)
	IDFilter().filterIdentificationsByThreshold(identification, 1.3, identification2);
	peptide_hits = identification2.getHits();
	TEST_EQUAL(identification2.getScoreType() , "Mascot")

	TEST_EQUAL(peptide_hits.size(), 0)
	IDFilter().filterIdentificationsByThreshold(identification, 1.0, identification2);
	peptide_hits = identification2.getHits();
	TEST_EQUAL(peptide_hits.size(), 5)
	TEST_REAL_SIMILAR(peptide_hits[0].getScore() , 40)
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL((identification2.getHits()[0].getSequence()=="FINFGVNVEVLSRFQTK" && identification2.getHits()[1].getSequence()=="MSLLSNMISIVKVGYNAR") || (identification2.getHits()[0].getSequence()=="MSLLSNMISIVKVGYNAR" && identification2.getHits()[1].getSequence()=="FINFGVNVEVLSRFQTK") , true)
	TEST_REAL_SIMILAR(peptide_hits[1].getScore() , 40)
	TEST_EQUAL(peptide_hits[1].getRank() , 1)
	TEST_EQUAL(peptide_hits[2].getSequence() , "THPYGHAIVAGIERYPSK")
	TEST_REAL_SIMILAR(peptide_hits[2].getScore() , 39)
	TEST_EQUAL(peptide_hits[2].getRank() , 2)
	TEST_EQUAL(peptide_hits[3].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_SIMILAR(peptide_hits[3].getScore() , 34.85)
	TEST_EQUAL(peptide_hits[3].getRank() , 3)
	TEST_EQUAL(peptide_hits[4].getSequence() , "MRSLGYVAVISAVATDTDK")
	TEST_REAL_SIMILAR(peptide_hits[4].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[4].getRank() , 4)
END_SECTION

START_SECTION((template <class IdentificationType> void filterIdentificationsByScore(const IdentificationType &identification, DoubleReal threshold_score, IdentificationType &filtered_identification)))
	PeptideIdentification identification2;
	vector<PeptideHit> peptide_hits;

	TEST_EQUAL(identification.getHits().size(), 10)
	IDFilter().filterIdentificationsByScore(identification, 41, identification2);
	peptide_hits = identification2.getHits();
	TEST_EQUAL(identification2.getScoreType() , "Mascot")

	TEST_EQUAL(peptide_hits.size(), 0)
	IDFilter().filterIdentificationsByScore(identification, 33, identification2);
	peptide_hits = identification2.getHits();
	TEST_EQUAL(peptide_hits.size(), 5)
	TEST_REAL_SIMILAR(peptide_hits[0].getScore() , 40)
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL((identification2.getHits()[0].getSequence()=="FINFGVNVEVLSRFQTK" && identification2.getHits()[1].getSequence()=="MSLLSNMISIVKVGYNAR") || (identification2.getHits()[0].getSequence()=="MSLLSNMISIVKVGYNAR" && identification2.getHits()[1].getSequence()=="FINFGVNVEVLSRFQTK") , true)
	TEST_REAL_SIMILAR(peptide_hits[1].getScore() , 40)
	TEST_EQUAL(peptide_hits[1].getRank() , 1)
	TEST_EQUAL(peptide_hits[2].getSequence() , "THPYGHAIVAGIERYPSK")
	TEST_REAL_SIMILAR(peptide_hits[2].getScore() , 39)
	TEST_EQUAL(peptide_hits[2].getRank() , 2)
	TEST_EQUAL(peptide_hits[3].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_SIMILAR(peptide_hits[3].getScore() , 34.85)
	TEST_EQUAL(peptide_hits[3].getRank() , 3)
	TEST_EQUAL(peptide_hits[4].getSequence() , "MRSLGYVAVISAVATDTDK")
	TEST_REAL_SIMILAR(peptide_hits[4].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[4].getRank() , 4)
END_SECTION

START_SECTION((void filterIdentificationsByLength(const PeptideIdentification &identification, PeptideIdentification &filtered_identification, Size min_length, Size max_length)))
	PeptideIdentification identification_len(identification);
	AASequence eighter("OKTAMERR");
	AASequence niner("NONAMERRR");
	AASequence tener("DECAMERRRR");
	identification_len.insertHit(PeptideHit (99.99, 1, 2, eighter));
	identification_len.insertHit(PeptideHit (99.99, 1, 2, niner));
	identification_len.insertHit(PeptideHit (99.99, 1, 2, tener)) ;
	PeptideIdentification identification2;
	vector<PeptideHit> peptide_hits;

	TEST_EQUAL(identification_len.getHits().size(), 13)
	IDFilter().filterIdentificationsByLength(identification_len,identification2, 10);
	peptide_hits = identification2.getHits();
	TEST_EQUAL(peptide_hits.size(), 11)
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	for (Size i = 0; i < peptide_hits.size(); ++i)
	{
		TEST_EQUAL(peptide_hits[i].getSequence().size()>=10,true)
	}

	PeptideIdentification identification3;
	IDFilter().filterIdentificationsByLength(identification_len,identification3, 9,10);
	peptide_hits = identification3.getHits();
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	for (Size i = 0; i < peptide_hits.size(); ++i)
	{
		TEST_EQUAL(peptide_hits[i].getSequence().size()<=10,true)
	}
	for (Size i = 0; i < peptide_hits.size(); ++i)
	{
		TEST_EQUAL(peptide_hits[i].getSequence().size()>=9,true)
	}

	PeptideIdentification identification4;
	IDFilter().filterIdentificationsByLength(identification_len,identification4, 9,8);
	peptide_hits = identification4.getHits();
	TEST_EQUAL(peptide_hits.size(), 12)
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	for (Size i = 0; i < peptide_hits.size(); ++i)
	{
		TEST_EQUAL(peptide_hits[i].getSequence().size()>=9,true)
	}

END_SECTION

START_SECTION((void filterIdentificationsByExclusionPeptides(const PeptideIdentification &identification, const std::set< String > &peptides, PeptideIdentification &filtered_identification)))
	PeptideIdentification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;
	set<String> peptides;

	peptides.insert("LHASGITVTEIPVTATNFK");
	peptides.insert("MRSLGYVAVISAVATDTDK");
	peptides.insert("EGASTDFAALRTFLAEDGK");
	peptides.insert("DLEPGTDYEVTVSTLFGR");
	peptides.insert("FINFGVNVEVLSRFQTK");
	peptides.insert("MSLLSNMISIVKVGYNAR");
	peptides.insert("THPYGHAIVAGIERYPSK");
	peptides.insert("AITSDFANQAKTVLQNFK");

	IDFilter().filterIdentificationsByExclusionPeptides(identification, peptides, identification2);
	peptide_hits = identification2.getHits();
	TEST_EQUAL(identification2.getScoreType() , "Mascot")

	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence() , "TGCDTWGQGTLVTVSSASTK")
	TEST_REAL_SIMILAR(peptide_hits[0].getScore() , 10.93)
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL(peptide_hits[1].getSequence() , "TLCHHDATFDNLVWTPK")
	TEST_REAL_SIMILAR(peptide_hits[1].getScore() , 10.37)
	TEST_EQUAL(peptide_hits[1].getRank() , 2)
	protein_hits = protein_identification.getHits();
END_SECTION

START_SECTION((template<class PeakT> void filterIdentificationsByProteins(MSExperiment< PeakT > &experiment, const std::vector<FASTAFile::FASTAEntry> &proteins)))

	MSExperiment<> experiment;
  vector< FASTAFile::FASTAEntry > proteins;
  vector< PeptideIdentification > ids;
	PeptideIdentification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

  ids.push_back(identification);

  proteins.push_back(FASTAFile::FASTAEntry("Q824A5", "first desription", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(FASTAFile::FASTAEntry("Q872T5", "second description", "THPYGHAIVAGIERYPSK"));

	for (Size i = 0; i < 5; ++i)
	{
		experiment.addSpectrum(MSSpectrum<>());
	}
	experiment[3].setMSLevel(2);
	experiment[3].setPeptideIdentifications(ids);

	IDFilter().filterIdentificationsByProteins(experiment, proteins);

	identification2 = experiment[3].getPeptideIdentifications()[0];
	peptide_hits = identification2.getHits();
	TEST_EQUAL(identification2.getScoreType() , "Mascot")

	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_SIMILAR(peptide_hits[0].getScore() , 34.85)
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL(peptide_hits[1].getSequence() , "MRSLGYVAVISAVATDTDK")
	TEST_REAL_SIMILAR(peptide_hits[1].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[1].getRank() , 2)

END_SECTION

START_SECTION((void filterIdentificationsByBestHits(const PeptideIdentification& identification, PeptideIdentification& filtered_identification, bool strict = false)))
	PeptideIdentification identification2;

	//strict
	IDFilter().filterIdentificationsByBestHits(identification, identification2, true);
	TEST_EQUAL(identification2.getHits().size(), 0)
	TEST_EQUAL(identification2.getScoreType() , "Mascot")

	//not strict
	IDFilter().filterIdentificationsByBestHits(identification, identification2);
	TEST_EQUAL(identification2.getScoreType() , "Mascot")
	TEST_EQUAL(identification2.getHits().size(), 2)
	TEST_REAL_SIMILAR(identification2.getHits()[0].getScore() , 40)
	TEST_EQUAL(identification2.getHits()[0].getRank() , 1)
	TEST_REAL_SIMILAR(identification2.getHits()[1].getScore() , 40)
	TEST_EQUAL(identification2.getHits()[1].getRank() , 1)
	TEST_EQUAL((identification2.getHits()[0].getSequence()=="FINFGVNVEVLSRFQTK" && identification2.getHits()[1].getSequence()=="MSLLSNMISIVKVGYNAR") || (identification2.getHits()[0].getSequence()=="MSLLSNMISIVKVGYNAR" && identification2.getHits()[1].getSequence()=="FINFGVNVEVLSRFQTK") , true)
END_SECTION

START_SECTION((template <class PeakT> void filterIdentificationsByThresholds(MSExperiment< PeakT > &experiment, DoubleReal peptide_threshold_fraction, DoubleReal protein_threshold_fraction)))


	MSExperiment<> experiment;
  vector< PeptideIdentification > ids;
	PeptideIdentification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

  ids.push_back(identification);

	for (Size i = 0; i < 5; ++i)
	{
		experiment.addSpectrum(MSSpectrum<>());
	}
	experiment[3].setMSLevel(2);
	experiment[3].setPeptideIdentifications(ids);

	IDFilter().filterIdentificationsByThresholds(experiment, 1.0, 1.0);
	identification2 = experiment[3].getPeptideIdentifications()[0];
	peptide_hits = identification2.getHits();
	TEST_EQUAL(identification2.getScoreType() , "Mascot")

	TEST_EQUAL(peptide_hits.size(), 5)
	TEST_REAL_SIMILAR(peptide_hits[0].getScore() , 40)
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL((identification2.getHits()[0].getSequence()=="FINFGVNVEVLSRFQTK" && identification2.getHits()[1].getSequence()=="MSLLSNMISIVKVGYNAR") || (identification2.getHits()[0].getSequence()=="MSLLSNMISIVKVGYNAR" && identification2.getHits()[1].getSequence()=="FINFGVNVEVLSRFQTK") , true)
	TEST_REAL_SIMILAR(peptide_hits[1].getScore() , 40)
	TEST_EQUAL(peptide_hits[1].getRank() , 1)
	TEST_EQUAL(peptide_hits[2].getSequence() , "THPYGHAIVAGIERYPSK")
	TEST_REAL_SIMILAR(peptide_hits[2].getScore() , 39)
	TEST_EQUAL(peptide_hits[2].getRank() , 2)
	TEST_EQUAL(peptide_hits[3].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_SIMILAR(peptide_hits[3].getScore() , 34.85)
	TEST_EQUAL(peptide_hits[3].getRank() , 3)
	TEST_EQUAL(peptide_hits[4].getSequence() , "MRSLGYVAVISAVATDTDK")
	TEST_REAL_SIMILAR(peptide_hits[4].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[4].getRank() , 4)
END_SECTION

START_SECTION((template <class PeakT> void filterIdentificationsByScores(MSExperiment< PeakT > &experiment, DoubleReal peptide_threshold_score, DoubleReal protein_threshold_score)))


	MSExperiment<> experiment;
  vector< PeptideIdentification > ids;
	PeptideIdentification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

  ids.push_back(identification);

	for (Size i = 0; i < 5; ++i)
	{
		experiment.addSpectrum(MSSpectrum<>());
	}
	experiment[3].setMSLevel(2);
	experiment[3].setPeptideIdentifications(ids);

	IDFilter().filterIdentificationsByScores(experiment, 31.8621, 0);
	identification2 = experiment[3].getPeptideIdentifications()[0];
	peptide_hits = identification2.getHits();
	TEST_EQUAL(identification2.getScoreType() , "Mascot")

	TEST_EQUAL(peptide_hits.size(), 5)
	TEST_REAL_SIMILAR(peptide_hits[0].getScore() , 40)
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL((identification2.getHits()[0].getSequence()=="FINFGVNVEVLSRFQTK" && identification2.getHits()[1].getSequence()=="MSLLSNMISIVKVGYNAR") || (identification2.getHits()[0].getSequence()=="MSLLSNMISIVKVGYNAR" && identification2.getHits()[1].getSequence()=="FINFGVNVEVLSRFQTK") , true)
	TEST_REAL_SIMILAR(peptide_hits[1].getScore() , 40)
	TEST_EQUAL(peptide_hits[1].getRank() , 1)
	TEST_EQUAL(peptide_hits[2].getSequence() , "THPYGHAIVAGIERYPSK")
	TEST_REAL_SIMILAR(peptide_hits[2].getScore() , 39)
	TEST_EQUAL(peptide_hits[2].getRank() , 2)
	TEST_EQUAL(peptide_hits[3].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_SIMILAR(peptide_hits[3].getScore() , 34.85)
	TEST_EQUAL(peptide_hits[3].getRank() , 3)
	TEST_EQUAL(peptide_hits[4].getSequence() , "MRSLGYVAVISAVATDTDK")
	TEST_REAL_SIMILAR(peptide_hits[4].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[4].getRank() , 4)
END_SECTION

START_SECTION((template < class PeakT > void filterIdentificationsByBestNHits(MSExperiment< PeakT > &experiment, Size n)))
	MSExperiment<> experiment;
  vector< PeptideIdentification > ids;
	PeptideIdentification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

  ids.push_back(identification);

	for (Size i = 0; i < 5; ++i)
	{
		experiment.addSpectrum(MSSpectrum<>());
	}
	experiment[3].setMSLevel(2);
	experiment[3].setPeptideIdentifications(ids);

	IDFilter().filterIdentificationsByBestNHits(experiment, 3);
	identification2 = experiment[3].getPeptideIdentifications()[0];
	peptide_hits = identification2.getHits();
	TEST_EQUAL(identification2.getScoreType() , "Mascot")

	TEST_EQUAL(peptide_hits.size(), 3)
	TEST_REAL_SIMILAR(peptide_hits[0].getScore() , 40)
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL((identification2.getHits()[0].getSequence()=="FINFGVNVEVLSRFQTK" && identification2.getHits()[1].getSequence()=="MSLLSNMISIVKVGYNAR") || (identification2.getHits()[0].getSequence()=="MSLLSNMISIVKVGYNAR" && identification2.getHits()[1].getSequence()=="FINFGVNVEVLSRFQTK") , true)
	TEST_REAL_SIMILAR(peptide_hits[1].getScore() , 40)
	TEST_EQUAL(peptide_hits[1].getRank() , 1)
	TEST_EQUAL(peptide_hits[2].getSequence() , "THPYGHAIVAGIERYPSK")
	TEST_REAL_SIMILAR(peptide_hits[2].getScore() , 39)
	TEST_EQUAL(peptide_hits[2].getRank() , 2)

END_SECTION

START_SECTION((template < class IdentificationType > void filterIdentificationsByBestNHits(const IdentificationType &identification, Size n, IdentificationType &filtered_identification)))
	PeptideIdentification identification2;
	vector<PeptideHit> peptide_hits;

	IDFilter().filterIdentificationsByBestNHits(identification, 3, identification2);
	peptide_hits = identification2.getHits();
	TEST_EQUAL(identification2.getScoreType() , "Mascot")

	TEST_EQUAL(peptide_hits.size(), 3)
	TEST_REAL_SIMILAR(peptide_hits[0].getScore() , 40)
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL((identification2.getHits()[0].getSequence()=="FINFGVNVEVLSRFQTK" && identification2.getHits()[1].getSequence()=="MSLLSNMISIVKVGYNAR") || (identification2.getHits()[0].getSequence()=="MSLLSNMISIVKVGYNAR" && identification2.getHits()[1].getSequence()=="FINFGVNVEVLSRFQTK") , true)
	TEST_REAL_SIMILAR(peptide_hits[1].getScore() , 40)
	TEST_EQUAL(peptide_hits[1].getRank() , 1)
	TEST_EQUAL(peptide_hits[2].getSequence() , "THPYGHAIVAGIERYPSK")
	TEST_REAL_SIMILAR(peptide_hits[2].getScore() , 39)
	TEST_EQUAL(peptide_hits[2].getRank() , 2)

END_SECTION

START_SECTION((void filterIdentificationsByRTPValues(const PeptideIdentification &identification, PeptideIdentification &filtered_identification, DoubleReal p_value=0.05)))
	PeptideIdentification filtered_identification;
	String document_id;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test2.idXML"), protein_identifications, identifications, document_id);
	PeptideIdentification identification2 = identifications[0];
	ProteinIdentification protein_identification2 = protein_identifications[0];
	IDFilter().filterIdentificationsByRTPValues(identification2 , filtered_identification, 0.08);

	vector<PeptideHit> hits = filtered_identification.getHits();

	TEST_EQUAL(hits.size(), 4)
	TEST_EQUAL(hits[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(hits[1].getSequence(), "DLEPGTDYEVTVSTLFGR")
	TEST_EQUAL(hits[2].getSequence(), "FINFGVNVEVLSRFQTK")
	TEST_EQUAL(hits[3].getSequence(), "MSLLSNMISIVKVGYNAR")
END_SECTION

START_SECTION((void filterIdentificationsByRTFirstDimPValues(const PeptideIdentification &identification, PeptideIdentification &filtered_identification, DoubleReal p_value=0.05)))
	PeptideIdentification filtered_identification;
	String document_id;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test3.idXML"), protein_identifications, identifications, document_id);
	PeptideIdentification identification2 = identifications[0];
	ProteinIdentification protein_identification2 = protein_identifications[0];
	IDFilter().filterIdentificationsByRTFirstDimPValues(identification2 , filtered_identification, 0.08);

	vector<PeptideHit> hits = filtered_identification.getHits();

	TEST_EQUAL(hits.size(), 4)
	TEST_EQUAL(hits[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(hits[1].getSequence(), "DLEPGTDYEVTVSTLFGR")
	TEST_EQUAL(hits[2].getSequence(), "FINFGVNVEVLSRFQTK")
	TEST_EQUAL(hits[3].getSequence(), "MSLLSNMISIVKVGYNAR")
END_SECTION

START_SECTION((void removeUnreferencedProteinHits(const ProteinIdentification &identification, const std::vector< PeptideIdentification > peptide_identifications, ProteinIdentification &filtered_identification)))
	String document_id;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test4.idXML"), protein_identifications, identifications, document_id);

	ProteinIdentification protein_identification;
  IDFilter().removeUnreferencedProteinHits(protein_identifications[0], identifications, protein_identification);

  TEST_EQUAL(protein_identification.getHits().size(), 3)
  TEST_EQUAL(protein_identification.getHits()[0].getAccession(), "Q824A5")
  TEST_EQUAL(protein_identification.getHits()[1].getAccession(), "S53854")
  TEST_EQUAL(protein_identification.getHits()[2].getAccession(), "Q872T5")
END_SECTION

START_SECTION((void filterIdentificationsUnique(const PeptideIdentification &identification, PeptideIdentification &filtered_identification)))
	PeptideIdentification id, id2;
	vector<PeptideHit> hits;
	PeptideHit hit;
	hit.setSequence("DFPIANGER");
	hit.setCharge(1);
	hit.setScore(0.3);
	hits.push_back(hit);
	hit.setCharge(2);
	hits.push_back(hit);
	hit.setScore(0.5);
	hits.push_back(hit);
	hit.setSequence("DFPIANGEK");
	hits.push_back(hit);
	hits.push_back(hit);
	hits.push_back(hit);
	hit.setCharge(5);
	hits.push_back(hit);
	IDFilter id_filter;
	TEST_EQUAL(hits.size(), 7)
	id.setHits(hits);

	id_filter.filterIdentificationsUnique(id, id2);
	TEST_EQUAL(id2.getHits().size(), 5)
	TEST_STRING_EQUAL(id2.getHits()[3].getSequence().toString(), "DFPIANGEK")
	TEST_EQUAL(id2.getHits()[3].getCharge(), 2)
	TEST_STRING_EQUAL(id2.getHits()[4].getSequence().toString(), "DFPIANGEK")
	TEST_EQUAL(id2.getHits()[4].getCharge(), 5)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
