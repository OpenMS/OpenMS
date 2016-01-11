// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Authors: Nico Pfeifer, Mathias Walzer, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

START_TEST(IDFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// load input data
vector<ProteinIdentification> proteins;
vector<PeptideIdentification> peptides;
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test.idXML"), proteins,
                 peptides);
peptides[0].sort();

IDFilter* ptr = 0;
IDFilter* nullPointer = 0;

START_SECTION((IDFilter()))
  ptr = new IDFilter();
  TEST_NOT_EQUAL(ptr, nullPointer);
END_SECTION

START_SECTION((~IDFilter()))
  delete ptr;
END_SECTION

START_SECTION((template <class IdentificationType> static void keepHitsMatchingProteins(std::vector<IdentificationType>& ids, const set<String> accessions)))
{
  set<String> accessions;
  accessions.insert("Q824A5");
  accessions.insert("Q872T5");

  vector<ProteinIdentification> proteins_copy = proteins;
  IDFilter::keepHitsMatchingProteins(proteins_copy, accessions);

  TEST_EQUAL(proteins_copy[0].getScoreType(), "Mascot");
  TEST_EQUAL(proteins_copy[0].getHits().size(), 2);
  TEST_EQUAL(proteins_copy[0].getHits()[0].getAccession(), "Q824A5");
  TEST_EQUAL(proteins_copy[0].getHits()[1].getAccession(), "Q872T5");

  vector<PeptideIdentification> peptides_copy = peptides;
  IDFilter::keepHitsMatchingProteins(peptides_copy, accessions);

  TEST_EQUAL(peptides_copy[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptides_copy[0].getHits().size(), 2);
  TEST_EQUAL(peptides_copy[0].getHits()[0].getSequence(),
             AASequence::fromString("LHASGITVTEIPVTATNFK"));
  TEST_EQUAL(peptides_copy[0].getHits()[1].getSequence(),
             AASequence::fromString("MRSLGYVAVISAVATDTDK"));
}
END_SECTION

START_SECTION((template <class IdentificationType> static void filterHitsBySignificance(vector<IdentificationType>& ids, double threshold_fraction = 1.0)))
{
  vector<PeptideIdentification> peptides_copy = peptides;
  vector<PeptideHit>& peptide_hits = peptides_copy[0].getHits();
  TEST_EQUAL(peptide_hits.size(), 11);

  IDFilter::filterHitsBySignificance(peptides_copy, 1.0);
  TEST_EQUAL(peptide_hits.size(), 5);
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40);
  TEST_EQUAL(peptide_hits[0].getSequence(),
             AASequence::fromString("FINFGVNVEVLSRFQTK"));
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40);
  TEST_EQUAL(peptide_hits[1].getSequence(),
             AASequence::fromString("MSLLSNMISIVKVGYNAR"));
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 39);
  TEST_EQUAL(peptide_hits[2].getSequence(),
             AASequence::fromString("THPYGHAIVAGIERYPSK"));
  TEST_REAL_SIMILAR(peptide_hits[3].getScore(), 34.85);
  TEST_EQUAL(peptide_hits[3].getSequence(),
             AASequence::fromString("LHASGITVTEIPVTATNFK"));
  TEST_REAL_SIMILAR(peptide_hits[4].getScore(), 33.85);
  TEST_EQUAL(peptide_hits[4].getSequence(),
             AASequence::fromString("MRSLGYVAVISAVATDTDK"));

  IDFilter::filterHitsBySignificance(peptides_copy, 1.3);
  TEST_EQUAL(peptides_copy[0].getScoreType() , "Mascot")
  TEST_EQUAL(peptide_hits.size(), 0);
}
END_SECTION

START_SECTION((template <class IdentificationType> static void filterHitsByScore(vector<IdentificationType>& ids, double threshold_score)))
{
  vector<PeptideIdentification> peptides_copy = peptides;
  vector<PeptideHit>& peptide_hits = peptides_copy[0].getHits();
  TEST_EQUAL(peptide_hits.size(), 11);

  IDFilter::filterHitsByScore(peptides_copy, 33);
  TEST_EQUAL(peptide_hits.size(), 5);
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40);
  TEST_EQUAL(peptide_hits[0].getSequence(),
             AASequence::fromString("FINFGVNVEVLSRFQTK"));
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40);
  TEST_EQUAL(peptide_hits[1].getSequence(),
             AASequence::fromString("MSLLSNMISIVKVGYNAR"));
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 39);
  TEST_EQUAL(peptide_hits[2].getSequence(),
             AASequence::fromString("THPYGHAIVAGIERYPSK"));
  TEST_REAL_SIMILAR(peptide_hits[3].getScore(), 34.85);
  TEST_EQUAL(peptide_hits[3].getSequence(),
             AASequence::fromString("LHASGITVTEIPVTATNFK"));
  TEST_REAL_SIMILAR(peptide_hits[4].getScore(), 33.85);
  TEST_EQUAL(peptide_hits[4].getSequence(),
             AASequence::fromString("MRSLGYVAVISAVATDTDK"));

  IDFilter::filterHitsByScore(peptides_copy, 41);
  TEST_EQUAL(peptides_copy[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 0);
}
END_SECTION

START_SECTION((static void filterPeptidesByLength(vector<PeptideIdentification>& peptides, Size min_length, Size max_length = UINT_MAX)))
{
  vector<PeptideIdentification> peptides_copy = peptides;
  AASequence eighter = AASequence::fromString("OKTAMERR");
  AASequence niner = AASequence::fromString("NONAMERRR");
  AASequence tener = AASequence::fromString("DECAMERRRR");
  peptides_copy[0].insertHit(PeptideHit(99.99, 1, 2, eighter));
  peptides_copy[0].insertHit(PeptideHit(99.99, 1, 2, niner));
  peptides_copy[0].insertHit(PeptideHit(99.99, 1, 2, tener));
  TEST_EQUAL(peptides_copy[0].getHits().size(), 14);
  
  vector<PeptideIdentification> peptides_copy2 = peptides_copy;
  vector<PeptideHit>& peptide_hits = peptides_copy2[0].getHits();
  IDFilter::filterPeptidesByLength(peptides_copy2, 10);
  TEST_EQUAL(peptide_hits.size(), 12)
  for (Size i = 0; i < peptide_hits.size(); ++i)
  {
    TEST_EQUAL(peptide_hits[i].getSequence().size() >= 10, true);
  }

  peptides_copy2 = peptides_copy;
  IDFilter::filterPeptidesByLength(peptides_copy2, 9, 10);
  TEST_EQUAL(peptide_hits.size(), 2);
  for (Size i = 0; i < peptide_hits.size(); ++i)
  {
    TEST_EQUAL(peptide_hits[i].getSequence().size() >= 9, true);
    TEST_EQUAL(peptide_hits[i].getSequence().size() <= 10, true);
  }

  peptides_copy2 = peptides_copy;
  IDFilter::filterPeptidesByLength(peptides_copy2, 9, 8);
  TEST_EQUAL(peptide_hits.size(), 13)
  for (Size i = 0; i < peptide_hits.size(); ++i)
  {
    TEST_EQUAL(peptide_hits[i].getSequence().size() >= 9, true);
  }
}
END_SECTION

START_SECTION((static void removePeptidesWithMatchingSequences(vector<PeptideIdentification>& peptides, const vector<PeptideIdentification>& bad_peptides, bool ignore_mods)))
{
  vector<PeptideIdentification> peptides_copy = peptides;
  vector<PeptideHit>& peptide_hits = peptides_copy[0].getHits();
  vector<PeptideIdentification> bad_peptides(1);
  vector<PeptideHit>& bad_hits = bad_peptides[0].getHits();
  bad_hits.resize(8);
  bad_hits[0].setSequence(AASequence::fromString("LHASGITVTEIPVTATNFK"));
  bad_hits[1].setSequence(AASequence::fromString("MRSLGYVAVISAVATDTDK"));
  bad_hits[2].setSequence(AASequence::fromString("EGASTDFAALRTFLAEDGK"));
  bad_hits[3].setSequence(AASequence::fromString("DLEPGTDYEVTVSTLFGR"));
  bad_hits[4].setSequence(AASequence::fromString("FINFGVNVEVLSRFQTK"));
  bad_hits[5].setSequence(AASequence::fromString("MSLLSNMISIVKVGYNAR"));
  bad_hits[6].setSequence(AASequence::fromString("THPYGHAIVAGIERYPSK"));
  bad_hits[7].setSequence(AASequence::fromString("AITSDFANQAKTVLQNFK"));

  // modification-aware filtering:
  IDFilter::removePeptidesWithMatchingSequences(peptides_copy, bad_peptides,
                                                false);
  TEST_EQUAL(peptides_copy[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 3);
  TEST_EQUAL(peptide_hits[0].getSequence(),
             AASequence::fromString("TGCDTWGQGTLVTVSSASTK"));
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 10.93);
  TEST_EQUAL(peptide_hits[1].getSequence(),
             AASequence::fromString("TLCHHDATFDNLVWTPK"));
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 10.37);
  TEST_EQUAL(peptide_hits[2].getSequence(),
             AASequence::fromString("MSLLSNM(Oxidation)ISIVKVGYNAR"));
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 10);

  // modification-unaware filtering:
  IDFilter::removePeptidesWithMatchingSequences(peptides_copy, bad_peptides,
                                                true);
  TEST_EQUAL(peptides_copy[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 2);
  TEST_EQUAL(peptide_hits[0].getSequence(),
             AASequence::fromString("TGCDTWGQGTLVTVSSASTK"));
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 10.93);
  TEST_EQUAL(peptide_hits[1].getSequence(),
             AASequence::fromString("TLCHHDATFDNLVWTPK"));
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 10.37);
}
END_SECTION

START_SECTION((template<class PeakT> static void keepHitsMatchingProteins(MSExperiment<PeakT>& experiment, const vector<FASTAFile::FASTAEntry>& proteins)))
{
  MSExperiment<> experiment;
  vector<FASTAFile::FASTAEntry> proteins;
  vector<PeptideIdentification> peptides_copy = peptides;

  proteins.push_back(FASTAFile::FASTAEntry("Q824A5", "first desription",
                                           "LHASGITVTEIPVTATNFK"));
  proteins.push_back(FASTAFile::FASTAEntry("Q872T5", "second description",
                                           "THPYGHAIVAGIERYPSK"));

  for (Size i = 0; i < 5; ++i)
  {
    experiment.addSpectrum(MSSpectrum<>());
  }
  experiment[3].setMSLevel(2);
  experiment[3].setPeptideIdentifications(peptides_copy);

  IDFilter::keepHitsMatchingProteins(experiment, proteins);
  TEST_EQUAL(experiment[3].getPeptideIdentifications()[0].getScoreType(),
             "Mascot");

  vector<PeptideHit>& peptide_hits =
    experiment[3].getPeptideIdentifications()[0].getHits();
  TEST_EQUAL(peptide_hits.size(), 2);
  TEST_EQUAL(peptide_hits[0].getSequence(),
             AASequence::fromString("LHASGITVTEIPVTATNFK"));
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 34.85);
  TEST_EQUAL(peptide_hits[0].getRank(), 1);
  TEST_EQUAL(peptide_hits[1].getSequence(),
             AASequence::fromString("MRSLGYVAVISAVATDTDK"));
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 33.85);
  TEST_EQUAL(peptide_hits[1].getRank(), 2);
}
END_SECTION

START_SECTION((static void keepBestPeptideHits(vector<PeptideIdentification>& peptides, bool strict = false)))
{
  vector<PeptideIdentification> peptides_copy = peptides;
  vector<PeptideHit>& peptide_hits = peptides_copy[0].getHits();

  // not strict:
  IDFilter::keepBestPeptideHits(peptides_copy);
  TEST_EQUAL(peptides_copy[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 2);
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40);
  TEST_EQUAL(peptide_hits[0].getSequence(),
             AASequence::fromString("FINFGVNVEVLSRFQTK"));
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40);
  TEST_EQUAL(peptide_hits[1].getSequence(),
             AASequence::fromString("MSLLSNMISIVKVGYNAR"));

  // strict:
  IDFilter::keepBestPeptideHits(peptides_copy, true);
  TEST_EQUAL(peptides_copy[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 0);
}
END_SECTION

START_SECTION((template <class PeakT> static void filterHitsBySignificance(MSExperiment<PeakT>& experiment, double peptide_threshold_fraction, double protein_threshold_fraction)))
{
  MSExperiment<> experiment;
  vector<PeptideIdentification> ids(1, peptides[0]);

  ids[0].assignRanks();

  for (Size i = 0; i < 5; ++i)
  {
    experiment.addSpectrum(MSSpectrum<>());
  }
  experiment[3].setMSLevel(2);
  experiment[3].setPeptideIdentifications(ids);

  IDFilter::filterHitsBySignificance(experiment, 1.0, 1.0);
  PeptideIdentification& identification = experiment[3].getPeptideIdentifications()[0];
  TEST_EQUAL(identification.getScoreType(), "Mascot")

  vector<PeptideHit>& peptide_hits = identification.getHits();
  TEST_EQUAL(peptide_hits.size(), 5)
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40)
  TEST_EQUAL(peptide_hits[0].getRank(), 1)
  TEST_EQUAL((peptide_hits[0].getSequence() == AASequence::fromString("FINFGVNVEVLSRFQTK") && peptide_hits[1].getSequence() == AASequence::fromString("MSLLSNMISIVKVGYNAR")) || 
             (peptide_hits[0].getSequence() == AASequence::fromString("MSLLSNMISIVKVGYNAR") && peptide_hits[1].getSequence() == AASequence::fromString("FINFGVNVEVLSRFQTK")), true)

  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40)
  TEST_EQUAL(peptide_hits[1].getRank(), 1)
  TEST_EQUAL(peptide_hits[2].getSequence(), AASequence::fromString("THPYGHAIVAGIERYPSK"))
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 39)
  TEST_EQUAL(peptide_hits[2].getRank(), 2)
  TEST_EQUAL(peptide_hits[3].getSequence(), AASequence::fromString("LHASGITVTEIPVTATNFK"))
  TEST_REAL_SIMILAR(peptide_hits[3].getScore() , 34.85)
  TEST_EQUAL(peptide_hits[3].getRank(), 3)
  TEST_EQUAL(peptide_hits[4].getSequence(), AASequence::fromString("MRSLGYVAVISAVATDTDK"))
  TEST_REAL_SIMILAR(peptide_hits[4].getScore(), 33.85)
  TEST_EQUAL(peptide_hits[4].getRank(), 4)
}
END_SECTION

START_SECTION((template <class PeakT> static void filterHitsByScore(MSExperiment<PeakT>& experiment, double peptide_threshold_score, double protein_threshold_score)))
{
  MSExperiment<> experiment;
  vector<PeptideIdentification> ids(1, peptides[0]);

  ids[0].assignRanks();

  for (Size i = 0; i < 5; ++i)
  {
    experiment.addSpectrum(MSSpectrum<>());
  }
  experiment[3].setMSLevel(2);
  experiment[3].setPeptideIdentifications(ids);

  IDFilter::filterHitsByScore(experiment, 31.8621, 0);
  PeptideIdentification& identification = experiment[3].getPeptideIdentifications()[0];
  TEST_EQUAL(identification.getScoreType(), "Mascot")

  vector<PeptideHit>& peptide_hits = identification.getHits();
  TEST_EQUAL(peptide_hits.size(), 5)
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40)
  TEST_EQUAL(peptide_hits[0].getRank(), 1)
  TEST_EQUAL((peptide_hits[0].getSequence() == AASequence::fromString("FINFGVNVEVLSRFQTK") && peptide_hits[1].getSequence() == AASequence::fromString("MSLLSNMISIVKVGYNAR")) ||
             (peptide_hits[0].getSequence() == AASequence::fromString("MSLLSNMISIVKVGYNAR") && peptide_hits[1].getSequence() == AASequence::fromString("FINFGVNVEVLSRFQTK")), true)

  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40)
  TEST_EQUAL(peptide_hits[1].getRank(), 1)
  TEST_EQUAL(peptide_hits[2].getSequence(), AASequence::fromString("THPYGHAIVAGIERYPSK"))
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 39)
  TEST_EQUAL(peptide_hits[2].getRank(), 2)
  TEST_EQUAL(peptide_hits[3].getSequence(), AASequence::fromString("LHASGITVTEIPVTATNFK"))
  TEST_REAL_SIMILAR(peptide_hits[3].getScore(), 34.85)
  TEST_EQUAL(peptide_hits[3].getRank(), 3)
  TEST_EQUAL(peptide_hits[4].getSequence(), AASequence::fromString("MRSLGYVAVISAVATDTDK"))
  TEST_REAL_SIMILAR(peptide_hits[4].getScore(), 33.85)
  TEST_EQUAL(peptide_hits[4].getRank(), 4)
}
END_SECTION

START_SECTION((template <class PeakT> static void keepNBestHits(MSExperiment<PeakT>& experiment, Size n)))
{
  MSExperiment<> experiment;
  vector<PeptideIdentification> ids(1, peptides[0]);

  ids[0].assignRanks();

  for (Size i = 0; i < 5; ++i)
  {
    experiment.addSpectrum(MSSpectrum<>());
  }
  experiment[3].setMSLevel(2);
  experiment[3].setPeptideIdentifications(ids);

  IDFilter::keepNBestHits(experiment, 3);
  PeptideIdentification& identification = experiment[3].getPeptideIdentifications()[0];
  TEST_EQUAL(identification.getScoreType(), "Mascot")

  vector<PeptideHit>& peptide_hits = identification.getHits();
  TEST_EQUAL(peptide_hits.size(), 3)
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40)
  TEST_EQUAL(peptide_hits[0].getRank(), 1)
  TEST_EQUAL((peptide_hits[0].getSequence() == AASequence::fromString("FINFGVNVEVLSRFQTK") && peptide_hits[1].getSequence() == AASequence::fromString("MSLLSNMISIVKVGYNAR")) ||
             (peptide_hits[0].getSequence() == AASequence::fromString("MSLLSNMISIVKVGYNAR") && peptide_hits[1].getSequence() == AASequence::fromString("FINFGVNVEVLSRFQTK")), true)
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40)
  TEST_EQUAL(peptide_hits[1].getRank(), 1)
  TEST_EQUAL(peptide_hits[2].getSequence(), AASequence::fromString("THPYGHAIVAGIERYPSK"))
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 39)
  TEST_EQUAL(peptide_hits[2].getRank(), 2)
}
END_SECTION

START_SECTION((template <class IdentificationType> static void keepNBestHits(vector<IdentificationType>& ids, Size n)))
{
  vector<PeptideIdentification> peptides_copy = peptides;
  vector<PeptideHit>& peptide_hits = peptides_copy[0].getHits();

  IDFilter::keepNBestHits(peptides_copy, 3);
  TEST_EQUAL(peptides_copy[0].getScoreType(), "Mascot");

  TEST_EQUAL(peptide_hits.size(), 3);
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40);
  TEST_EQUAL(peptide_hits[0].getSequence(),
             AASequence::fromString("FINFGVNVEVLSRFQTK"));
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40);
  TEST_EQUAL(peptide_hits[1].getSequence(),
             AASequence::fromString("MSLLSNMISIVKVGYNAR"));
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 39);
  TEST_EQUAL(peptide_hits[2].getSequence(),
             AASequence::fromString("THPYGHAIVAGIERYPSK"));
}
END_SECTION

START_SECTION((static void filterPeptidesByRTPredictPValue(vector<PeptideIdentification>& peptides, const String& metavalue_key, double threshold = 0.05)))
{
  { // RT prediction:
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test2.idXML"),
                     proteins, peptides);
    IDFilter::filterPeptidesByRTPredictPValue(peptides, "predicted_RT_p_value",
                                              0.08);
    vector<PeptideHit>& hits = peptides[0].getHits();

    TEST_EQUAL(hits.size(), 4);
    TEST_EQUAL(hits[0].getSequence(),
               AASequence::fromString("LHASGITVTEIPVTATNFK"));
    TEST_EQUAL(hits[1].getSequence(),
               AASequence::fromString("DLEPGTDYEVTVSTLFGR"));
    TEST_EQUAL(hits[2].getSequence(),
               AASequence::fromString("FINFGVNVEVLSRFQTK"));
    TEST_EQUAL(hits[3].getSequence(),
               AASequence::fromString("MSLLSNMISIVKVGYNAR"));
  }
  { // first dim. RT prediction:
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test3.idXML"),
                     proteins, peptides);
    IDFilter::filterPeptidesByRTPredictPValue(peptides,
                                              "predicted_RT_p_value_first_dim",
                                              0.08);
    vector<PeptideHit>& hits = peptides[0].getHits();

    TEST_EQUAL(hits.size(), 4);
    TEST_EQUAL(hits[0].getSequence(),
               AASequence::fromString("LHASGITVTEIPVTATNFK"));
    TEST_EQUAL(hits[1].getSequence(),
               AASequence::fromString("DLEPGTDYEVTVSTLFGR"));
    TEST_EQUAL(hits[2].getSequence(),
               AASequence::fromString("FINFGVNVEVLSRFQTK"));
    TEST_EQUAL(hits[3].getSequence(),
               AASequence::fromString("MSLLSNMISIVKVGYNAR"));
  }
}
END_SECTION

START_SECTION((static void removeUnreferencedProteins(vector<ProteinIdentification>& proteins, vector<PeptideIdentification>& peptides)))
{
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test4.idXML"), proteins,
                   peptides);
  IDFilter::removeUnreferencedProteins(proteins, peptides);
  vector<ProteinHit>& hits = proteins[0].getHits();

  TEST_EQUAL(hits.size(), 3);
  TEST_EQUAL(hits[0].getAccession(), "Q824A5");
  TEST_EQUAL(hits[1].getAccession(), "S53854");
  TEST_EQUAL(hits[2].getAccession(), "Q872T5");
}
END_SECTION

// test for "updateProteinReferences" is missing!

START_SECTION((static void removeDuplicatePeptideHits(vector<PeptideIdentification>& peptides)))
{
  peptides.resize(1);
  vector<PeptideHit>& hits = peptides[0].getHits();
  hits.clear();
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("DFPIANGER"));
  hit.setCharge(1);
  hit.setScore(0.3);
  hits.push_back(hit);
  hit.setCharge(2);
  hits.push_back(hit);
  hit.setScore(0.5);
  hits.push_back(hit);
  hit.setSequence(AASequence::fromString("DFPIANGEK"));
  hits.push_back(hit);
  hits.push_back(hit);
  hits.push_back(hit);
  hit.setCharge(5);
  hits.push_back(hit);
  TEST_EQUAL(hits.size(), 7);

  IDFilter::removeDuplicatePeptideHits(peptides);
  TEST_EQUAL(hits.size(), 5);
  TEST_STRING_EQUAL(hits[3].getSequence().toString(), "DFPIANGEK");
  TEST_EQUAL(hits[3].getCharge(), 2);
  TEST_STRING_EQUAL(hits[4].getSequence().toString(), "DFPIANGEK");
  TEST_EQUAL(hits[4].getCharge(), 5);
}
END_SECTION

START_SECTION((bool updateProteinGroups(vector<ProteinIdentification::ProteinGroup>& groups, const vector<ProteinHit>& hits)))
{
  vector<ProteinIdentification::ProteinGroup> groups(2);
  groups[0].accessions.push_back("A");
  groups[0].probability = 0.1;
  groups[1].accessions.push_back("B");
  groups[1].accessions.push_back("C");
  groups[1].probability = 0.2;

  vector<ProteinHit> hits(3);
  hits[0].setAccession("C");
  hits[1].setAccession("B");
  hits[2].setAccession("A");

  vector<ProteinIdentification::ProteinGroup> groups_copy = groups;

  // no protein to remove:
  bool valid = IDFilter::updateProteinGroups(groups_copy, hits);
  TEST_EQUAL(valid, true);
  TEST_EQUAL(groups_copy.size(), 2);
  TEST_EQUAL(groups_copy == groups, true);
  
  // remove full protein group:
  hits.pop_back();
  valid = IDFilter::updateProteinGroups(groups_copy, hits);
  TEST_EQUAL(valid, true);
  TEST_EQUAL(groups_copy.size(), 1);
  TEST_EQUAL(groups_copy[0].accessions.size(), 2);
  TEST_EQUAL(groups_copy[0].accessions[0], "B");
  TEST_EQUAL(groups_copy[0].accessions[1], "C");
  TEST_EQUAL(groups_copy[0].probability, 0.2);
  
  // remove part of a protein group:
  hits.pop_back();
  valid = IDFilter::updateProteinGroups(groups_copy, hits);
  TEST_EQUAL(valid, false);
  TEST_EQUAL(groups_copy.size(), 1);
  TEST_EQUAL(groups_copy[0].accessions.size(), 1);
  TEST_EQUAL(groups_copy[0].accessions[0], "C");
  TEST_EQUAL(groups_copy[0].probability, 0.2);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

#pragma clang diagnostic pop

