// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/PercolatorOutfile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PercolatorOutfile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PercolatorOutfile* ptr = nullptr;
PercolatorOutfile* null_pointer = nullptr;
PercolatorOutfile file;

START_SECTION(PercolatorOutfile())
{
  ptr = new PercolatorOutfile();
  TEST_NOT_EQUAL(ptr, null_pointer);
}
END_SECTION

START_SECTION(~PercolatorOutfile())
{
  delete ptr;
}
END_SECTION


START_SECTION(enum ScoreType getScoreType(String score_type_name))
{
  TEST_EQUAL(PercolatorOutfile::getScoreType("qvalue"),
             PercolatorOutfile::QVALUE);
  TEST_EQUAL(PercolatorOutfile::getScoreType("q-value"),
             PercolatorOutfile::QVALUE);
  TEST_EQUAL(PercolatorOutfile::getScoreType("PEP"),
             PercolatorOutfile::POSTERRPROB);
  TEST_EQUAL(PercolatorOutfile::getScoreType("Posterior Error Probability"),
             PercolatorOutfile::POSTERRPROB);
  TEST_EQUAL(PercolatorOutfile::getScoreType("score"),
             PercolatorOutfile::SCORE);
}
END_SECTION

START_SECTION(void load(const String& filename, ProteinIdentification& proteins,
                        vector<PeptideIdentification>& peptides, 
                        SpectrumMetaDataLookup& lookup,
                        enum ScoreType output_score))
{
  // mock-up raw data like those used for the search:
  vector<MSSpectrum> spectra(3);
  double rt = 2.0;
  for (vector<MSSpectrum>::iterator it = spectra.begin();
       it != spectra.end(); ++it, rt += 1.0)
  {
    it->setMSLevel(2);
    it->setRT(rt);
    Precursor precursor;
    precursor.setCharge(Int(rt));
    precursor.setMZ(rt * 111.1);
    it->getPrecursors().push_back(precursor);
  }
  SpectrumMetaDataLookup lookup;
  lookup.readSpectra(spectra, ""); // no native IDs set, so don't parse them

  String filename = OPENMS_GET_TEST_DATA_PATH("PercolatorOutfile_test.psms");
  ProteinIdentification proteins;
  vector<PeptideIdentification> peptides;
  file.load(filename, proteins, peptides, lookup, PercolatorOutfile::SCORE);

  TEST_EQUAL(proteins.getHits().size(), 3);
  TEST_STRING_EQUAL(proteins.getHits()[0].getAccession(), "Protein1");
  TEST_STRING_EQUAL(proteins.getHits()[1].getAccession(), "Protein2");
  TEST_STRING_EQUAL(proteins.getHits()[2].getAccession(), "UniProt_P01834");

  TEST_EQUAL(proteins.getSearchParameters().fixed_modifications.size(), 1);
  TEST_STRING_EQUAL(proteins.getSearchParameters().fixed_modifications[0],
                    "Carbamidomethyl (C)");

  TEST_EQUAL(peptides.size(), 3);
  TEST_REAL_SIMILAR(peptides[0].getRT(), 2.0);
  TEST_REAL_SIMILAR(peptides[1].getRT(), 3.0);
  TEST_REAL_SIMILAR(peptides[2].getRT(), 4.0);
  TEST_REAL_SIMILAR(peptides[0].getMZ(), 222.2);
  TEST_REAL_SIMILAR(peptides[1].getMZ(), 333.3);
  TEST_REAL_SIMILAR(peptides[2].getMZ(), 444.4);
  TEST_EQUAL(peptides[0].getHits().size(), 1);
  TEST_EQUAL(peptides[1].getHits().size(), 1);
  TEST_EQUAL(peptides[2].getHits().size(), 1);
  TEST_EQUAL(peptides[0].getHits()[0].getCharge(), 2);
  TEST_EQUAL(peptides[1].getHits()[0].getCharge(), 3);
  TEST_EQUAL(peptides[2].getHits()[0].getCharge(), 4);
  TEST_REAL_SIMILAR(peptides[0].getHits()[0].getScore(), 6.77991);
  TEST_REAL_SIMILAR(peptides[1].getHits()[0].getScore(), 6.57945);
  TEST_REAL_SIMILAR(peptides[2].getHits()[0].getScore(), 6.50586);
  TEST_STRING_EQUAL(peptides[0].getHits()[0].getSequence().toString(),
                    "VDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSK");
  TEST_STRING_EQUAL(peptides[1].getHits()[0].getSequence().toString(),
                    "VDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSK");
  TEST_STRING_EQUAL(peptides[2].getHits()[0].getSequence().toString(),
                    "VTLSC(Carbamidomethyl)TGSSSNLGAGYDVHWYQQLPGTAPK");
  TEST_REAL_SIMILAR(peptides[0].getHits()[0].getMetaValue("Percolator_score"),
                    6.77991);
  TEST_REAL_SIMILAR(peptides[1].getHits()[0].getMetaValue("Percolator_qvalue"),
                    0.0);
  TEST_REAL_SIMILAR(peptides[2].getHits()[0].getMetaValue("Percolator_PEP"),
                    1.8014e-14);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
