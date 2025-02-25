// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>
//#include <OpenMS/KERNEL/MSSpectrum.h>
//#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>

using namespace OpenMS;

START_TEST(OPXLDataStructs, "$Id$")

  OPXLDataStructs::ProteinProteinCrossLink cross_link;
  AASequence alpha = AASequence::fromString("PEPTIDE");
  AASequence beta = AASequence::fromString("EDEPITPEPE");
  cross_link.alpha = &alpha;
  cross_link.beta = &beta;
  cross_link.cross_link_position = std::make_pair<SignedSize, SignedSize>(3, 5);
  cross_link.cross_linker_mass = 150.0;
  cross_link.cross_linker_name = "NOTDSS";
  cross_link.term_spec_alpha = ResidueModification::N_TERM;
  cross_link.term_spec_beta = ResidueModification::ANYWHERE;

START_SECTION(ProteinProteinCrossLink())

  TEST_EQUAL(cross_link.getType(), OPXLDataStructs::CROSS)

  cross_link.beta = nullptr;
  TEST_EQUAL(cross_link.getType(), OPXLDataStructs::LOOP)

  cross_link.cross_link_position = std::make_pair<SignedSize, SignedSize>(3, -1);
  TEST_EQUAL(cross_link.getType(), OPXLDataStructs::MONO)

END_SECTION


START_SECTION(XLPrecursor())

  std::vector<OPXLDataStructs::XLPrecursor> precursors;
  for (Size i = 20; i > 1; --i)
  {
    OPXLDataStructs::XLPrecursor prec;
    prec.precursor_mass = i * 3.33;
    prec.alpha_index = 1;
    prec.beta_index = 2;
    precursors.push_back(prec);
  }

  // sorting using the XLPrecursorComparator
  std::sort(precursors.begin(), precursors.end(), OPXLDataStructs::XLPrecursorComparator());

  for (Size i = 0; i < precursors.size()-1; ++i)
  {
    TEST_EQUAL(precursors[i].precursor_mass < precursors[i+1].precursor_mass, true)
  }

  // searching for a precursor mass using a double value
  std::vector< OPXLDataStructs::XLPrecursor >::const_iterator low_it;
  low_it = lower_bound(precursors.begin(), precursors.end(), 9 * 3.33 - 1, OPXLDataStructs::XLPrecursorComparator());
  TEST_REAL_SIMILAR((*low_it).precursor_mass, 9 * 3.33)

END_SECTION

START_SECTION(AASeqWithMass())

  std::vector<OPXLDataStructs::AASeqWithMass> peptides;

  OPXLDataStructs::AASeqWithMass pep;
  pep.position = OPXLDataStructs::INTERNAL;

  pep.peptide_seq = AASequence::fromString("TESTEE");
  pep.peptide_mass = pep.peptide_seq.getMonoWeight();
  peptides.push_back(pep);

  pep.peptide_seq = AASequence::fromString("TESTEEE");
  pep.peptide_mass = pep.peptide_seq.getMonoWeight();
  peptides.push_back(pep);

  pep.peptide_seq = AASequence::fromString("TESTEEEEEEEEEEEE");
  pep.peptide_mass = pep.peptide_seq.getMonoWeight();
  peptides.push_back(pep);

  pep.peptide_seq = AASequence::fromString("TESTEEEEE");
  pep.peptide_mass = pep.peptide_seq.getMonoWeight();
  peptides.push_back(pep);

  pep.peptide_seq = AASequence::fromString("TES");
  pep.peptide_mass = pep.peptide_seq.getMonoWeight();
  peptides.push_back(pep);

  // sorting using the AASeqWithMassComparator
  std::sort(peptides.begin(), peptides.end(), OPXLDataStructs::AASeqWithMassComparator());

  for (Size i = 0; i < peptides.size()-1; ++i)
  {
    TEST_EQUAL(peptides[i].peptide_mass < peptides[i+1].peptide_mass, true)
  }

  // searching for a peptide mass using a double value
  std::vector< OPXLDataStructs::AASeqWithMass >::const_iterator low_it;
  low_it = lower_bound(peptides.begin(), peptides.end(), AASequence::fromString("TESTEEE").getMonoWeight() - 0.1, OPXLDataStructs::AASeqWithMassComparator());
  TEST_REAL_SIMILAR((*low_it).peptide_mass, AASequence::fromString("TESTEEE").getMonoWeight())

END_SECTION

END_TEST
