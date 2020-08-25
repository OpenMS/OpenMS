// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
