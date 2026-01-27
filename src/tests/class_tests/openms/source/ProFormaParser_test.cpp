// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/ProForma.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <fstream>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

// Type aliases for ProForma nested types
using Peptidoform = ProForma::Peptidoform;
using PeptidoformIon = ProForma::PeptidoformIon;
using SequenceElement = ProForma::SequenceElement;
using SequenceSection = ProForma::SequenceSection;
using Modification = ProForma::Modification;
using ModificationTag = ProForma::ModificationTag;
using CvAccession = ProForma::CvAccession;
using CvDatabase = ProForma::CvDatabase;
using NamedMod = ProForma::NamedMod;
using MassDelta = ProForma::MassDelta;
using FormulaTag = ProForma::FormulaTag;
using GlycanComposition = ProForma::GlycanComposition;
using InfoTag = ProForma::InfoTag;
using PositionConstraint = ProForma::PositionConstraint;
using Label = ProForma::Label;
using UnlocalisedMod = ProForma::UnlocalisedMod;
using LabileModification = ProForma::LabileModification;
using GlobalModEntry = ProForma::GlobalModEntry;
using GlobalModification = ProForma::GlobalModification;
using IsotopeReplacement = ProForma::IsotopeReplacement;
using AmbiguousRegion = ProForma::AmbiguousRegion;
using ModifiedRange = ProForma::ModifiedRange;
using ChargeState = ProForma::ChargeState;
using AdductIon = ProForma::AdductIon;
using ConversionIssue = ProForma::ConversionIssue;
using ConversionIssueType = ProForma::ConversionIssueType;
using ConversionPolicy = ProForma::ConversionPolicy;
using WriteMode = ProForma::WriteMode;
using ErrorCode = ProForma::ErrorCode;

///////////////////////////

// Helper function to load test cases from fixture file
vector<string> loadTestCases(const string& filename)
{
  vector<string> cases;
  ifstream file(filename);
  if (!file.is_open())
  {
    throw std::runtime_error("Failed to open test fixture file: " + filename);
  }
  string line;
  while (getline(file, line))
  {
    // Skip empty lines and comments
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    cases.push_back(line);
  }
  return cases;
}

///////////////////////////

START_TEST(ProFormaParser, "$Id$")

/////////////////////////////////////////////////////////////
// Parser tests (Tokenizer is now internal to ProForma.cpp)
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - simple sequences)
{
  // Simple sequence without modifications
  Peptidoform pf1 = ProForma::parse("PEPTIDE");
  TEST_EQUAL(pf1.sequence.size(), 7)

  // Single amino acid
  Peptidoform pf2 = ProForma::parse("A");
  TEST_EQUAL(pf2.sequence.size(), 1)

  // Two amino acids
  Peptidoform pf3 = ProForma::parse("AA");
  TEST_EQUAL(pf3.sequence.size(), 2)
}
END_SECTION

START_SECTION(ProForma::parse - CV accessions)
{
  // UNIMOD accession
  Peptidoform pf1 = ProForma::parse("EM[UNIMOD:35]K");
  TEST_EQUAL(pf1.sequence.size(), 3)

  // Get the modification on M (second residue)
  auto* seq_elem_ptr = std::get_if<SequenceElement>(&pf1.sequence[1]);
  TEST_NOT_EQUAL(seq_elem_ptr, nullptr)
  if (seq_elem_ptr)
  {
    TEST_EQUAL(seq_elem_ptr->amino_acid, 'M')
    TEST_EQUAL(seq_elem_ptr->modifications.size(), 1)

    if (!seq_elem_ptr->modifications.empty())
    {
      auto& mod = seq_elem_ptr->modifications[0];
      TEST_EQUAL(mod.alternatives.size(), 1)

      if (!mod.alternatives.empty())
      {
        // Check it's a CvAccession
        auto* cv = std::get_if<CvAccession>(&mod.alternatives[0].first);
        TEST_NOT_EQUAL(cv, nullptr)
        if (cv)
        {
          TEST_EQUAL(cv->database, CvDatabase::UNIMOD)
          TEST_EQUAL(cv->accession, "35")
        }
      }
    }
  }
}
END_SECTION

START_SECTION(ProForma::parse - named modifications)
{
  Peptidoform pf = ProForma::parse("EM[Oxidation]K");
  TEST_EQUAL(pf.sequence.size(), 3)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[1]);
  auto& mod = seq_elem.modifications[0];

  auto* named = std::get_if<NamedMod>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(named, nullptr)
  TEST_EQUAL(named->name, "Oxidation")
}
END_SECTION

START_SECTION(ProForma::parse - mass deltas)
{
  Peptidoform pf = ProForma::parse("A[+15.9949]A");
  TEST_EQUAL(pf.sequence.size(), 2)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[0]);
  auto& mod = seq_elem.modifications[0];

  auto* mass = std::get_if<MassDelta>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(mass, nullptr)
  TEST_REAL_SIMILAR(mass->mass, 15.9949)
  TEST_EQUAL(mass->original_text, "+15.9949")
}
END_SECTION

START_SECTION(ProForma::parse - formula tags)
{
  Peptidoform pf = ProForma::parse("SEQUEN[Formula:C12H20O2]CE");
  TEST_EQUAL(pf.sequence.size(), 8)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[5]); // N
  auto& mod = seq_elem.modifications[0];

  auto* formula = std::get_if<FormulaTag>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(formula, nullptr)
  TEST_EQUAL(formula->formula_string, "C12H20O2")
}
END_SECTION

START_SECTION(ProForma::parse - glycan compositions)
{
  Peptidoform pf = ProForma::parse("SEQUEN[Glycan:HexNAc]CE");
  TEST_EQUAL(pf.sequence.size(), 8)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[5]); // N
  auto& mod = seq_elem.modifications[0];

  auto* glycan = std::get_if<GlycanComposition>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(glycan, nullptr)
  TEST_EQUAL(glycan->components.size(), 1)
}
END_SECTION

START_SECTION(ProForma::parse - N-terminal modifications)
{
  Peptidoform pf = ProForma::parse("[Acetyl]-PEPTIDE");
  TEST_EQUAL(pf.sequence.size(), 7)
  TEST_EQUAL(pf.n_term_mods.size(), 1)

  auto* named = std::get_if<NamedMod>(&pf.n_term_mods[0].alternatives[0].first);
  TEST_NOT_EQUAL(named, nullptr)
  TEST_EQUAL(named->name, "Acetyl")
}
END_SECTION

START_SECTION(ProForma::parse - C-terminal modifications)
{
  Peptidoform pf = ProForma::parse("PEPTIDE-[Amidated]");
  TEST_EQUAL(pf.sequence.size(), 7)
  TEST_EQUAL(pf.c_term_mods.size(), 1)

  auto* named = std::get_if<NamedMod>(&pf.c_term_mods[0].alternatives[0].first);
  TEST_NOT_EQUAL(named, nullptr)
  TEST_EQUAL(named->name, "Amidated")
}
END_SECTION

START_SECTION(ProForma::parse - unlocalised modifications)
{
  Peptidoform pf = ProForma::parse("[Phospho]?PEPTIDE");
  TEST_EQUAL(pf.sequence.size(), 7)
  TEST_EQUAL(pf.unlocalised_mods.size(), 1)
}
END_SECTION

START_SECTION(ProForma::parse - labile modifications)
{
  Peptidoform pf = ProForma::parse("{Glycan:Hex}PEPTIDE");
  TEST_EQUAL(pf.sequence.size(), 7)
  TEST_EQUAL(pf.labile_mods.size(), 1)
}
END_SECTION

START_SECTION(ProForma::parse - global modifications)
{
  Peptidoform pf = ProForma::parse("<13C>PEPTIDE");
  TEST_EQUAL(pf.sequence.size(), 7)
  TEST_EQUAL(pf.global_mods.size(), 1)

  auto* isotope = std::get_if<IsotopeReplacement>(&pf.global_mods[0]);
  TEST_NOT_EQUAL(isotope, nullptr)
  TEST_EQUAL(isotope->isotope, "13C")
}
END_SECTION

START_SECTION(ProForma::parse - cross-link labels)
{
  // EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1] has 12 amino acids
  Peptidoform pf = ProForma::parse("EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]");
  TEST_EQUAL(pf.sequence.size(), 12)

  // Check first cross-link site at K (index 5)
  auto& elem1 = std::get<SequenceElement>(pf.sequence[5]); // K
  TEST_EQUAL(elem1.amino_acid, 'K')
  TEST_EQUAL(elem1.modifications.size(), 1)
  auto& label1 = elem1.modifications[0].alternatives[0].second;
  TEST_EQUAL(label1.has_value(), true)
  TEST_EQUAL(label1->identifier, "XL1")

  // Check second cross-link site at K (index 11)
  auto& elem2 = std::get<SequenceElement>(pf.sequence[11]); // K
  TEST_EQUAL(elem2.amino_acid, 'K')
  TEST_EQUAL(elem2.modifications.size(), 1)

  // Verify roundtrip - label-only [#XL1] should be preserved
  String roundtrip = ProForma::toString(pf);
  TEST_EQUAL(roundtrip, "EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]")
}
END_SECTION

START_SECTION(ProForma::parse - modification alternatives)
{
  Peptidoform pf = ProForma::parse("ELVIS[Phospho|+79.966331]K");
  TEST_EQUAL(pf.sequence.size(), 6)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[4]); // S
  auto& mod = seq_elem.modifications[0];
  TEST_EQUAL(mod.alternatives.size(), 2)

  // First alternative: NamedMod
  auto* named = std::get_if<NamedMod>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(named, nullptr)
  TEST_EQUAL(named->name, "Phospho")

  // Second alternative: MassDelta
  auto* mass = std::get_if<MassDelta>(&mod.alternatives[1].first);
  TEST_NOT_EQUAL(mass, nullptr)
}
END_SECTION

START_SECTION(ProForma::parse - ambiguous regions)
{
  Peptidoform pf = ProForma::parse("(?DQ)PEPTIDE");
  TEST_EQUAL(pf.sequence.size(), 8) // ambiguous region counts as 1

  auto* ambig = std::get_if<AmbiguousRegion>(&pf.sequence[0]);
  TEST_NOT_EQUAL(ambig, nullptr)
  TEST_EQUAL(ambig->elements.size(), 2)
}
END_SECTION

START_SECTION(ProForma::parse - modified ranges)
{
  Peptidoform pf = ProForma::parse("PROT(EOSFORMS)[+19.0523]ISK");
  TEST_EQUAL(pf.sequence.size(), 8) // (EOSFORMS) counts as 1

  auto* range = std::get_if<ModifiedRange>(&pf.sequence[4]);
  TEST_NOT_EQUAL(range, nullptr)
  TEST_EQUAL(range->elements.size(), 8)
  TEST_EQUAL(range->modifications.size(), 1)
}
END_SECTION

START_SECTION(ProForma::parseIon - charge states)
{
  PeptidoformIon ion1 = ProForma::parseIon("PEPTIDE/2");
  TEST_EQUAL(ion1.chains.size(), 1)
  TEST_EQUAL(ion1.charge.has_value(), true)

  auto* simple_charge = std::get_if<int>(&ion1.charge.value());
  TEST_NOT_EQUAL(simple_charge, nullptr)
  TEST_EQUAL(*simple_charge, 2)
}
END_SECTION

START_SECTION(ProForma::parseIon - multiple chains)
{
  PeptidoformIon ion = ProForma::parseIon("PEPTIDE//SEQUENCE");
  TEST_EQUAL(ion.chains.size(), 2)
  TEST_EQUAL(ion.chains[0].sequence.size(), 7)
  TEST_EQUAL(ion.chains[1].sequence.size(), 8)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Additional CV database tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - MOD database accessions)
{
  Peptidoform pf = ProForma::parse("EM[MOD:00719]EVEES[MOD:00046]PEK");
  TEST_EQUAL(pf.sequence.size(), 10)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[1]); // M
  auto& mod = seq_elem.modifications[0];
  auto* cv = std::get_if<CvAccession>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(cv, nullptr)
  TEST_EQUAL(cv->database, CvDatabase::MOD)
  TEST_EQUAL(cv->accession, "00719")
}
END_SECTION

START_SECTION(ProForma::parse - RESID database accessions)
{
  Peptidoform pf = ProForma::parse("EM[RESID:AA0581]EVEES[RESID:AA0037]PEK");
  TEST_EQUAL(pf.sequence.size(), 10)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[1]); // M
  auto& mod = seq_elem.modifications[0];
  auto* cv = std::get_if<CvAccession>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(cv, nullptr)
  TEST_EQUAL(cv->database, CvDatabase::RESID)
  TEST_EQUAL(cv->accession, "AA0581")
}
END_SECTION

START_SECTION(ProForma::parse - GNO database accessions)
{
  Peptidoform pf = ProForma::parse("NEEYN[GNO:G59626AS]K");
  TEST_EQUAL(pf.sequence.size(), 6)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[4]); // N
  auto& mod = seq_elem.modifications[0];
  auto* cv = std::get_if<CvAccession>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(cv, nullptr)
  TEST_EQUAL(cv->database, CvDatabase::GNO)
  TEST_EQUAL(cv->accession, "G59626AS")
}
END_SECTION

START_SECTION(ProForma::parse - XLMOD database accessions)
{
  Peptidoform pf = ProForma::parse("EMEVTK[XLMOD:02001]SESPEK");
  TEST_EQUAL(pf.sequence.size(), 12)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[5]); // K
  auto& mod = seq_elem.modifications[0];
  auto* cv = std::get_if<CvAccession>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(cv, nullptr)
  TEST_EQUAL(cv->database, CvDatabase::XLMOD)
  TEST_EQUAL(cv->accession, "02001")
}
END_SECTION

/////////////////////////////////////////////////////////////
// Isotope formula tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - isotope formulas with brackets)
{
  Peptidoform pf = ProForma::parse("SEQUEN[Formula:[13C2][12C-2]H2N]CE");
  TEST_EQUAL(pf.sequence.size(), 8)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[5]); // N
  auto& mod = seq_elem.modifications[0];
  auto* formula = std::get_if<FormulaTag>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(formula, nullptr)
  TEST_EQUAL(formula->formula_string, "[13C2][12C-2]H2N")
}
END_SECTION

START_SECTION(ProForma::parse - multiple isotope labels)
{
  Peptidoform pf = ProForma::parse("<13C><15N>ATPEILTVNSIGQLK");
  TEST_EQUAL(pf.sequence.size(), 15)
  TEST_EQUAL(pf.global_mods.size(), 2)

  auto* iso1 = std::get_if<IsotopeReplacement>(&pf.global_mods[0]);
  TEST_NOT_EQUAL(iso1, nullptr)
  TEST_EQUAL(iso1->isotope, "13C")

  auto* iso2 = std::get_if<IsotopeReplacement>(&pf.global_mods[1]);
  TEST_NOT_EQUAL(iso2, nullptr)
  TEST_EQUAL(iso2->isotope, "15N")
}
END_SECTION

START_SECTION(ProForma::parse - deuterium label)
{
  Peptidoform pf = ProForma::parse("<D>ATPEILTVNSIGQLK");
  TEST_EQUAL(pf.sequence.size(), 15)
  TEST_EQUAL(pf.global_mods.size(), 1)

  auto* iso = std::get_if<IsotopeReplacement>(&pf.global_mods[0]);
  TEST_NOT_EQUAL(iso, nullptr)
  TEST_EQUAL(iso->isotope, "D")
}
END_SECTION

/////////////////////////////////////////////////////////////
// Localization score tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - localization scores)
{
  // E-M-E-V-T-S-E-S-P-E-K = 11 amino acids
  Peptidoform pf = ProForma::parse("EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK");
  TEST_EQUAL(pf.sequence.size(), 11)

  // Check score on T at position 4
  auto& elem_t = std::get<SequenceElement>(pf.sequence[4]); // T
  TEST_EQUAL(elem_t.modifications.size(), 1)
  auto& label_t = elem_t.modifications[0].alternatives[0].second;
  TEST_EQUAL(label_t.has_value(), true)
  TEST_EQUAL(label_t->identifier, "g1")
  TEST_EQUAL(label_t->score.has_value(), true)
  TEST_REAL_SIMILAR(label_t->score.value(), 0.01)

  // Check score on S at position 5
  auto& elem_s1 = std::get<SequenceElement>(pf.sequence[5]); // S
  auto& label_s1 = elem_s1.modifications[0].alternatives[0].second;
  TEST_EQUAL(label_s1.has_value(), true)
  TEST_REAL_SIMILAR(label_s1->score.value(), 0.09)

  // Check score on S at position 7 (has the Phospho mod)
  auto& elem_s2 = std::get<SequenceElement>(pf.sequence[7]); // S
  auto& label_s2 = elem_s2.modifications[0].alternatives[0].second;
  TEST_EQUAL(label_s2.has_value(), true)
  TEST_REAL_SIMILAR(label_s2->score.value(), 0.90)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Multiplicity tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - multiplicity with caret)
{
  // E-M-E-V-T-S-E-S-P-E-K = 11 amino acids
  Peptidoform pf = ProForma::parse("[Phospho]^2?[Acetyl]-EM[Oxidation]EVTSESPEK");
  TEST_EQUAL(pf.sequence.size(), 11)
  TEST_EQUAL(pf.unlocalised_mods.size(), 1)
  TEST_EQUAL(pf.n_term_mods.size(), 1)

  // Check occurrence (multiplicity) on the unlocalised Phospho
  auto& unloc = pf.unlocalised_mods[0];
  TEST_EQUAL(unloc.occurrence.has_value(), true)
  TEST_EQUAL(unloc.occurrence.value(), 2)
}
END_SECTION

START_SECTION(ProForma::parse - multiple modifications on range)
{
  Peptidoform pf = ProForma::parse("MPGLVDSNPAPPESQEKKPLK(PCCACPETKKARDACIIEKGEEHCGHLIEAHKECMRALGFKI)[Oxidation][Oxidation][half cystine][half cystine]");
  TEST_EQUAL(pf.sequence.size(), 22) // 21 AA + 1 modified range

  auto* range = std::get_if<ModifiedRange>(&pf.sequence[21]);
  TEST_NOT_EQUAL(range, nullptr)
  TEST_EQUAL(range->modifications.size(), 4)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Chimeric spectra tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parseIon - chimeric spectra with plus)
{
  PeptidoformIon ion = ProForma::parseIon("EMEVEESPEK+ELVISLIVER");
  TEST_EQUAL(ion.chains.size(), 2)
  TEST_EQUAL(ion.is_chimeric, true)
  TEST_EQUAL(ion.chains[0].sequence.size(), 10)
  TEST_EQUAL(ion.chains[1].sequence.size(), 10)

  // Verify roundtrip
  String output = ProForma::toString(ion);
  TEST_EQUAL(output, "EMEVEESPEK+ELVISLIVER")
}
END_SECTION

START_SECTION(ProForma::parseIon - chimeric spectra with charges)
{
  PeptidoformIon ion = ProForma::parseIon("EMEVEESPEK/2+ELVISLIVER/3");
  TEST_EQUAL(ion.chains.size(), 2)
  TEST_EQUAL(ion.is_chimeric, true)

  // Check per-chain charges
  TEST_EQUAL(ion.chains[0].charge.has_value(), true)
  auto* charge0 = std::get_if<int>(&ion.chains[0].charge.value());
  TEST_NOT_EQUAL(charge0, nullptr)
  TEST_EQUAL(*charge0, 2)

  TEST_EQUAL(ion.chains[1].charge.has_value(), true)
  auto* charge1 = std::get_if<int>(&ion.chains[1].charge.value());
  TEST_NOT_EQUAL(charge1, nullptr)
  TEST_EQUAL(*charge1, 3)

  // Verify roundtrip
  String output = ProForma::toString(ion);
  TEST_EQUAL(output, "EMEVEESPEK/2+ELVISLIVER/3")
}
END_SECTION

/////////////////////////////////////////////////////////////
// Adduct/charge notation tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parseIon - adduct charge notation)
{
  PeptidoformIon ion = ProForma::parseIon("PEPTIDE/[Na:z+1]");
  TEST_EQUAL(ion.chains.size(), 1)
  TEST_EQUAL(ion.charge.has_value(), true)

  // ChargeState is variant<int, vector<AdductIon>>
  auto* adducts = std::get_if<std::vector<AdductIon>>(&ion.charge.value());
  TEST_NOT_EQUAL(adducts, nullptr)
  TEST_EQUAL(adducts->size(), 1)
  TEST_EQUAL(adducts->at(0).formula, "Na")
  TEST_EQUAL(adducts->at(0).charge, 1)
}
END_SECTION

START_SECTION(ProForma::parseIon - multiple adducts)
{
  PeptidoformIon ion = ProForma::parseIon("PEPTIDE/[Na:z+1,H:z+1]");
  TEST_EQUAL(ion.chains.size(), 1)

  auto* adducts = std::get_if<std::vector<AdductIon>>(&ion.charge.value());
  TEST_NOT_EQUAL(adducts, nullptr)
  TEST_EQUAL(adducts->size(), 2)
}
END_SECTION

START_SECTION(ProForma::parseIon - adduct with count)
{
  PeptidoformIon ion = ProForma::parseIon("PEPTIDE/[Na:z+1^2]");
  TEST_EQUAL(ion.chains.size(), 1)

  auto* adducts = std::get_if<std::vector<AdductIon>>(&ion.charge.value());
  TEST_NOT_EQUAL(adducts, nullptr)
  TEST_EQUAL(adducts->at(0).occurrence.has_value(), true)
  TEST_EQUAL(adducts->at(0).occurrence.value(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Gene/protein prefix tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parseIon - gene prefix)
{
  PeptidoformIon ion = ProForma::parseIon("(>Trypsin)AANSIPYQVSLNS+(>Keratin)AKEQFERQTA");
  TEST_EQUAL(ion.chains.size(), 2)
  TEST_EQUAL(ion.is_chimeric, true)

  // Gene/protein prefix is stored in the 'name' field of Peptidoform
  TEST_EQUAL(ion.chains[0].name.has_value(), true)
  TEST_EQUAL(ion.chains[0].name.value(), "Trypsin")

  TEST_EQUAL(ion.chains[1].name.has_value(), true)
  TEST_EQUAL(ion.chains[1].name.value(), "Keratin")

  // Verify roundtrip
  String output = ProForma::toString(ion);
  TEST_EQUAL(output, "(>Trypsin)AANSIPYQVSLNS+(>Keratin)AKEQFERQTA")
}
END_SECTION

START_SECTION(ProForma::parse - gene prefix single chain)
{
  Peptidoform pf = ProForma::parse("(>sp|P12345|PROT_HUMAN)PEPTIDE");
  TEST_EQUAL(pf.name.has_value(), true)
  TEST_EQUAL(pf.name.value(), "sp|P12345|PROT_HUMAN")
  TEST_EQUAL(pf.sequence.size(), 7)

  // Verify roundtrip
  String output = ProForma::toString(pf);
  TEST_EQUAL(output, "(>sp|P12345|PROT_HUMAN)PEPTIDE")
}
END_SECTION

/////////////////////////////////////////////////////////////
// Cation modification tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - cation modifications)
{
  // Cation:Mg[II] is parsed as a named modification
  Peptidoform pf = ProForma::parse("EM[Oxidation]EVE[Cation:Mg[II]]ES[Phospho]PEK");
  TEST_EQUAL(pf.sequence.size(), 10)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[4]); // E
  TEST_EQUAL(seq_elem.modifications.size(), 1)
  // The modification is parsed - specific type depends on implementation
}
END_SECTION

START_SECTION(ProForma::parse - aluminum cation)
{
  // Cation:Al[III] is parsed as a named modification
  Peptidoform pf = ProForma::parse("PE[Cation:Al[III]]PTIDE");
  TEST_EQUAL(pf.sequence.size(), 7)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[1]); // E
  TEST_EQUAL(seq_elem.modifications.size(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Half cystine tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - half cystine)
{
  Peptidoform pf = ProForma::parse("EVTSEKC[half cystine]LEMSC[half cystine]EFD");
  TEST_EQUAL(pf.sequence.size(), 15)

  auto& elem1 = std::get<SequenceElement>(pf.sequence[6]); // C
  auto& mod1 = elem1.modifications[0];
  auto* named1 = std::get_if<NamedMod>(&mod1.alternatives[0].first);
  TEST_NOT_EQUAL(named1, nullptr)
  TEST_EQUAL(named1->name, "half cystine")

  auto& elem2 = std::get<SequenceElement>(pf.sequence[11]); // C
  auto& mod2 = elem2.modifications[0];
  auto* named2 = std::get_if<NamedMod>(&mod2.alternatives[0].first);
  TEST_NOT_EQUAL(named2, nullptr)
  TEST_EQUAL(named2->name, "half cystine")
}
END_SECTION

/////////////////////////////////////////////////////////////
// BRANCH cross-link tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parseIon - branch cross-links)
{
  PeptidoformIon ion = ProForma::parseIon("ETFGD[MOD:00093#BRANCH]//R[#BRANCH]ATER");
  TEST_EQUAL(ion.chains.size(), 2)

  // Check first chain has BRANCH label
  auto& elem1 = std::get<SequenceElement>(ion.chains[0].sequence[4]); // D
  auto& label1 = elem1.modifications[0].alternatives[0].second;
  TEST_EQUAL(label1.has_value(), true)
  TEST_EQUAL(label1->identifier, "BRANCH")

  // Check second chain references BRANCH
  auto& elem2 = std::get<SequenceElement>(ion.chains[1].sequence[0]); // R
  auto& label2 = elem2.modifications[0].alternatives[0].second;
  TEST_EQUAL(label2.has_value(), true)
  TEST_EQUAL(label2->identifier, "BRANCH")
}
END_SECTION

START_SECTION(ProForma::parseIon - branch with C-terminal)
{
  PeptidoformIon ion = ProForma::parseIon("AVTKYTSSK[MOD:00134#BRANCH]//AGKQLEDGRTLSDYNIQKESTLHLVLRLRG-[#BRANCH]");
  TEST_EQUAL(ion.chains.size(), 2)

  // Second chain should have BRANCH on C-terminus
  TEST_EQUAL(ion.chains[1].c_term_mods.size(), 1)
  auto& label = ion.chains[1].c_term_mods[0].alternatives[0].second;
  TEST_EQUAL(label.has_value(), true)
  TEST_EQUAL(label->identifier, "BRANCH")
}
END_SECTION

/////////////////////////////////////////////////////////////
// INFO tag tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - INFO tags)
{
  Peptidoform pf = ProForma::parse("ELV[INFO:AnyString]IS");
  TEST_EQUAL(pf.sequence.size(), 5)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[2]); // V
  auto& mod = seq_elem.modifications[0];
  auto* info = std::get_if<InfoTag>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(info, nullptr)
  TEST_EQUAL(info->text, "AnyString")
}
END_SECTION

START_SECTION(ProForma::parse - INFO tags case insensitive)
{
  Peptidoform pf = ProForma::parse("ELV[info:AnyString]IS");
  TEST_EQUAL(pf.sequence.size(), 5)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[2]); // V
  auto& mod = seq_elem.modifications[0];
  auto* info = std::get_if<InfoTag>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(info, nullptr)
  TEST_EQUAL(info->text, "AnyString")
}
END_SECTION

START_SECTION(ProForma::parse - modification with INFO alternative)
{
  Peptidoform pf = ProForma::parse("ELVIS[Phospho|INFO:newly discovered]K");
  TEST_EQUAL(pf.sequence.size(), 6)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[4]); // S
  auto& mod = seq_elem.modifications[0];
  TEST_EQUAL(mod.alternatives.size(), 2)

  auto* named = std::get_if<NamedMod>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(named, nullptr)
  TEST_EQUAL(named->name, "Phospho")

  auto* info = std::get_if<InfoTag>(&mod.alternatives[1].first);
  TEST_NOT_EQUAL(info, nullptr)
  TEST_EQUAL(info->text, "newly discovered")
}
END_SECTION

START_SECTION(ProForma::parse - multiple INFO tags)
{
  Peptidoform pf = ProForma::parse("ELVIS[Phospho|INFO:newly discovered|INFO:really awesome]K");
  TEST_EQUAL(pf.sequence.size(), 6)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[4]); // S
  auto& mod = seq_elem.modifications[0];
  TEST_EQUAL(mod.alternatives.size(), 3)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Position constraint tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - position constraints on modifications)
{
  // Position constraints are part of modification alternatives
  Peptidoform pf = ProForma::parse("PEPTI(MERMERMERM)[Oxidation|Position:M][Oxidation|Position:M]DE");
  TEST_EQUAL(pf.sequence.size(), 8) // PEPTI (5) + (range) (1) + DE (2)

  auto* range = std::get_if<ModifiedRange>(&pf.sequence[5]);
  TEST_NOT_EQUAL(range, nullptr)
  TEST_EQUAL(range->modifications.size(), 2)

  // Check that each modification has two alternatives: Oxidation and Position:M
  for (const auto& mod : range->modifications)
  {
    TEST_EQUAL(mod.alternatives.size(), 2)

    // First alternative should be Oxidation (NamedMod)
    auto* named = std::get_if<NamedMod>(&mod.alternatives[0].first);
    TEST_NOT_EQUAL(named, nullptr)
    TEST_EQUAL(named->name, "Oxidation")

    // Second alternative should be Position:M (PositionConstraint)
    auto* pos = std::get_if<PositionConstraint>(&mod.alternatives[1].first);
    TEST_NOT_EQUAL(pos, nullptr)
    TEST_EQUAL(pos->residues.size(), 1)
    TEST_EQUAL(pos->residues[0], 'M')
  }

  // Test roundtrip
  String output = ProForma::toString(pf);
  TEST_EQUAL(output, "PEPTI(MERMERMERM)[Oxidation|Position:M][Oxidation|Position:M]DE")
}
END_SECTION

START_SECTION(ProForma::parse - position constraints with multiple residues)
{
  // Position constraint with multiple allowed residues
  Peptidoform pf = ProForma::parse("PEPTIDE[Phospho|Position:STY]");
  TEST_EQUAL(pf.sequence.size(), 7)

  // Check modification on the last E
  auto* elem = std::get_if<SequenceElement>(&pf.sequence[6]);
  TEST_NOT_EQUAL(elem, nullptr)
  TEST_EQUAL(elem->modifications.size(), 1)
  TEST_EQUAL(elem->modifications[0].alternatives.size(), 2)

  // Second alternative should be Position:STY
  auto* pos = std::get_if<PositionConstraint>(&elem->modifications[0].alternatives[1].first);
  TEST_NOT_EQUAL(pos, nullptr)
  TEST_EQUAL(pos->residues.size(), 3)
  TEST_EQUAL(pos->residues[0], 'S')
  TEST_EQUAL(pos->residues[1], 'T')
  TEST_EQUAL(pos->residues[2], 'Y')

  // Test roundtrip
  String output = ProForma::toString(pf);
  TEST_EQUAL(output, "PEPTIDE[Phospho|Position:STY]")
}
END_SECTION

START_SECTION(ProForma::parse - position constraints on unlocalised)
{
  // Position constraints can be specified on unlocalised mods
  Peptidoform pf = ProForma::parse("[Oxidation|CoMKP]?PEPT[Phospho]IDE");
  TEST_EQUAL(pf.sequence.size(), 7)
  TEST_EQUAL(pf.unlocalised_mods.size(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Multiple terminal modifications tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - multiple N-terminal modifications)
{
  Peptidoform pf = ProForma::parse("[Acetyl][Carbamyl]-QPEPTIDE");
  TEST_EQUAL(pf.sequence.size(), 8)
  TEST_EQUAL(pf.n_term_mods.size(), 2)

  auto* named1 = std::get_if<NamedMod>(&pf.n_term_mods[0].alternatives[0].first);
  TEST_NOT_EQUAL(named1, nullptr)
  TEST_EQUAL(named1->name, "Acetyl")

  auto* named2 = std::get_if<NamedMod>(&pf.n_term_mods[1].alternatives[0].first);
  TEST_NOT_EQUAL(named2, nullptr)
  TEST_EQUAL(named2->name, "Carbamyl")
}
END_SECTION

START_SECTION(ProForma::parse - multiple C-terminal modifications)
{
  Peptidoform pf = ProForma::parse("PEPTIDEG-[Methyl][Amidated]");
  TEST_EQUAL(pf.sequence.size(), 8)
  TEST_EQUAL(pf.c_term_mods.size(), 2)

  auto* named1 = std::get_if<NamedMod>(&pf.c_term_mods[0].alternatives[0].first);
  TEST_NOT_EQUAL(named1, nullptr)
  TEST_EQUAL(named1->name, "Methyl")

  auto* named2 = std::get_if<NamedMod>(&pf.c_term_mods[1].alternatives[0].first);
  TEST_NOT_EQUAL(named2, nullptr)
  TEST_EQUAL(named2->name, "Amidated")
}
END_SECTION

/////////////////////////////////////////////////////////////
// Observed mass prefix tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - observed mass prefix)
{
  Peptidoform pf = ProForma::parse("EM[U:+15.995]EVEES[Obs:+79.978]PEK");
  TEST_EQUAL(pf.sequence.size(), 10)

  // Check Obs prefix on S - MassDelta uses 'source' enum
  auto& seq_elem = std::get<SequenceElement>(pf.sequence[6]); // S
  auto& mod = seq_elem.modifications[0];
  auto* mass = std::get_if<MassDelta>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(mass, nullptr)
  TEST_EQUAL(mass->source, MassDelta::Source::OBS)
  TEST_REAL_SIMILAR(mass->mass, 79.978)
}
END_SECTION

START_SECTION(ProForma::parse - UNIMOD prefix on mass)
{
  Peptidoform pf = ProForma::parse("EM[U:+15.9949]EVEES[U:+79.9663]PEK");
  TEST_EQUAL(pf.sequence.size(), 10)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[1]); // M
  auto& mod = seq_elem.modifications[0];
  auto* mass = std::get_if<MassDelta>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(mass, nullptr)
  TEST_EQUAL(mass->source, MassDelta::Source::U)
}
END_SECTION

/////////////////////////////////////////////////////////////
// C-term targeted global modifications tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - C-term targeted global mod)
{
  Peptidoform pf = ProForma::parse("<[Oxidation]@W,C-term:G>QATPEILTWCNSIGCLKG");
  TEST_EQUAL(pf.sequence.size(), 18)
  TEST_EQUAL(pf.global_mods.size(), 1)

  // GlobalModEntry is variant<IsotopeReplacement, GlobalModification>
  auto* global_mod = std::get_if<GlobalModification>(&pf.global_mods[0]);
  TEST_NOT_EQUAL(global_mod, nullptr)
  // Check targets include W and C-term:G
  TEST_EQUAL(global_mod->locations.size() >= 2, true)
}
END_SECTION

START_SECTION(ProForma::parse - N-term targeted global mod)
{
  Peptidoform pf = ProForma::parse("<[TMT6plex]@K,N-term>ATPEILTCNSIGCLK");
  TEST_EQUAL(pf.sequence.size(), 15)
  TEST_EQUAL(pf.global_mods.size(), 1)

  auto* global_mod = std::get_if<GlobalModification>(&pf.global_mods[0]);
  TEST_NOT_EQUAL(global_mod, nullptr)
}
END_SECTION

START_SECTION(ProForma::parse - N-term with specific residue)
{
  Peptidoform pf = ProForma::parse("<[TMT6plex]@K,N-term:A>ATPEILTCNSIGCLK");
  TEST_EQUAL(pf.sequence.size(), 15)
  TEST_EQUAL(pf.global_mods.size(), 1)

  auto* global_mod = std::get_if<GlobalModification>(&pf.global_mods[0]);
  TEST_NOT_EQUAL(global_mod, nullptr)
}
END_SECTION

START_SECTION(ProForma::parse - amidated C-term global)
{
  Peptidoform pf = ProForma::parse("<[Amidated]@C-term>QATPEILTWCNSIGCLKG");
  TEST_EQUAL(pf.sequence.size(), 18)
  TEST_EQUAL(pf.global_mods.size(), 1)

  auto* global_mod = std::get_if<GlobalModification>(&pf.global_mods[0]);
  TEST_NOT_EQUAL(global_mod, nullptr)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Formula with charge tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - formula with charge)
{
  Peptidoform pf = ProForma::parse("SEQUEN[Formula:Zn1:z+2]CE");
  TEST_EQUAL(pf.sequence.size(), 8)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[5]); // N
  auto& mod = seq_elem.modifications[0];
  auto* formula = std::get_if<FormulaTag>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(formula, nullptr)
  TEST_EQUAL(formula->formula_string, "Zn1")
  TEST_EQUAL(formula->charge.has_value(), true)
  TEST_EQUAL(formula->charge.value(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Ion fragment tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - a-type ion)
{
  Peptidoform pf = ProForma::parse("PEPTID-[a-type-ion]");
  TEST_EQUAL(pf.sequence.size(), 6)
  TEST_EQUAL(pf.c_term_mods.size(), 1)

  auto* named = std::get_if<NamedMod>(&pf.c_term_mods[0].alternatives[0].first);
  TEST_NOT_EQUAL(named, nullptr)
  TEST_EQUAL(named->name, "a-type-ion")
}
END_SECTION

START_SECTION(ProForma::parse - d-ion with formula)
{
  Peptidoform pf = ProForma::parse("PEPTID[Formula:H-1C-1O-2|Info:d-ion]-[a-type-ion]");
  TEST_EQUAL(pf.sequence.size(), 6)

  // Check the D residue has the formula modification
  auto& seq_elem = std::get<SequenceElement>(pf.sequence[5]); // D
  TEST_EQUAL(seq_elem.modifications.size(), 1)
  TEST_EQUAL(seq_elem.modifications[0].alternatives.size(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Cross-link with multiple chains and XL tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parseIon - complex cross-links with chimeric separator)
{
  // Test complex cross-links with + (chimeric) separator combining multiple cross-linked complexes
  // Complex 1: A[X:DSS#XL1]//B[#XL1] (two chains cross-linked via XL1)
  // Complex 2: C[X:DSS#XL2]//D[#XL2] (two chains cross-linked via XL2)
  PeptidoformIon ion = ProForma::parseIon("A[X:DSS#XL1]//B[#XL1]+C[X:DSS#XL2]//D[#XL2]");
  TEST_EQUAL(ion.chains.size(), 4)
  TEST_EQUAL(ion.is_chimeric, true)

  // Check chain 0 (A with DSS modification and XL1 label)
  auto& chain0_elem = std::get<SequenceElement>(ion.chains[0].sequence[0]);
  TEST_EQUAL(chain0_elem.amino_acid, 'A')
  TEST_EQUAL(chain0_elem.modifications.size(), 1)
  auto& mod0_label = chain0_elem.modifications[0].alternatives[0].second;
  TEST_EQUAL(mod0_label.has_value(), true)
  TEST_EQUAL(mod0_label->identifier, "XL1")

  // Check chain 1 (B with XL1 label only)
  auto& chain1_elem = std::get<SequenceElement>(ion.chains[1].sequence[0]);
  TEST_EQUAL(chain1_elem.amino_acid, 'B')
  auto& mod1_label = chain1_elem.modifications[0].alternatives[0].second;
  TEST_EQUAL(mod1_label.has_value(), true)
  TEST_EQUAL(mod1_label->identifier, "XL1")

  // Check chain 2 (C with DSS modification and XL2 label)
  auto& chain2_elem = std::get<SequenceElement>(ion.chains[2].sequence[0]);
  TEST_EQUAL(chain2_elem.amino_acid, 'C')
  auto& mod2_label = chain2_elem.modifications[0].alternatives[0].second;
  TEST_EQUAL(mod2_label.has_value(), true)
  TEST_EQUAL(mod2_label->identifier, "XL2")

  // Check chain 3 (D with XL2 label only)
  auto& chain3_elem = std::get<SequenceElement>(ion.chains[3].sequence[0]);
  TEST_EQUAL(chain3_elem.amino_acid, 'D')
  auto& mod3_label = chain3_elem.modifications[0].alternatives[0].second;
  TEST_EQUAL(mod3_label.has_value(), true)
  TEST_EQUAL(mod3_label->identifier, "XL2")
}
END_SECTION

START_SECTION(ProForma::parseIon - disulfide cross-links)
{
  PeptidoformIon ion = ProForma::parseIon("EVTSEKC[XLMOD:02009#XL1]LEMSC[#XL1]EFD");
  TEST_EQUAL(ion.chains.size(), 1)

  auto& elem1 = std::get<SequenceElement>(ion.chains[0].sequence[6]); // C
  auto& label1 = elem1.modifications[0].alternatives[0].second;
  TEST_EQUAL(label1.has_value(), true)
  TEST_EQUAL(label1->identifier, "XL1")
}
END_SECTION

/////////////////////////////////////////////////////////////
// Glycan composition tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - complex glycan composition)
{
  Peptidoform pf = ProForma::parse("SEQUEN[Glycan:HexNAc1Hex2]CE");
  TEST_EQUAL(pf.sequence.size(), 8)

  auto& seq_elem = std::get<SequenceElement>(pf.sequence[5]); // N
  auto& mod = seq_elem.modifications[0];
  auto* glycan = std::get_if<GlycanComposition>(&mod.alternatives[0].first);
  TEST_NOT_EQUAL(glycan, nullptr)
  TEST_EQUAL(glycan->components.size(), 2)
}
END_SECTION

START_SECTION(ProForma::parse - multiple labile glycans)
{
  Peptidoform pf = ProForma::parse("{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK");
  TEST_EQUAL(pf.sequence.size(), 10)
  TEST_EQUAL(pf.labile_mods.size(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Unlocalised modification score tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::parse - unlocalised with scores)
{
  // E-M-E-V-T-S-E-S-P-E-K = 11 amino acids
  Peptidoform pf = ProForma::parse("[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]S[#s1(0.09)]ES[#s1(0.90)]PEK");
  TEST_EQUAL(pf.sequence.size(), 11)
  TEST_EQUAL(pf.unlocalised_mods.size(), 1)

  // UnlocalisedMod has 'modifications' field (vector of Modification)
  auto& unloc = pf.unlocalised_mods[0];
  TEST_EQUAL(unloc.modifications.size() >= 1, true)

  // Check that the first modification has a label with identifier "s1"
  auto& mod = unloc.modifications[0];
  TEST_EQUAL(mod.alternatives.size() >= 1, true)
  auto& label = mod.alternatives[0].second;
  TEST_EQUAL(label.has_value(), true)
  TEST_EQUAL(label->identifier, "s1")
}
END_SECTION

/////////////////////////////////////////////////////////////
// Writer tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::toString - simple sequences)
{
  Peptidoform pf = ProForma::parse("PEPTIDE");
  String result = ProForma::toString(pf);
  TEST_EQUAL(result, "PEPTIDE")
}
END_SECTION

START_SECTION(ProForma::toString - modifications)
{
  Peptidoform pf = ProForma::parse("EM[UNIMOD:35]K");
  String result = ProForma::toString(pf);
  TEST_EQUAL(result, "EM[UNIMOD:35]K")
}
END_SECTION

START_SECTION(ProForma::toString - terminal modifications)
{
  Peptidoform pf = ProForma::parse("[Acetyl]-PEPTIDE-[Amidated]");
  String result = ProForma::toString(pf);
  TEST_EQUAL(result, "[Acetyl]-PEPTIDE-[Amidated]")
}
END_SECTION

/////////////////////////////////////////////////////////////
// Roundtrip tests - parse then write should produce equivalent output
/////////////////////////////////////////////////////////////

START_SECTION(Roundtrip tests from positive fixture)
{
  // Test a selection of cases from the positive test fixtures
  vector<string> test_cases = {
    "AA",
    "A[+1]",
    "EM[Oxidation]EVEES[Phospho]PEK",
    "EM[UNIMOD:35]K",
    "[Acetyl]-PEPTIDE",
    "PEPTIDE-[Amidated]",
    "[Phospho]?PEPTIDE",
    "{Glycan:Hex}PEPTIDE",
    "<13C>PEPTIDE",
    "SEQUEN[Formula:C12H20O2]CE",
    "SEQUEN[Glycan:HexNAc]CE"
  };

  for (const auto& test : test_cases)
  {
    Peptidoform pf = ProForma::parse(test);
    String result = ProForma::toString(pf);
    // Re-parse and compare structure
    Peptidoform pf2 = ProForma::parse(result);
    TEST_EQUAL(pf.sequence.size(), pf2.sequence.size())
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// Error handling tests
/////////////////////////////////////////////////////////////

START_SECTION(ProFormaParser error handling - unclosed bracket)
{
  TEST_EXCEPTION(ProForma::ParseError, ProForma::parse("A[+1"))
}
END_SECTION

START_SECTION(ProFormaParser error handling - empty sequence)
{
  TEST_EXCEPTION(ProForma::ParseError, ProForma::parse(""))
}
END_SECTION

START_SECTION(ProFormaParser error handling - invalid global mod)
{
  TEST_EXCEPTION(ProForma::ParseError, ProForma::parse("<[TMT6plex]>AA"))
}
END_SECTION

START_SECTION(ProForma::ParseError - position clamping for noexcept safety)
{
  // Test 1: Position way beyond input.size() should be clamped
  {
    String input = "ABC";
    size_t out_of_bounds_position = 100;

    ProForma::ParseError err(
      __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      ErrorCode::UNEXPECTED_CHARACTER,
      out_of_bounds_position,
      input,
      "Test error"
    );

    TEST_EQUAL(err.getPosition(), input.size())
    String formatted = err.getFormattedMessage();
    TEST_EQUAL(formatted.empty(), false)
    // Should show end of input marker since position is at end
    TEST_EQUAL(formatted.hasSubstring("END OF INPUT"), true)
  }

  // Test 2: Position exactly at input.size() (boundary case)
  {
    String input = "ABCDEF";
    size_t boundary_position = input.size(); // Position 6, exactly at end

    ProForma::ParseError err(
      __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      ErrorCode::UNEXPECTED_END_OF_INPUT,
      boundary_position,
      input,
      "Unexpected end"
    );

    TEST_EQUAL(err.getPosition(), input.size())
    TEST_EQUAL(err.getContextBefore(), "ABCDEF")
    TEST_EQUAL(err.getContextAfter(), "")
  }

  // Test 3: Empty input with position 0
  {
    String input = "";
    size_t position = 0;

    ProForma::ParseError err(
      __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      ErrorCode::EMPTY_SEQUENCE,
      position,
      input,
      "Empty sequence"
    );

    TEST_EQUAL(err.getPosition(), 0)
    TEST_EQUAL(err.getContextBefore(), "")
    TEST_EQUAL(err.getContextAfter(), "")
  }

  // Test 4: Different error code (verify error code is preserved)
  {
    String input = "PEP[";
    size_t position = 3;

    ProForma::ParseError err(
      __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      ErrorCode::UNCLOSED_BRACKET,
      position,
      input,
      "Unclosed bracket"
    );

    TEST_EQUAL(err.getErrorCode(), ErrorCode::UNCLOSED_BRACKET)
    String formatted = err.getFormattedMessage();
    TEST_EQUAL(formatted.hasSubstring("Unclosed bracket"), true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// JSON serialization tests
/////////////////////////////////////////////////////////////

START_SECTION(JSON serialization - Peptidoform)
{
  Peptidoform pf = ProForma::parse("EM[UNIMOD:35]K");

  // Serialize to JSON
  String json = ProForma::toJSON(pf);
  TEST_EQUAL(json.empty(), false)

  // Deserialize back
  Peptidoform pf2 = ProForma::peptidoformFromJSON(json);
  TEST_EQUAL(pf2.sequence.size(), pf.sequence.size())
}
END_SECTION

START_SECTION(JSON serialization - PeptidoformIon)
{
  PeptidoformIon ion = ProForma::parseIon("PEPTIDE//SEQUENCE/2");

  // Serialize to JSON
  String json = ProForma::toJSON(ion);
  TEST_EQUAL(json.empty(), false)

  // Deserialize back
  PeptidoformIon ion2 = ProForma::peptidoformIonFromJSON(json);
  TEST_EQUAL(ion2.chains.size(), ion.chains.size())
}
END_SECTION

/////////////////////////////////////////////////////////////
// WriteMode tests (lossless vs canonical)
/////////////////////////////////////////////////////////////

START_SECTION(WriteMode - lossless preserves original mass format)
{
  // Parse a mass delta with specific formatting
  Peptidoform pf = ProForma::parse("EM[+15.99]K");

  // Lossless mode should preserve the original text
  String lossless = ProForma::toString(pf, WriteMode::LOSSLESS);
  TEST_EQUAL(lossless, "EM[+15.99]K")
}
END_SECTION

START_SECTION(WriteMode - canonical uses fixed precision)
{
  // Parse a mass delta with specific formatting
  Peptidoform pf = ProForma::parse("EM[+15.99]K");

  // Canonical mode should use fixed 4 decimal places
  String canonical = ProForma::toString(pf, WriteMode::CANONICAL);
  TEST_EQUAL(canonical, "EM[+15.9900]K")
}
END_SECTION

START_SECTION(WriteMode - canonical normalizes mass precision)
{
  // Parse with many decimal places
  Peptidoform pf = ProForma::parse("EM[+15.99491234]K");

  // Canonical mode should normalize to 4 decimal places
  String canonical = ProForma::toString(pf, WriteMode::CANONICAL);
  TEST_EQUAL(canonical, "EM[+15.9949]K")
}
END_SECTION

START_SECTION(WriteMode - CV accessions in both modes)
{
  // CV accessions should be the same in both modes
  Peptidoform pf = ProForma::parse("EM[UNIMOD:35]K");

  String lossless = ProForma::toString(pf, WriteMode::LOSSLESS);
  String canonical = ProForma::toString(pf, WriteMode::CANONICAL);

  // Both should produce the same output for CV accessions
  TEST_EQUAL(lossless, "EM[UNIMOD:35]K")
  TEST_EQUAL(canonical, "EM[UNIMOD:35]K")
}
END_SECTION

START_SECTION(WriteMode - named modifications in both modes)
{
  // Named modifications should be the same in both modes
  Peptidoform pf = ProForma::parse("EM[Oxidation]K");

  String lossless = ProForma::toString(pf, WriteMode::LOSSLESS);
  String canonical = ProForma::toString(pf, WriteMode::CANONICAL);

  // Both should produce the same output for named mods
  TEST_EQUAL(lossless, "EM[Oxidation]K")
  TEST_EQUAL(canonical, "EM[Oxidation]K")
}
END_SECTION

START_SECTION(WriteMode - default is LOSSLESS)
{
  // Default mode should be lossless
  Peptidoform pf = ProForma::parse("EM[+15.99]K");

  String default_output = ProForma::toString(pf);
  String lossless = ProForma::toString(pf, WriteMode::LOSSLESS);

  TEST_EQUAL(default_output, lossless)
}
END_SECTION

/////////////////////////////////////////////////////////////
// AASequence Conversion Tests
/////////////////////////////////////////////////////////////

START_SECTION(resolveModifications - UNIMOD lookup)
{
  // Parse a ProForma string with UNIMOD modification
  Peptidoform pf = ProForma::parse("EM[UNIMOD:35]K");

  // Resolve modifications
  ProForma::resolveModifications(pf);

  // Check that the modification was resolved
  TEST_EQUAL(pf.sequence.size(), 3)

  // The second element (M) should have a resolved modification
  const SequenceElement* m_elem = std::get_if<SequenceElement>(&pf.sequence[1]);
  TEST_NOT_EQUAL(m_elem, nullptr)
  TEST_EQUAL(m_elem->modifications.size(), 1)

  // Check that resolved_mod is not null (UNIMOD:35 = Oxidation)
  TEST_NOT_EQUAL(m_elem->modifications[0].resolved_mod, nullptr)
}
END_SECTION

START_SECTION(resolveModifications - named modification lookup)
{
  // Parse a ProForma string with named modification
  Peptidoform pf = ProForma::parse("EM[Oxidation]K");

  // Resolve modifications
  ProForma::resolveModifications(pf);

  // The second element (M) should have a resolved modification
  const SequenceElement* m_elem = std::get_if<SequenceElement>(&pf.sequence[1]);
  TEST_NOT_EQUAL(m_elem, nullptr)
  TEST_EQUAL(m_elem->modifications.size(), 1)

  // Oxidation should be resolved
  TEST_NOT_EQUAL(m_elem->modifications[0].resolved_mod, nullptr)
}
END_SECTION

START_SECTION(isRepresentableAsAASequence - simple sequence)
{
  // Simple unmodified sequence should be representable
  Peptidoform pf = ProForma::parse("PEPTIDE");
  TEST_EQUAL(ProForma::isRepresentableAsAASequence(pf), true)

  // Sequence with UNIMOD modification should be representable
  Peptidoform pf2 = ProForma::parse("EM[UNIMOD:35]K");
  TEST_EQUAL(ProForma::isRepresentableAsAASequence(pf2), true)
}
END_SECTION

START_SECTION(isRepresentableAsAASequence - unsupported features)
{
  // Unlocalised modification is not representable
  Peptidoform pf_unloc = ProForma::parse("[Phospho]?PEPTIDE");
  TEST_EQUAL(ProForma::isRepresentableAsAASequence(pf_unloc), false)

  // Labile modification is not representable
  Peptidoform pf_labile = ProForma::parse("{Glycan:Hex}PEPTIDE");
  TEST_EQUAL(ProForma::isRepresentableAsAASequence(pf_labile), false)

  // Ambiguous region is not representable
  Peptidoform pf_ambig = ProForma::parse("PEP(?DQ)IDE");
  TEST_EQUAL(ProForma::isRepresentableAsAASequence(pf_ambig), false)
}
END_SECTION

START_SECTION(getAASequenceConversionIssues)
{
  // Test that conversion issues are properly reported

  // Unlocalised modification
  Peptidoform pf_unloc = ProForma::parse("[Phospho]?PEPTIDE");
  auto issues_unloc = ProForma::getAASequenceConversionIssues(pf_unloc);
  TEST_EQUAL(issues_unloc.empty(), false)
  bool found_unloc_issue = false;
  for (const auto& issue : issues_unloc)
  {
    if (issue.type == ConversionIssueType::UNLOCALISED_MOD)
    {
      found_unloc_issue = true;
      break;
    }
  }
  TEST_EQUAL(found_unloc_issue, true)

  // Labile modification
  Peptidoform pf_labile = ProForma::parse("{Glycan:Hex}PEPTIDE");
  auto issues_labile = ProForma::getAASequenceConversionIssues(pf_labile);
  TEST_EQUAL(issues_labile.empty(), false)
  bool found_labile_issue = false;
  for (const auto& issue : issues_labile)
  {
    if (issue.type == ConversionIssueType::LABILE_MOD)
    {
      found_labile_issue = true;
      break;
    }
  }
  TEST_EQUAL(found_labile_issue, true)
}
END_SECTION

START_SECTION(toAASequence - simple sequence)
{
  // Convert simple unmodified sequence
  Peptidoform pf = ProForma::parse("PEPTIDE");
  AASequence seq = ProForma::toAASequence(pf);

  TEST_EQUAL(seq.toUnmodifiedString(), "PEPTIDE")
  TEST_EQUAL(seq.size(), 7)
}
END_SECTION

START_SECTION(toAASequence - with UNIMOD modification)
{
  // Convert sequence with UNIMOD modification (Oxidation on M)
  Peptidoform pf = ProForma::parse("EM[UNIMOD:35]K");
  AASequence seq = ProForma::toAASequence(pf);

  TEST_EQUAL(seq.toUnmodifiedString(), "EMK")
  TEST_EQUAL(seq.size(), 3)

  // Check that M is modified
  TEST_EQUAL(seq[1].isModified(), true)
}
END_SECTION

START_SECTION(toAASequence - with N-terminal modification)
{
  // Convert sequence with N-terminal acetylation
  Peptidoform pf = ProForma::parse("[UNIMOD:1]-PEPTIDE");
  AASequence seq = ProForma::toAASequence(pf);

  TEST_EQUAL(seq.toUnmodifiedString(), "PEPTIDE")
  TEST_EQUAL(seq.hasNTerminalModification(), true)
}
END_SECTION

START_SECTION(toAASequence - FAIL_ON_LOSS policy throws on unsupported)
{
  // FAIL_ON_LOSS policy should throw on unlocalised modifications
  Peptidoform pf = ProForma::parse("[Phospho]?PEPTIDE");

  TEST_EXCEPTION(Exception::ConversionError,
    ProForma::toAASequence(pf, ConversionPolicy::FAIL_ON_LOSS))
}
END_SECTION

START_SECTION(toAASequence - BEST_EFFORT policy ignores unsupported)
{
  // BEST_EFFORT policy should not throw
  Peptidoform pf = ProForma::parse("[Phospho]?PEPTIDE");

  AASequence seq = ProForma::toAASequence(pf, ConversionPolicy::BEST_EFFORT);
  TEST_EQUAL(seq.toUnmodifiedString(), "PEPTIDE")
}
END_SECTION

START_SECTION(fromAASequence - simple sequence)
{
  // Create AASequence and convert to Peptidoform
  AASequence seq = AASequence::fromString("PEPTIDE");
  Peptidoform pf = ProForma::fromAASequence(seq);

  TEST_EQUAL(pf.sequence.size(), 7)

  // Check amino acids
  for (size_t i = 0; i < 7; ++i)
  {
    const SequenceElement* elem = std::get_if<SequenceElement>(&pf.sequence[i]);
    TEST_NOT_EQUAL(elem, nullptr)
  }
}
END_SECTION

START_SECTION(fromAASequence - with modification)
{
  // Create AASequence with oxidation
  AASequence seq = AASequence::fromString("EM(Oxidation)K");
  Peptidoform pf = ProForma::fromAASequence(seq);

  TEST_EQUAL(pf.sequence.size(), 3)

  // Check that M has a modification
  const SequenceElement* m_elem = std::get_if<SequenceElement>(&pf.sequence[1]);
  TEST_NOT_EQUAL(m_elem, nullptr)
  TEST_EQUAL(m_elem->amino_acid, 'M')
  TEST_EQUAL(m_elem->modifications.empty(), false)
}
END_SECTION

START_SECTION(fromAASequence - roundtrip)
{
  // Test roundtrip: AASequence -> Peptidoform -> toString -> parse -> toAASequence
  AASequence orig = AASequence::fromString("EM(Oxidation)K");

  // Convert to Peptidoform
  Peptidoform pf = ProForma::fromAASequence(orig);

  // Convert to string
  String proforma_str = ProForma::toString(pf);

  // Parse back
  Peptidoform pf2 = ProForma::parse(proforma_str);

  // Convert back to AASequence
  AASequence result = ProForma::toAASequence(pf2);

  // Compare
  TEST_EQUAL(result.toUnmodifiedString(), orig.toUnmodifiedString())
  TEST_EQUAL(result.size(), orig.size())
}
END_SECTION

START_SECTION(ProForma to AASequence roundtrip)
{
  // Test roundtrip: ProForma -> Peptidoform -> AASequence -> Peptidoform -> ProForma
  String orig = "EM[UNIMOD:35]K";

  // Parse
  Peptidoform pf = ProForma::parse(orig);

  // Convert to AASequence
  AASequence seq = ProForma::toAASequence(pf);

  // Convert back to Peptidoform
  Peptidoform pf2 = ProForma::fromAASequence(seq);

  // The result should be similar (may use different notation for the same mod)
  TEST_EQUAL(pf2.sequence.size(), pf.sequence.size())
}
END_SECTION

/////////////////////////////////////////////////////////////
// Load and test from fixture files
/////////////////////////////////////////////////////////////

START_SECTION(Positive test cases from fixture file)
{
  string fixture_path = OPENMS_GET_TEST_DATA_PATH("ProFormaParser_positive_tests.txt");
  vector<string> positive_tests = loadTestCases(fixture_path);

  // Helper to detect if a test case needs parseIon() instead of parse()
  // Cases with '//', '+', or '/[digit]' need parseIon()
  auto needsParseIon = [](const std::string& s) -> bool
  {
    // Check for multi-chain separator
    if (s.find("//") != std::string::npos) return true;
    // Check for chimeric spectra separator (but not inside brackets)
    // Simple heuristic: '+' not preceded by '[' context
    size_t pos = 0;
    int bracket_depth = 0;
    while (pos < s.size())
    {
      if (s[pos] == '[') bracket_depth++;
      else if (s[pos] == ']') bracket_depth--;
      else if (s[pos] == '+' && bracket_depth == 0 && pos > 0)
      {
        // Check it's not a mass delta like [+1.5]
        if (pos > 0 && s[pos-1] != '[' && s[pos-1] != ':' && s[pos-1] != '|')
        {
          return true;
        }
      }
      else if (s[pos] == '/' && bracket_depth == 0 && pos + 1 < s.size())
      {
        // Check for charge state: /[digit] or /+ or /- or /[
        char next = s[pos + 1];
        if (std::isdigit(next) || next == '+' || next == '-' || next == '[')
        {
          return true;
        }
      }
      pos++;
    }
    return false;
  };

  // Test all cases (file must exist, loadTestCases throws if it doesn't)
  if (!positive_tests.empty())  // Handle case where file only contains comments
  {
    int passed = 0;
    int failed = 0;

    for (const auto& test_case : positive_tests)
    {
      try
      {
        if (needsParseIon(test_case))
        {
          PeptidoformIon ion = ProForma::parseIon(test_case);
        }
        else
        {
          Peptidoform pf = ProForma::parse(test_case);
        }
        passed++;
      }
      catch (const ProForma::ParseError& e)
      {
        // Some complex cases may not be fully implemented yet
        failed++;
        // Uncomment for debugging:
        // std::cerr << "FAIL (should pass): " << test_case << " - " << e.what() << std::endl;
      }
      catch (const std::exception& e)
      {
        failed++;
        // Uncomment for debugging:
        // std::cerr << "FAIL (should pass): " << test_case << " - " << e.what() << std::endl;
      }
    }

    // All positive test cases must parse successfully
    TEST_EQUAL(failed, 0)
  }
}
END_SECTION

START_SECTION(Negative test cases from fixture file)
{
  string fixture_path = OPENMS_GET_TEST_DATA_PATH("ProFormaParser_negative_tests.txt");
  vector<string> negative_tests = loadTestCases(fixture_path);

  // Test that parsing fails for all cases (file must exist, loadTestCases throws if it doesn't)
  if (!negative_tests.empty())  // Handle case where file only contains comments
  {
    int correctly_rejected = 0;
    int incorrectly_accepted = 0;

    for (const auto& test_case : negative_tests)
    {
      try
      {
        Peptidoform pf = ProForma::parse(test_case);
        // If we get here, parsing succeeded when it shouldn't have
        incorrectly_accepted++;
        // Uncomment for debugging:
        // std::cerr << "FAIL (should reject): " << test_case << std::endl;
      }
      catch (const ProForma::ParseError&)
      {
        // Good - parsing should fail for negative cases
        correctly_rejected++;
      }
      catch (const std::exception&)
      {
        // Other exceptions also count as rejection
        correctly_rejected++;
      }
    }

    // All negative test cases must be rejected
    TEST_EQUAL(incorrectly_accepted, 0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// Component-level grammar tests
// Tests individual ProForma grammar components extracted from
// the component test fixture files
/////////////////////////////////////////////////////////////

START_SECTION(Component-level grammar tests)
{
  int passed = 0;
  int failed = 0;

  // Helper lambda to test that a string parses successfully as peptidoform
  auto expectPass = [&](const std::string& input, const std::string& context = "")
  {
    try
    {
      ProForma::parse(input);
      passed++;
    }
    catch (const std::exception& e)
    {
      failed++;
      std::cerr << "FAIL (should parse): " << input;
      if (!context.empty()) std::cerr << " [" << context << "]";
      std::cerr << " - " << e.what() << std::endl;
    }
  };

  // Helper lambda to test that a string parses successfully as ion (with charge)
  auto expectPassIon = [&](const std::string& input, const std::string& context = "")
  {
    try
    {
      ProForma::parseIon(input);
      passed++;
    }
    catch (const std::exception& e)
    {
      failed++;
      std::cerr << "FAIL (should parse): " << input;
      if (!context.empty()) std::cerr << " [" << context << "]";
      std::cerr << " - " << e.what() << std::endl;
    }
  };

  // Helper lambda to test that a string fails to parse as peptidoform
  auto expectFail = [&](const std::string& input, const std::string& context = "")
  {
    try
    {
      ProForma::parse(input);
      failed++;
      std::cerr << "FAIL (should reject): " << input;
      if (!context.empty()) std::cerr << " [" << context << "]";
      std::cerr << std::endl;
    }
    catch (const std::exception&)
    {
      passed++;
    }
  };

  // Helper lambda to test that a string fails to parse as ion
  auto expectFailIon = [&](const std::string& input, const std::string& context = "")
  {
    try
    {
      ProForma::parseIon(input);
      failed++;
      std::cerr << "FAIL (should reject): " << input;
      if (!context.empty()) std::cerr << " [" << context << "]";
      std::cerr << std::endl;
    }
    catch (const std::exception&)
    {
      passed++;
    }
  };

  // =========================================================
  // sequenceElement tests - single amino acids and modifications
  // =========================================================
  // Positive cases
  expectPass("A", "sequenceElement: single AA");
  expectPass("U", "sequenceElement: selenocysteine");
  expectPass("X", "sequenceElement: unknown AA");
  expectPass("A[#g1]", "sequenceElement: AA with label");

  // Note: "Asn" now parses as three lowercase amino acids (A-s-n) with lowercase support
  expectPass("Asn", "sequenceElement: lowercase letters now allowed");

  // Negative cases
  expectFail("A[#g1", "sequenceElement: unclosed bracket");

  // =========================================================
  // INT tests - integer parsing in charge context (use parseIon)
  // =========================================================
  // Positive cases
  expectPassIon("A/0", "INT: zero charge");
  expectPassIon("A/1", "INT: single digit");
  expectPassIon("A/2", "INT: double charge");

  // Negative cases
  expectFailIon("A/+1i", "INT: letter in number");

  // =========================================================
  // peptidoformCharge tests - charge state specifications
  // =========================================================
  // Positive cases (use parseIon for charge)
  expectPassIon("A/1", "charge: simple positive");
  expectPassIon("A/+1", "charge: explicit positive");
  expectPassIon("A/-1", "charge: negative");
  expectPassIon("A/2", "charge: double");
  expectPassIon("A/+2", "charge: explicit double positive");
  expectPassIon("A/-10101", "charge: large negative");
  expectPassIon("A/[Na:z+1]", "charge: single adduct");
  expectPassIon("A/[Na:z+1^2]", "charge: adduct with count");
  expectPassIon("A/[Na:z+1,Zn:z+2,H:z+1]", "charge: multiple adducts");

  // Negative cases
  expectFailIon("A/1/1", "charge: double charge spec");

  // =========================================================
  // adductIon tests - adduct ion specifications
  // =========================================================
  // Positive cases
  expectPassIon("A/[Na:z+1]", "adduct: basic");
  expectPassIon("A/[Na:z+1^2]", "adduct: with occurrence");

  // Negative cases
  expectFailIon("A/[Na^1]", "adduct: missing z");

  // =========================================================
  // formula tests - chemical formula parsing
  // =========================================================
  // Positive cases
  expectPass("A[Formula:N1H3]", "formula: explicit counts");
  expectPass("A[Formula:NH3]", "formula: implicit counts");
  expectPass("A[Formula:C-1O-1]", "formula: negative counts");
  expectPass("A[Formula:H1H1H1N1]", "formula: repeated elements");
  expectPass("A[Formula:[15N]H3]", "formula: isotope bracket");
  expectPass("A[Formula:[15N1]H3]", "formula: isotope with count");
  expectPass("A[Formula:C12H20O2]", "formula: standard");
  expectPass("A[Formula:HN-1O2]", "formula: negative in middle");
  expectPass("A[Formula:[13C2][12C-2]H2N]", "formula: multiple isotopes");
  expectPass("A[Formula:[13C2]C-2H2N]", "formula: isotope and negative");
  expectPass("A[Formula:Zn1]", "formula: zinc");

  // Negative cases
  expectFail("A[Formula:[15NH3]", "formula: unclosed isotope bracket");

  // =========================================================
  // modFormula tests - formula modifications with charge
  // =========================================================
  expectPass("A[Formula:Zn:z+2]", "modFormula: with charge");
  expectPass("A[Formula:Zn1:z+2]", "modFormula: explicit count with charge");
  expectPass("A[Formula:Na:z+1]", "modFormula: sodium with charge");

  // =========================================================
  // modMass tests - mass delta modifications
  // =========================================================
  // Positive cases
  expectPass("A[+23]", "modMass: positive integer");
  expectPass("A[-23.0]", "modMass: negative decimal");
  expectPass("A[+23.0]", "modMass: positive decimal");
  expectPass("A[U:+23]", "modMass: UNIMOD prefix");
  expectPass("A[Obs:+23]", "modMass: observed prefix");

  // =========================================================
  // mod tests - general modifications
  // =========================================================
  expectPass("A[+1]", "mod: simple mass");
  expectPass("A[U:+1]", "mod: UNIMOD mass");
  expectPass("A[+23.092]", "mod: decimal mass");
  expectPass("A[-23.092]", "mod: negative decimal");
  expectPass("A[14|Obs:+14|UNIMOD:0034|U:Methyl]", "mod: alternatives");
  expectPass("A[Formula:C-1O-1]", "mod: formula");
  expectPass("A[Glycan:Hex1HexNAc2]", "mod: glycan");
  expectPass("A[Formula:AlH-3:z+1]", "mod: formula with charge");
  expectPass("A[Formula:H-1:z-1]", "mod: negative charge");

  // =========================================================
  // modGlobal tests - global modification prefixes
  // =========================================================
  expectPass("<[TMT6plex]@K,N-term>A", "modGlobal: K and N-term");
  expectPass("<[TMT6plex]@K,N-terM>A", "modGlobal: case variation 1");
  expectPass("<[TMT6plex]@K,n-terM>A", "modGlobal: case variation 2");
  expectPass("<[TMT6plex]@K,n-tErM>A", "modGlobal: case variation 3");
  expectPass("<[TMT6plex]@K,N-term:A>AA", "modGlobal: with protein suffix");
  expectPass("<[TMT6plex]@K,N-term:A,N-term:B>AA", "modGlobal: multiple locations");

  // =========================================================
  // modGlycan tests - glycan modifications
  // =========================================================
  expectPass("A[Glycan:Hex2HexNAc]", "modGlycan: basic");
  expectPass("A[Glycan:Hex2HexNAc1]", "modGlycan: explicit count");

  // =========================================================
  // NAMETEXT tests - various name/text patterns
  // =========================================================
  // Positive cases (sequences that should parse)
  expectPass("III", "nametext: multiple I");

  // Negative cases
  expectFail("(KKK", "nametext: unclosed paren");

  // =========================================================
  // Summary
  // =========================================================
  TEST_EQUAL(failed, 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Mass Calculation tests
/////////////////////////////////////////////////////////////

START_SECTION(ProForma::canCalculateMass - simple sequence)
{
  Peptidoform pf = ProForma::parse("PEPTIDE");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)
}
END_SECTION

START_SECTION(ProForma::canCalculateMass - with UNIMOD modification)
{
  // Oxidation should resolve and allow mass calculation
  Peptidoform pf = ProForma::parse("PEM[UNIMOD:35]TIDE");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)
}
END_SECTION

START_SECTION(ProForma::canCalculateMass - with mass delta)
{
  // Mass delta should always allow mass calculation
  Peptidoform pf = ProForma::parse("PEM[+15.9949]TIDE");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)
}
END_SECTION

START_SECTION(ProForma::canCalculateMass - with formula)
{
  // Formula should allow mass calculation
  Peptidoform pf = ProForma::parse("PEM[Formula:O]TIDE");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)
}
END_SECTION

START_SECTION(ProForma::getMassCalculationIssues - unresolved modification)
{
  // Unknown modification name that won't resolve
  Peptidoform pf = ProForma::parse("PEM[UnknownMod12345]TIDE");
  auto issues = ProForma::getMassCalculationIssues(pf);
  TEST_EQUAL(issues.empty(), false)
  TEST_EQUAL(issues[0].type, ConversionIssueType::UNRESOLVED_MOD)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - simple sequence)
{
  Peptidoform pf = ProForma::parse("PEPTIDE");
  double mass = ProForma::getMonoWeight(pf);

  // Compare to AASequence calculation
  AASequence aas = AASequence::fromString("PEPTIDE");
  double expected = aas.getMonoWeight();

  TEST_REAL_SIMILAR(mass, expected)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - with UNIMOD modification)
{
  Peptidoform pf = ProForma::parse("PEM[UNIMOD:35]TIDE");
  double mass = ProForma::getMonoWeight(pf);

  // Compare to AASequence calculation
  AASequence aas = AASequence::fromString("PEM(Oxidation)TIDE");
  double expected = aas.getMonoWeight();

  TEST_REAL_SIMILAR(mass, expected)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - with mass delta)
{
  Peptidoform pf = ProForma::parse("PEM[+15.9949]TIDE");
  double mass = ProForma::getMonoWeight(pf);

  // Compare to unmodified + mass delta
  AASequence aas = AASequence::fromString("PEMTIDE");
  double expected = aas.getMonoWeight() + 15.9949;

  TEST_REAL_SIMILAR(mass, expected)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - with terminal modifications)
{
  Peptidoform pf = ProForma::parse("[UNIMOD:1]-PEPTIDE-[UNIMOD:2]");
  double mass = ProForma::getMonoWeight(pf);

  // Compare to AASequence with Acetyl (N-term) and Amidated (C-term)
  AASequence aas = AASequence::fromString("(Acetyl)PEPTIDE(Amidated)");
  double expected = aas.getMonoWeight();

  TEST_REAL_SIMILAR(mass, expected)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - cross-linked single chain)
{
  // Cross-link within single chain using mass delta (disulfide bond: -2.0156 Da for loss of 2H)
  // The cross-link mass should only be counted once, not at both sites
  Peptidoform pf = ProForma::parse("EVTSEKC[-2.0156#XL1]LEMSC[#XL1]EFD");

  // Should be able to calculate mass
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)

  double mass = ProForma::getMonoWeight(pf);

  // Calculate expected mass: sequence + one disulfide delta (not two)
  // Sequence is EVTSEKCLEMSCEFD (15 amino acids)
  AASequence unmod_seq = AASequence::fromString("EVTSEKCLEMSCEFD");
  double expected = unmod_seq.getMonoWeight() - 2.0156;  // One cross-link, not two

  TEST_REAL_SIMILAR(mass, expected)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - multi-chain cross-link)
{
  // Cross-linked peptides using mass delta for DSS/BS3 cross-linker (+138.06807961 Da)
  // This is the exact mass from CrossLinksDB_test.cpp
  PeptidoformIon ion = ProForma::parseIon("PEPTIDEK[+138.06807961#XL1]//SEQUENCEK[#XL1]");

  TEST_EQUAL(ProForma::canCalculateMass(ion), true)

  double mass = ProForma::getMonoWeight(ion);

  // Expected: sum of both peptide masses + one cross-linker mass
  AASequence seq1 = AASequence::fromString("PEPTIDEK");
  AASequence seq2 = AASequence::fromString("SEQUENCEK");
  double expected = seq1.getMonoWeight() + seq2.getMonoWeight() + 138.06807961;

  TEST_REAL_SIMILAR(mass, expected)
}
END_SECTION

START_SECTION(ProForma::getMZ - with charge state)
{
  PeptidoformIon ion = ProForma::parseIon("PEPTIDE/2");
  double mz = ProForma::getMZ(ion);

  // Compare to AASequence calculation
  AASequence aas = AASequence::fromString("PEPTIDE");
  double expected = aas.getMZ(2);

  TEST_REAL_SIMILAR(mz, expected)
}
END_SECTION

START_SECTION(ProForma::getMZ - Peptidoform with charge)
{
  Peptidoform pf = ProForma::parse("PEPTIDE");
  double mz = ProForma::getMZ(pf, 2);

  // Compare to AASequence calculation
  AASequence aas = AASequence::fromString("PEPTIDE");
  double expected = aas.getMZ(2);

  TEST_REAL_SIMILAR(mz, expected)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - throws on unresolved)
{
  Peptidoform pf = ProForma::parse("PEM[UnknownModification99999]TIDE");
  TEST_EXCEPTION(Exception::InvalidValue, ProForma::getMonoWeight(pf))
}
END_SECTION

START_SECTION(ProForma::getMZ - throws on zero charge)
{
  Peptidoform pf = ProForma::parse("PEPTIDE");
  TEST_EXCEPTION(Exception::InvalidValue, ProForma::getMZ(pf, 0))
}
END_SECTION

START_SECTION(ProForma::getMZ - throws on missing charge)
{
  PeptidoformIon ion = ProForma::parseIon("PEPTIDE");
  TEST_EXCEPTION(Exception::InvalidValue, ProForma::getMZ(ion))
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - reference peptide DFPIANGER)
{
  // DFPIANGER is a well-tested reference peptide in OpenMS (AASequence_test.cpp line 384)
  Peptidoform pf = ProForma::parse("DFPIANGER");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)

  double mass = ProForma::getMonoWeight(pf);
  AASequence aas = AASequence::fromString("DFPIANGER");
  TEST_REAL_SIMILAR(mass, aas.getMonoWeight())
  TEST_REAL_SIMILAR(mass, 1017.48796)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - phosphorylation)
{
  // Phosphorylation: +79.966331 Da (UNIMOD:21)
  Peptidoform pf = ProForma::parse("PEPTS[UNIMOD:21]IDE");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)

  double mass = ProForma::getMonoWeight(pf);
  AASequence aas = AASequence::fromString("PEPTS(Phospho)IDE");
  TEST_REAL_SIMILAR(mass, aas.getMonoWeight())

  // Also test with mass delta notation
  Peptidoform pf2 = ProForma::parse("PEPTS[+79.966331]IDE");
  TEST_REAL_SIMILAR(ProForma::getMonoWeight(pf2), mass)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - carbamidomethyl)
{
  // Carbamidomethyl on Cysteine: +57.0215 Da (UNIMOD:4)
  Peptidoform pf = ProForma::parse("PEPTC[UNIMOD:4]IDE");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)

  double mass = ProForma::getMonoWeight(pf);
  AASequence aas = AASequence::fromString("PEPTC(Carbamidomethyl)IDE");
  TEST_REAL_SIMILAR(mass, aas.getMonoWeight())
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - combined terminal mods)
{
  // N-term Acetyl (+42.010565) and C-term Amidated (-0.984016)
  Peptidoform pf = ProForma::parse("[UNIMOD:1]-PEPTIDE-[UNIMOD:2]");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)

  double mass = ProForma::getMonoWeight(pf);
  AASequence aas = AASequence::fromString("(Acetyl)PEPTIDE(Amidated)");
  TEST_REAL_SIMILAR(mass, aas.getMonoWeight())

  // Verify the mass difference from unmodified
  Peptidoform pf_unmod = ProForma::parse("PEPTIDE");
  double unmod_mass = ProForma::getMonoWeight(pf_unmod);
  // Acetyl adds ~42.01, Amidated subtracts ~0.98, net ~41.03
  TEST_REAL_SIMILAR(mass - unmod_mass, 42.010565 - 0.984016)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - multiple modifications)
{
  // Multiple modifications on same peptide
  Peptidoform pf = ProForma::parse("[UNIMOD:1]-PEM[UNIMOD:35]PTIDES[UNIMOD:21]K");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)

  double mass = ProForma::getMonoWeight(pf);
  AASequence aas = AASequence::fromString("(Acetyl)PEM(Oxidation)PTIDES(Phospho)K");
  TEST_REAL_SIMILAR(mass, aas.getMonoWeight())
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - formula tag)
{
  // Formula:C2H2O adds acetyl-equivalent mass
  Peptidoform pf = ProForma::parse("[Formula:C2H2O]-PEPTIDE");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)

  double mass = ProForma::getMonoWeight(pf);

  // C2H2O has same composition as Acetyl
  Peptidoform pf_acetyl = ProForma::parse("[UNIMOD:1]-PEPTIDE");
  TEST_REAL_SIMILAR(mass, ProForma::getMonoWeight(pf_acetyl))
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - unlocalised modification)
{
  // Unlocalised phosphorylation: [Phospho]?PEPTIDE
  Peptidoform pf = ProForma::parse("[+79.966331]?PEPTIDE");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)

  double mass = ProForma::getMonoWeight(pf);
  Peptidoform pf_base = ProForma::parse("PEPTIDE");
  double base_mass = ProForma::getMonoWeight(pf_base);

  // Mass should include the unlocalised phospho
  TEST_REAL_SIMILAR(mass, base_mass + 79.966331)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - labile modification)
{
  // Labile modification: {Glycan:Hex}PEPTIDE
  // Using mass delta since Glycan may not resolve
  Peptidoform pf = ProForma::parse("{+162.0528}PEPTIDE");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)

  double mass = ProForma::getMonoWeight(pf);
  Peptidoform pf_base = ProForma::parse("PEPTIDE");
  double base_mass = ProForma::getMonoWeight(pf_base);

  // Labile mods are included in mass calculation
  TEST_REAL_SIMILAR(mass, base_mass + 162.0528)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - global modification)
{
  // Global Carbamidomethyl on all C residues: <[+57.0215]@C>PEPTCIDECK
  Peptidoform pf = ProForma::parse("<[+57.0215]@C>PEPTCIDECK");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)

  double mass = ProForma::getMonoWeight(pf);

  // Should be equivalent to explicit modifications on both C residues
  Peptidoform pf_explicit = ProForma::parse("PEPTC[+57.0215]IDEC[+57.0215]K");
  TEST_REAL_SIMILAR(mass, ProForma::getMonoWeight(pf_explicit))

  // Verify: base + 2 * carbamidomethyl
  Peptidoform pf_base = ProForma::parse("PEPTCIDECK");
  double base_mass = ProForma::getMonoWeight(pf_base);
  TEST_REAL_SIMILAR(mass, base_mass + 2 * 57.0215)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - disulfide cross-link)
{
  // Disulfide bond: -2.01565 Da (loss of 2H)
  // Single chain with intramolecular disulfide
  Peptidoform pf = ProForma::parse("PEPTC[-2.01565#XL1]IDEC[#XL1]K");
  TEST_EQUAL(ProForma::canCalculateMass(pf), true)

  double mass = ProForma::getMonoWeight(pf);
  Peptidoform pf_base = ProForma::parse("PEPTCIDECK");
  double base_mass = ProForma::getMonoWeight(pf_base);

  // Disulfide removes 2H, counted once
  TEST_REAL_SIMILAR(mass, base_mass - 2.01565)
}
END_SECTION

START_SECTION(ProForma::getMonoWeight - throws on chimeric)
{
  // Chimeric spectra have ambiguous mass - must calculate per chain
  PeptidoformIon ion = ProForma::parseIon("PEPTIDE/2+EDITPEP/3");
  TEST_EQUAL(ion.is_chimeric, true)
  TEST_EXCEPTION(Exception::InvalidValue, ProForma::getMonoWeight(ion))

  // Individual chains should work fine
  TEST_EQUAL(ProForma::canCalculateMass(ion.chains[0]), true)
  TEST_EQUAL(ProForma::canCalculateMass(ion.chains[1]), true)
  double mass1 = ProForma::getMonoWeight(ion.chains[0]);
  double mass2 = ProForma::getMonoWeight(ion.chains[1]);
  TEST_REAL_SIMILAR(mass1, 799.3599) // PEPTIDE
  TEST_REAL_SIMILAR(mass2, 799.3599) // EDITPEP (same residues)
}
END_SECTION

START_SECTION(ProForma::tryGetMonoWeight - basic usage)
{
  // Simple case - should succeed
  Peptidoform pf = ProForma::parse("PEPTIDE");
  auto mass = ProForma::tryGetMonoWeight(pf);
  TEST_EQUAL(mass.has_value(), true)
  TEST_REAL_SIMILAR(*mass, 799.3599)

  // Compare with throwing version
  TEST_REAL_SIMILAR(*mass, ProForma::getMonoWeight(pf))
}
END_SECTION

START_SECTION(ProForma::tryGetMonoWeight - with issues)
{
  // Valid peptide - no issues
  Peptidoform pf = ProForma::parse("PEM[UNIMOD:35]TIDE");
  std::vector<ConversionIssue> issues;
  auto mass = ProForma::tryGetMonoWeight(pf, issues);

  TEST_EQUAL(mass.has_value(), true)
  TEST_EQUAL(issues.empty(), true)

  // Invalid modification - should return nullopt and populate issues
  Peptidoform pf_bad = ProForma::parse("PEM[UnknownMod999]TIDE");
  issues.clear();
  auto mass_bad = ProForma::tryGetMonoWeight(pf_bad, issues);

  TEST_EQUAL(mass_bad.has_value(), false)
  TEST_EQUAL(issues.empty(), false)
  // Should have at least one issue about unresolved modification
  TEST_EQUAL(issues[0].type, ConversionIssueType::UNRESOLVED_MOD)
}
END_SECTION

START_SECTION(ProForma::tryGetMonoWeight - PeptidoformIon)
{
  // Cross-linked peptides - should work
  PeptidoformIon ion = ProForma::parseIon("PEPTIDEK[+138.068#XL1]//SEQUENCEK[#XL1]");
  auto mass = ProForma::tryGetMonoWeight(ion);
  TEST_EQUAL(mass.has_value(), true)
  TEST_REAL_SIMILAR(*mass, ProForma::getMonoWeight(ion))

  // Chimeric - should return nullopt
  PeptidoformIon chimeric = ProForma::parseIon("PEPTIDE/2+EDITPEP/3");
  std::vector<ConversionIssue> issues;
  auto mass_chimeric = ProForma::tryGetMonoWeight(chimeric, issues);
  TEST_EQUAL(mass_chimeric.has_value(), false)
  TEST_EQUAL(issues.empty(), false)
}
END_SECTION

START_SECTION(ProForma::tryGetMZ - basic usage)
{
  // Peptidoform with explicit charge
  Peptidoform pf = ProForma::parse("PEPTIDE");
  auto mz = ProForma::tryGetMZ(pf, 2);
  TEST_EQUAL(mz.has_value(), true)
  TEST_REAL_SIMILAR(*mz, ProForma::getMZ(pf, 2))

  // Zero charge - should return nullopt
  auto mz_zero = ProForma::tryGetMZ(pf, 0);
  TEST_EQUAL(mz_zero.has_value(), false)

  // With issues output
  std::vector<ConversionIssue> issues;
  auto mz_issues = ProForma::tryGetMZ(pf, 0, issues);
  TEST_EQUAL(mz_issues.has_value(), false)
  TEST_EQUAL(issues.empty(), false)
}
END_SECTION

START_SECTION(ProForma::tryGetMZ - PeptidoformIon)
{
  // PeptidoformIon with charge
  PeptidoformIon ion = ProForma::parseIon("PEPTIDE/2");
  auto mz = ProForma::tryGetMZ(ion);
  TEST_EQUAL(mz.has_value(), true)
  TEST_REAL_SIMILAR(*mz, ProForma::getMZ(ion))

  // Without charge - should return nullopt
  PeptidoformIon ion_no_charge = ProForma::parseIon("PEPTIDE");
  std::vector<ConversionIssue> issues;
  auto mz_no_charge = ProForma::tryGetMZ(ion_no_charge, issues);
  TEST_EQUAL(mz_no_charge.has_value(), false)
  TEST_EQUAL(issues.empty(), false)
}
END_SECTION

START_SECTION(ProForma::tryGetMonoWeight - efficiency test)
{
  // Verify that tryGetMonoWeight does single-pass (doesn't throw, returns same result)
  // This is more of a usage pattern demonstration
  Peptidoform pf = ProForma::parse("PEM[UNIMOD:35]PTIDES[UNIMOD:21]K");

  // Efficient pattern: single call gets result or issues
  std::vector<ConversionIssue> issues;
  if (auto mass = ProForma::tryGetMonoWeight(pf, issues))
  {
    // Use the mass
    TEST_REAL_SIMILAR(*mass, ProForma::getMonoWeight(pf))
  }
  else
  {
    // Would examine issues here
    TEST_EQUAL(true, false) // Should not reach here for valid peptide
  }
}
END_SECTION

START_SECTION(ProForma::canGenerateSpectrum - Peptidoform)
{
  // Test canGenerateSpectrum for resolved peptide
  Peptidoform pf = ProForma::parse("PEPTIDE");
  ProForma::resolveModifications(pf);
  TEST_EQUAL(ProForma::canGenerateSpectrum(pf), true)

  // Unresolved modification should fail
  Peptidoform pf_unresolved = ProForma::parse("PEM[UnknownMod999]PTIDE");
  TEST_EQUAL(ProForma::canGenerateSpectrum(pf_unresolved), false)
}
END_SECTION

START_SECTION(ProForma::canGenerateSpectrum - PeptidoformIon)
{
  // Single chain should work
  PeptidoformIon pfi = ProForma::parseIon("PEPTIDE/2");
  ProForma::resolveModifications(pfi.chains[0]);
  TEST_EQUAL(ProForma::canGenerateSpectrum(pfi), true)

  // Cross-linked should work
  PeptidoformIon pfi_xl = ProForma::parseIon("PEPK[+138.068#XL1]IDE//ANOK[#XL1]THER");
  ProForma::resolveModifications(pfi_xl.chains[0]);
  ProForma::resolveModifications(pfi_xl.chains[1]);
  TEST_EQUAL(ProForma::canGenerateSpectrum(pfi_xl), true)

  // Chimeric (no cross-link) should fail
  PeptidoformIon pfi_chimeric = ProForma::parseIon("PEPTIDE//ANOTHER");
  TEST_EQUAL(ProForma::canGenerateSpectrum(pfi_chimeric), false)
}
END_SECTION

START_SECTION(ProForma::getSpectrumGenerationIssues - Peptidoform)
{
  // Resolved peptide should have no issues
  Peptidoform pf = ProForma::parse("PEPTIDE");
  ProForma::resolveModifications(pf);
  auto issues = ProForma::getSpectrumGenerationIssues(pf);
  TEST_EQUAL(issues.empty(), true)

  // Unresolved modification should have issues
  Peptidoform pf_unresolved = ProForma::parse("PEM[UnknownMod999]PTIDE");
  auto issues_unresolved = ProForma::getSpectrumGenerationIssues(pf_unresolved);
  TEST_EQUAL(issues_unresolved.empty(), false)
}
END_SECTION

START_SECTION(ProForma::getSpectrumGenerationIssues - PeptidoformIon)
{
  // Chimeric should have issues
  PeptidoformIon pfi_chimeric = ProForma::parseIon("PEPTIDE//ANOTHER");
  auto issues = ProForma::getSpectrumGenerationIssues(pfi_chimeric);
  TEST_EQUAL(issues.empty(), false)
}
END_SECTION

START_SECTION(ProForma::generateSpectrum - basic peptide)
{
  // Test spectrum generation for a simple peptide
  Peptidoform pf = ProForma::parse("PEPTIDE");
  ProForma::resolveModifications(pf);

  MSSpectrum spec = ProForma::generateSpectrum(pf, 1, 1, "by", false, true);
  // Should have fragments
  TEST_EQUAL(spec.size() > 0, true)
  // Should have b and y ions at minimum
  TEST_EQUAL(spec.size() >= 10, true)
}
END_SECTION

START_SECTION(ProForma::generateSpectrum - with modifications)
{
  // Test spectrum generation for modified peptide
  Peptidoform pf = ProForma::parse("PEM[UNIMOD:35]PTIDE");
  ProForma::resolveModifications(pf);

  MSSpectrum spec = ProForma::generateSpectrum(pf, 1, 2, "by", true, true);
  TEST_EQUAL(spec.size() > 0, true)
}
END_SECTION

START_SECTION(ProForma::generateSpectrum - ion types)
{
  // Test different ion types
  Peptidoform pf = ProForma::parse("PEPTIDE");
  ProForma::resolveModifications(pf);

  // Default by ions
  MSSpectrum spec_by = ProForma::generateSpectrum(pf, 1, 1, "by", false, true);
  // With precursor
  MSSpectrum spec_byM = ProForma::generateSpectrum(pf, 1, 1, "byM", false, true);
  // Precursor adds peaks
  TEST_EQUAL(spec_byM.size() >= spec_by.size(), true)
}
END_SECTION

START_SECTION(ProForma::generateSpectrum - PeptidoformIon single chain)
{
  // Test spectrum generation for PeptidoformIon with single chain
  PeptidoformIon pfi = ProForma::parseIon("PEPTIDE/2");
  ProForma::resolveModifications(pfi.chains[0]);

  MSSpectrum spec = ProForma::generateSpectrum(pfi, 1, 2, "by", false, true);
  TEST_EQUAL(spec.size() > 0, true)
}
END_SECTION

START_SECTION(ProForma::generateSpectrum - cross-linked peptides)
{
  // Test spectrum generation for cross-linked peptides
  // DSS cross-linker at +138.068 Da
  PeptidoformIon pfi = ProForma::parseIon("PEPK[+138.068#XL1]IDE//ANOK[#XL1]THER");
  ProForma::resolveModifications(pfi.chains[0]);
  ProForma::resolveModifications(pfi.chains[1]);

  MSSpectrum spec = ProForma::generateSpectrum(pfi, 1, 2, "aby", false, true);
  // Should succeed for properly formed cross-link
  TEST_EQUAL(spec.size() > 0, true)
}
END_SECTION

START_SECTION(ProForma::generateSpectrum - fails for unsupported cases)
{
  // Test that unresolved modifications cannot generate spectrum
  // (We use canGenerateSpectrum check rather than calling generateSpectrum
  // to avoid ModificationsDB error messages that fail the test framework)
  Peptidoform pf = ProForma::parse("PEM[UnknownMod999]PTIDE");
  TEST_EQUAL(ProForma::canGenerateSpectrum(pf), false)
  auto issues = ProForma::getSpectrumGenerationIssues(pf);
  TEST_EQUAL(issues.empty(), false)

  // Chimeric (multiple chains without cross-link) cannot generate spectrum
  PeptidoformIon pfi = ProForma::parseIon("PEPTIDE//ANOTHER");
  TEST_EQUAL(ProForma::canGenerateSpectrum(pfi), false)
  auto issues2 = ProForma::getSpectrumGenerationIssues(pfi);
  TEST_EQUAL(issues2.empty(), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
