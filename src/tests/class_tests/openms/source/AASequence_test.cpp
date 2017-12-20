// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <iostream>
#include <OpenMS/SYSTEM/StopWatch.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(AASequence, "$Id$")

/////////////////////////////////////////////////////////////

// AASequence seq = AASequence::fromString("PEPTM(Met->Hse)IDE");
// cout << seq.size() << endl;
// cout << ModificationsDB::getInstance()->getModification("Acetyl", "", ResidueModification::N_TERM).getOrigin() << endl;

AASequence* ptr = nullptr;
AASequence* nullPointer = nullptr;
START_SECTION(AASequence())
  ptr = new AASequence();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~AASequence())
  delete ptr;
END_SECTION

START_SECTION(AASequence(const AASequence& rhs))
  AASequence seq;
  seq = AASequence::fromString("AAA");
  AASequence seq2(seq);
  TEST_EQUAL(seq, seq2)
END_SECTION

START_SECTION(AASequence fromString(const String& s, bool permissive = true))
{
  AASequence seq = AASequence::fromString("CNARCKNCNCNARCDRE");
  TEST_EQUAL(seq.isModified(), false)
  TEST_EQUAL(seq.hasNTerminalModification(), false);
  TEST_EQUAL(seq.hasCTerminalModification(), false);
  TEST_EQUAL(seq.getResidue(4).getModification(), 0);

  AASequence seq2;
  seq2 = AASequence::fromString("CNARCKNCNCNARCDRE");
  TEST_EQUAL(seq, seq2);

  // test complex term-mods
  AASequence seq3 = AASequence::fromString("VPQVSTPTLVEVSRSLGK(Label:18O(2))");
  TEST_EQUAL(seq3.isModified(), true)
  TEST_EQUAL(seq3.hasNTerminalModification(), false);
  TEST_EQUAL(seq3.hasCTerminalModification(), true);
  TEST_EQUAL(seq3.getResidue(4).getModification(), 0);
  TEST_EQUAL(seq3.getCTerminalModificationName(), "Label:18O(2)");
  AASequence seq4 = AASequence::fromString("VPQVSTPTLVEVSRSLGK(Label:18O(2))");
  TEST_EQUAL(seq3, seq4);

  AASequence seq5 = AASequence::fromString("(ICPL:2H(4))CNARCNCNCN");
  TEST_EQUAL(seq5.hasNTerminalModification(), true);
  TEST_EQUAL(seq5.isModified(), true);
  TEST_EQUAL(seq5.getNTerminalModificationName(), "ICPL:2H(4)");

  AASequence seq6 = AASequence::fromString("CNARCK(Label:13C(6)15N(2))NCNCN");
  TEST_EQUAL(seq6.hasNTerminalModification(), false);
  TEST_EQUAL(seq6.hasCTerminalModification(), false);
  TEST_EQUAL(seq6.isModified(), true);
  TEST_EQUAL(seq6.getResidue(5).getModificationName(), "Label:13C(6)15N(2)");
  TEST_EQUAL(seq6.getResidue(4).getModificationName(), "");

  AASequence seq7 = AASequence::fromString("CNARCKNCNCNARCDRE(Amidated)");
  TEST_EQUAL(seq7.hasNTerminalModification(), false);
  TEST_EQUAL(seq7.hasCTerminalModification(), true);
  TEST_EQUAL(seq7.isModified(), true);
  TEST_EQUAL(seq7.getCTerminalModificationName(), "Amidated");

  // test square bracket modifications
  AASequence seq8 = AASequence::fromString("PEPTIDEK[136]");
  TEST_EQUAL(seq8.hasNTerminalModification(), false);
  TEST_EQUAL(seq8.hasCTerminalModification(), false);
  TEST_EQUAL(seq8.isModified(), true);
  TEST_STRING_EQUAL(seq8[7].getModificationName(), "Label:13C(6)15N(2)");

  AASequence seq9 = AASequence::fromString("PEPS[167]TIDEK");
  TEST_EQUAL(seq9.isModified(), true);
  TEST_STRING_EQUAL(seq9[3].getModificationName(), "Phospho");

  AASequence seq10 = AASequence::fromString("PEPC[160]TIDEK");
  TEST_EQUAL(seq10.isModified(), true);
  TEST_STRING_EQUAL(seq10[3].getModificationName(), "Carbamidomethyl");

  AASequence seq11 = AASequence::fromString("PEPM[147]TIDEK");
  TEST_EQUAL(seq11.isModified(), true);
  TEST_STRING_EQUAL(seq11[3].getModificationName(), "Oxidation");

  AASequence seq12 = AASequence::fromString("PEPT[181]TIDEK");
  TEST_EQUAL(seq12.isModified(), true);
  TEST_STRING_EQUAL(seq12[3].getModificationName(), "Phospho");

  AASequence seq13 = AASequence::fromString("PEPY[243]TIDEK");
  TEST_EQUAL(seq13.isModified(), true);
  TEST_STRING_EQUAL(seq13[3].getModificationName(), "Phospho");

  AASequence seq14 = AASequence::fromString("PEPR[166]TIDEK");
  TEST_EQUAL(seq14.isModified(), true);
  TEST_STRING_EQUAL(seq14[3].getModificationName(), "Label:13C(6)15N(4)");

  // invalid test case: "Pyro-carbamidomethyl" is only defined as N-terminal
  // AASequence seq15 = AASequence::fromString("PEPC[143]TIDEK");
  // TEST_EQUAL(seq15.isModified(), true);
  // TEST_STRING_EQUAL(seq15[3].getModificationName(), "Pyro-carbamidomethyl");

  // invalid test case: "Gln->pyro-Glu" is only defined as N-terminal
  // AASequence seq16 = AASequence::fromString("PEPQ[111]TIDEK");
  // TEST_EQUAL(seq16.isModified(), true);
  // TEST_STRING_EQUAL(seq16[3].getModificationName(), "Gln->pyro-Glu");

  // invalid test case: "Glu->pyro-Glu" is only defined as N-terminal
  // AASequence seq17 = AASequence::fromString("PEPE[111]TIDEK");
  // TEST_EQUAL(seq17.isModified(), true);
  // TEST_STRING_EQUAL(seq17[3].getModificationName(), "Glu->pyro-Glu");

  TEST_EXCEPTION(Exception::ParseError, AASequence::fromString("blDABCDEF"));
  TEST_EXCEPTION(Exception::ParseError, AASequence::fromString("a"));

  // test "permissive" option:
  AASequence seq18 = AASequence::fromString("PEP T*I#D+E", true);
  TEST_EQUAL(seq18.size(), 10);
  TEST_EQUAL(seq18.toString(), "PEPTXIXDXE");

  TEST_EXCEPTION(Exception::ParseError,
                 AASequence::fromString("PEP T*I#D+E", false));

  // prefer residue mod. over C-term mod.:
  AASequence seq19 = AASequence::fromString("PEPM(Oxidation)");
  TEST_EQUAL(seq19.hasCTerminalModification(), false);
  TEST_EQUAL(seq19[3].isModified(), true);
  TEST_STRING_EQUAL(seq19[3].getModificationName(), "Oxidation");
}
END_SECTION

START_SECTION(AASequence& operator=(const AASequence& rhs))
  AASequence seq = AASequence::fromString("AAA");
  AASequence seq2 = AASequence::fromString("AAA");
  TEST_EQUAL(seq, seq2)
END_SECTION

START_SECTION(([EXTRA]Test modifications with brackets))
  AASequence seq1 = AASequence::fromString("ANLVFK(Label:13C(6)15N(2))EIEK(Label:2H(4))");
  TEST_EQUAL(seq1.hasNTerminalModification(), false)
  TEST_EQUAL(seq1.hasCTerminalModification(), false)
  TEST_EQUAL(seq1.isModified(), true)
  AASequence seq2 = AASequence::fromString("ANLVFK(Label:13C(6)15N(2))EIEK(Label:2H(4))(Amidated)");
  TEST_EQUAL(seq2.hasNTerminalModification(), false)
  TEST_EQUAL(seq2.hasCTerminalModification(), true)
  TEST_EQUAL(seq2.isModified(), true)
END_SECTION

START_SECTION(bool operator==(const AASequence& rhs) const)
  AASequence seq1 = AASequence::fromString("(Acetyl)DFPIANGER");
  AASequence seq2 = AASequence::fromString("DFPIANGER");
  TEST_EQUAL(seq2 == AASequence::fromString("DFPIANGER"), true)
  TEST_EQUAL(seq1 == AASequence::fromString("(Acetyl)DFPIANGER"), true)

  AASequence seq3 = AASequence::fromString("DFPIANGER(ADP-Ribosyl)");
  AASequence seq4 = AASequence::fromString("DFPIANGER(Amidated)");
  TEST_EQUAL(seq3 == AASequence::fromString("DFPIANGER"), false)
  TEST_EQUAL(seq3 == AASequence::fromString("DFPIANGER(ADP-Ribosyl)"), true)
  TEST_EQUAL(seq4 == AASequence::fromString("DFPIANGER(Amidated)"), true)
  TEST_EQUAL(seq4 == AASequence::fromString("DFPIANGER"), false)

  AASequence seq5 = AASequence::fromString("DFBIANGER");
  TEST_EQUAL(seq5 == AASequence::fromString("DFPIANGER"), false)
  TEST_EQUAL(seq5 == AASequence::fromString("DFBIANGER"), true)
END_SECTION

START_SECTION(const Residue& getResidue(Size index) const)
  AASequence seq = AASequence::fromString("ACDEF");
  Size unsignedint(2);
  TEST_EQUAL(seq.getResidue(unsignedint).getOneLetterCode(), "D")
  TEST_EXCEPTION(Exception::IndexOverflow, seq.getResidue((Size)1000))
END_SECTION

START_SECTION((EmpiricalFormula getFormula(Residue::ResidueType type = Residue::Full, Int charge=0) const))
  AASequence seq = AASequence::fromString("ACDEF");
  TEST_EQUAL(seq.getFormula(), EmpiricalFormula("O10SH33N5C24"))
  TEST_EQUAL(seq.getFormula(Residue::Full, 1), EmpiricalFormula("O10SH33N5C24+"))
  TEST_EQUAL(seq.getFormula(Residue::BIon, 0), EmpiricalFormula("O9SH31N5C24"))
END_SECTION

START_SECTION((double getAverageWeight(Residue::ResidueType type = Residue::Full, Int charge=0) const))
  AASequence seq = AASequence::fromString("DFPIANGER");
  TOLERANCE_ABSOLUTE(0.01)
  TEST_REAL_SIMILAR(seq.getAverageWeight(), double(1018.08088))
  TEST_REAL_SIMILAR(seq.getAverageWeight(Residue::YIon, 1), double(1019.09))
END_SECTION

START_SECTION((double getMonoWeight(Residue::ResidueType type = Residue::Full, Int charge=0) const))
  TOLERANCE_ABSOLUTE(1e-6)
  TOLERANCE_RELATIVE(1.0 + 1e-6)

  // test if fragments of charged single amino acid sequences match the charged residue weight of the fragment ions
  EmpiricalFormula ala_res = EmpiricalFormula("C3H5NO");
  TEST_REAL_SIMILAR(AASequence::fromString("A").getMonoWeight(Residue::Internal, 0), ala_res.getMonoWeight());

  EmpiricalFormula ala_full = ala_res + EmpiricalFormula("H2O");
  TEST_REAL_SIMILAR(AASequence::fromString("A").getMonoWeight(Residue::Full, 0), ala_full.getMonoWeight());

  EmpiricalFormula ala_a_neutral = EmpiricalFormula("H")+ala_res-EmpiricalFormula("CHO");
  TEST_REAL_SIMILAR(AASequence::fromString("A").getMonoWeight(Residue::AIon, 1), ala_a_neutral.getMonoWeight()+Constants::PROTON_MASS_U);
  //44.04947

  EmpiricalFormula ala_b_neutral = EmpiricalFormula("H")+ala_res-EmpiricalFormula("H");
  TEST_REAL_SIMILAR(AASequence::fromString("A").getMonoWeight(Residue::BIon, 1), ala_b_neutral.getMonoWeight()+Constants::PROTON_MASS_U);
  //72.04439

  EmpiricalFormula ala_y_neutral = EmpiricalFormula("OH")+ala_res+EmpiricalFormula("H");
  TEST_REAL_SIMILAR(AASequence::fromString("A").getMonoWeight(Residue::YIon, 1), ala_y_neutral.getMonoWeight()+Constants::PROTON_MASS_U);
  //90.05496

  EmpiricalFormula ala_z_neutral = EmpiricalFormula("OH")+ala_res-EmpiricalFormula("NH2");
  TEST_REAL_SIMILAR(AASequence::fromString("A").getMonoWeight(Residue::ZIon, 1), ala_z_neutral.getMonoWeight()+Constants::PROTON_MASS_U);
  //73.02900


  TEST_REAL_SIMILAR(AASequence::fromString("DFPIANGER").getMonoWeight(), double(1017.48796))

  // test if direct calculation and calculation via empirical formula yield the same result
  TEST_REAL_SIMILAR(AASequence::fromString("DFPIANGER").getMonoWeight(Residue::YIon, 1), AASequence::fromString("DFPIANGER").getFormula(Residue::YIon, 1).getMonoWeight())

  TEST_REAL_SIMILAR(AASequence::fromString("DFPIANGER").getMonoWeight(Residue::YIon, 1), double(1018.4952))

  // test N-term modification
  AASequence seq2 = AASequence::fromString("(NIC)DFPIANGER");
  TEST_REAL_SIMILAR(seq2.getMonoWeight(), double(1122.51));

  // test old OpenMS NIC definition
  AASequence seq2a = AASequence::fromString("(MOD:09998)DFPIANGER");
  TEST_EQUAL(seq2 == seq2a, true)

  // test heavy modification
  AASequence seq3 = AASequence::fromString("(dNIC)DFPIANGER");
  TEST_REAL_SIMILAR(seq3.getMonoWeight(), double(1017.48796) + double(109.048119));

  // test old OpenMS dNIC definition
  AASequence seq3a = AASequence::fromString("(MOD:09999)DFPIANGER");
  TEST_EQUAL(seq3 == seq3a, true)

  TEST_REAL_SIMILAR(AASequence::fromString("TYQYS(Phospho)").getFormula().getMonoWeight(), AASequence::fromString("TYQYS(Phospho)").getMonoWeight());

  TEST_REAL_SIMILAR(AASequence::fromString("TYQYS(Phospho)").getFormula().getMonoWeight(), AASequence::fromString("TYQYS(Phospho)").getMonoWeight());
END_SECTION

START_SECTION(const Residue& operator[](Size index) const)
  AASequence seq = AASequence::fromString("DFPIANGER");
  Size index = 0;
  TEST_EQUAL(seq[index].getOneLetterCode(), "D")
  index = 20;
  TEST_EXCEPTION(Exception::IndexOverflow, seq[index])
END_SECTION

START_SECTION(AASequence operator+(const AASequence& peptide) const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("DFP");
  AASequence seq3 = AASequence::fromString("IANGER");
  TEST_EQUAL(seq1, seq2 + seq3);
END_SECTION

START_SECTION(AASequence operator+(const Residue* residue) const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("DFPIANGE");
  TEST_EQUAL(seq1, seq2 + ResidueDB::getInstance()->getResidue("R"))
END_SECTION

START_SECTION(AASequence& operator+=(const AASequence&))
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("DFP");
  AASequence seq3 = AASequence::fromString("IANGER");
  seq2 += seq3;
  TEST_EQUAL(seq1, seq2)
END_SECTION

START_SECTION(AASequence& operator+=(const Residue* residue))
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("DFPIANGE");
  seq2 += ResidueDB::getInstance()->getResidue("R");
  TEST_EQUAL(seq1, seq2)
END_SECTION

START_SECTION(Size size() const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  TEST_EQUAL(seq1.size(), 9)
END_SECTION

START_SECTION(AASequence getPrefix(Size index) const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("DFP");
  AASequence seq3 = AASequence::fromString("DFPIANGER");
  AASequence seq4 = AASequence::fromString("(TMT6plex)DFPIANGER");
  AASequence seq5 = AASequence::fromString("DFPIANGER(Label:18O(2))");
  AASequence seq6 = AASequence::fromString("DFPIANGERR(Label:18O(2))");
  TEST_EQUAL(seq2, seq1.getPrefix(3));
  TEST_EQUAL(seq3, seq1.getPrefix(9));
  TEST_NOT_EQUAL(seq4.getPrefix(3), seq1.getPrefix(3))
  TEST_NOT_EQUAL(seq5.getPrefix(9), seq1.getPrefix(9))
  TEST_EQUAL(seq6.getPrefix(9), seq1.getPrefix(9))
  TEST_EXCEPTION(Exception::IndexOverflow, seq1.getPrefix(10))
END_SECTION

START_SECTION(AASequence getSuffix(Size index) const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("GER");
  AASequence seq3 = AASequence::fromString("DFPIANGER");
  AASequence seq4 = AASequence::fromString("DFPIANGER(Label:18O(2))");
  AASequence seq5 = AASequence::fromString("(TMT6plex)DFPIANGER");
  AASequence seq6 = AASequence::fromString("(TMT6plex)DDFPIANGER");
  TEST_EQUAL(seq2, seq1.getSuffix(3));
  TEST_EQUAL(seq3, seq1.getSuffix(9));
  TEST_NOT_EQUAL(seq4.getSuffix(3), seq1.getSuffix(3))
  TEST_NOT_EQUAL(seq5.getSuffix(9), seq1.getSuffix(9))
  TEST_EQUAL(seq6.getSuffix(9), seq1.getSuffix(9))
  TEST_EXCEPTION(Exception::IndexOverflow, seq1.getSuffix(10))
END_SECTION

START_SECTION(AASequence getSubsequence(Size index, UInt number) const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("IAN");
  AASequence seq3 = AASequence::fromString("DFPIANGER");
  TEST_EQUAL(seq2, seq1.getSubsequence(3, 3))
  TEST_EQUAL(seq3, seq1.getSubsequence(0, 9))
  TEST_EXCEPTION(Exception::IndexOverflow, seq1.getSubsequence(0, 10))
END_SECTION

START_SECTION(bool has(const Residue& residue) const)
  AASequence seq = AASequence::fromString("DFPIANGER");
  TEST_EQUAL(seq.has(seq[(Size)0]), true)
  Residue res;
  TEST_NOT_EQUAL(seq.has(res), true)
END_SECTION

START_SECTION(bool hasSubsequence(const AASequence& peptide) const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("IANG");
  AASequence seq3 = AASequence::fromString("AIN");
  TEST_EQUAL(seq1.hasSubsequence(seq2), true)
  TEST_EQUAL(seq1.hasSubsequence(seq3), false)
END_SECTION

START_SECTION(bool hasPrefix(const AASequence& peptide) const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("DFP");
  AASequence seq3 = AASequence::fromString("AIN");
  AASequence seq4 = AASequence::fromString("(TMT6plex)DFP");
  AASequence seq5 = AASequence::fromString("DFPIANGER(Label:18O(2))");
  AASequence seq6 = AASequence::fromString("DFP(Label:18O(2))");
  TEST_EQUAL(seq1.hasPrefix(seq2), true)
  TEST_EQUAL(seq1.hasPrefix(seq3), false)
  TEST_EQUAL(seq1.hasPrefix(seq4), false)
  TEST_EQUAL(seq1.hasPrefix(seq5), false)
  TEST_EQUAL(seq1.hasPrefix(seq6), true)
END_SECTION

START_SECTION(bool hasSuffix(const AASequence& peptide) const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("GER");
  AASequence seq3 = AASequence::fromString("AIN");
  AASequence seq4 = AASequence::fromString("GER(Label:18O(2))");
  AASequence seq5 = AASequence::fromString("(TMT6plex)DFPIANGER");
  AASequence seq6 = AASequence::fromString("(TMT6plex)GER");
  TEST_EQUAL(seq1.hasSuffix(seq2), true)
  TEST_EQUAL(seq1.hasSuffix(seq3), false)
  TEST_EQUAL(seq1.hasSuffix(seq4), false)
  TEST_EQUAL(seq1.hasSuffix(seq5), false)
  TEST_EQUAL(seq1.hasSuffix(seq6), true)
END_SECTION

START_SECTION(ConstIterator begin() const)
  String result[] = { "D", "F", "P", "I", "A", "N", "G", "E", "R" };
  AASequence seq = AASequence::fromString("DFPIANGER");
  Size i = 0;
  for (AASequence::ConstIterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it).getOneLetterCode(), result[i])
  }
END_SECTION

START_SECTION(ConstIterator end() const)
  NOT_TESTABLE
END_SECTION

START_SECTION(Iterator begin())
  String result[] = { "D", "F", "P", "I", "A", "N", "G", "E", "R" };
  AASequence seq = AASequence::fromString("DFPIANGER");
  Size i = 0;
  for (AASequence::ConstIterator it = seq.begin(); it != seq.end(); ++it, ++i)
  {
    TEST_EQUAL((*it).getOneLetterCode(), result[i])
  }
END_SECTION

START_SECTION(Iterator end())
  NOT_TESTABLE
END_SECTION

//START_SECTION(friend std::ostream& operator << (std::ostream& os, const AASequence& peptide))
//  // TODO
//END_SECTION

//START_SECTION(friend std::istream& operator > (std::istream& is, const AASequence& peptide))
//  // TODO
//END_SECTION

START_SECTION(String toString() const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("(MOD:00051)DFPIANGER");
  AASequence seq3 = AASequence::fromString("DFPIAN(Deamidated)GER");

  TEST_STRING_EQUAL(seq1.toString(), "DFPIANGER")
  TEST_STRING_EQUAL(seq2.toString(), ".(MOD:00051)DFPIANGER")
  TEST_STRING_EQUAL(seq3.toString(), "DFPIAN(Deamidated)GER")
END_SECTION

START_SECTION(String toUnmodifiedString() const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("(MOD:00051)DFPIANGER");
  AASequence seq3 = AASequence::fromString("DFPIAN(Deamidated)GER");

  TEST_STRING_EQUAL(seq1.toUnmodifiedString(), "DFPIANGER")
  TEST_STRING_EQUAL(seq2.toUnmodifiedString(), "DFPIANGER")
  TEST_STRING_EQUAL(seq3.toUnmodifiedString(), "DFPIANGER")
END_SECTION

START_SECTION(String toUniModString() const)
  AASequence s = AASequence::fromString("PEPC(Carbamidomethyl)PEPM(Oxidation)PEPR");
  TEST_STRING_EQUAL(s.toUniModString(), "PEPC(UniMod:4)PEPM(UniMod:35)PEPR");
  s.setNTerminalModification("Acetyl (N-term)");
  s.setCTerminalModification("Amidated (C-term)");
  TEST_STRING_EQUAL(s.toUniModString(), ".(UniMod:1)PEPC(UniMod:4)PEPM(UniMod:35)PEPR.(UniMod:2)");
END_SECTION

START_SECTION(String toBracketString(const std::vector<String> & fixed_modifications = std::vector<String>()) const)
  AASequence s = AASequence::fromString("PEPC(Carbamidomethyl)PEPM(Oxidation)PEPR");
  TEST_STRING_EQUAL(s.toBracketString(), "PEPC[160]PEPM[147]PEPR");
  vector<String> fixed_mods;
  fixed_mods.push_back("Carbamidomethyl (C)");
  TEST_STRING_EQUAL(s.toBracketString(true, fixed_mods), "PEPCPEPM[147]PEPR");
  TEST_STRING_EQUAL(s.toBracketString(false, fixed_mods), "PEPCPEPM[147.0354000171]PEPR");
  s.setNTerminalModification("Acetyl (N-term)");
  s.setCTerminalModification("Amidated (C-term)");
  TEST_STRING_EQUAL(s.toBracketString(true, fixed_mods), "n[43]PEPCPEPM[147]PEPRc[16]");
  TEST_STRING_EQUAL(s.toBracketString(false, fixed_mods), "n[43.0183900319]PEPCPEPM[147.0354000171]PEPRc[16.0187240319]");
  fixed_mods.push_back("Acetyl (N-term)");
  fixed_mods.push_back("Amidated (C-term)");
  TEST_STRING_EQUAL(s.toBracketString(true, fixed_mods), "PEPCPEPM[147]PEPR");
  TEST_STRING_EQUAL(s.toBracketString(false, fixed_mods), "PEPCPEPM[147.0354000171]PEPR");
END_SECTION

START_SECTION(void setModification(Size index, const String &modification))
  AASequence seq1 = AASequence::fromString("ACDEFNK");
  seq1.setModification(5, "Deamidated");
  TEST_STRING_EQUAL(seq1[5].getModificationName(), "Deamidated");
  // remove modification
  seq1.setModification(5, "");
  TEST_STRING_EQUAL(seq1.toString(), "ACDEFNK") 
END_SECTION

START_SECTION(void setNTerminalModification(const String &modification))
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("(MOD:00051)DFPIANGER");
  TEST_EQUAL(seq1 == seq2, false)
  seq1.setNTerminalModification("MOD:00051");
  TEST_EQUAL(seq1 == seq2, true)

  AASequence seq3 = AASequence::fromString("DABCDEF");
  AASequence seq4 = AASequence::fromString("(MOD:00051)DABCDEF");
  TEST_EQUAL(seq3 == seq4, false)
  seq3.setNTerminalModification("MOD:00051");
  TEST_EQUAL(seq3.isModified(), true)
  TEST_EQUAL(seq4.isModified(), true)
  TEST_EQUAL(seq3 == seq4, true)
END_SECTION

START_SECTION(const String& getNTerminalModificationName() const)
  AASequence seq1 = AASequence::fromString("(MOD:00051)DFPIANGER");
  TEST_EQUAL(seq1.getNTerminalModificationName(), "MOD:00051");

  AASequence seq2 = AASequence::fromString("DFPIANGER");
  TEST_EQUAL(seq2.getNTerminalModificationName(), "");
END_SECTION



START_SECTION(const ResidueModification* getNTerminalModification() const)
  AASequence seq1 = AASequence::fromString("(Formyl)DFPIANGER");
  TEST_EQUAL(seq1.getNTerminalModification()->getId(), "Formyl");
  TEST_EQUAL(seq1.getNTerminalModification()->getFullId(), "Formyl (N-term)");

  AASequence seq2 = AASequence::fromString("DFPIANGER");
  TEST_EQUAL(seq2.getNTerminalModification(),  0);
END_SECTION

START_SECTION(const ResidueModification* getCTerminalModification() const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("DFPIANGER(Amidated)");

  TEST_EQUAL(seq2.getCTerminalModification()->getId(), "Amidated");
  TEST_EQUAL(seq2.getCTerminalModification()->getFullId(), "Amidated (C-term)");
  TEST_EQUAL(seq1.getCTerminalModification(),  0);
END_SECTION

START_SECTION(void setCTerminalModification(const String& modification))
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("DFPIANGER(Amidated)");

  TEST_EQUAL(seq1 == seq2, false)
  seq1.setCTerminalModification("Amidated");
  TEST_EQUAL(seq1 == seq2, true)

  AASequence seq3 = AASequence::fromString("DABCDER");
  AASequence seq4 = AASequence::fromString("DABCDER(Amidated)");
  TEST_EQUAL(seq3 == seq4, false)
  seq3.setCTerminalModification("Amidated");
  TEST_EQUAL(seq3.isModified(), true)
  TEST_EQUAL(seq4.isModified(), true)
  TEST_EQUAL(seq3 == seq4, true)

  AASequence seq5 = AASequence::fromString("DABCDER(MOD:00177)");
  AASequence seq6 = AASequence::fromString("DABCDER(MOD:00177)(Amidated)");
  TEST_EQUAL(seq5.isModified(), true)
  TEST_EQUAL(seq6.isModified(), true)
  seq5.setCTerminalModification("Amidated");
  TEST_EQUAL(seq5 == seq6, true)

  AASequence seq7 = AASequence::fromString("DFPIANGER(MOD:00177)");
  AASequence seq8 = AASequence::fromString("DFPIANGER(MOD:00177)(Amidated)");
  TEST_EQUAL(seq7.isModified(), true)
  TEST_EQUAL(seq8.isModified(), true)
  seq7.setCTerminalModification("Amidated");
  TEST_EQUAL(seq5 == seq6, true)
END_SECTION

START_SECTION(const String& getCTerminalModificationName() const)
  AASequence seq1 = AASequence::fromString("DFPIANGER(Amidated)");
  TEST_EQUAL(seq1.getCTerminalModificationName(), "Amidated");

  AASequence seq2 = AASequence::fromString("DFPIANGER");
  TEST_EQUAL(seq2.getCTerminalModificationName(), "");
END_SECTION

START_SECTION(bool hasNTerminalModification() const)
  AASequence seq1 = AASequence::fromString("(MOD:00051)DABCDEF");
  AASequence seq2 = AASequence::fromString("DABCDEF");

  TEST_EQUAL(seq1.hasNTerminalModification(), true)
  TEST_EQUAL(seq2.hasNTerminalModification(), false)

  AASequence seq3 = AASequence::fromString("(MOD:00051)DFPIANGER");
  AASequence seq4 = AASequence::fromString("DFPIANGER");
  TEST_EQUAL(seq3.hasNTerminalModification(), true)
  TEST_EQUAL(seq4.hasNTerminalModification(), false)
END_SECTION

START_SECTION(bool hasCTerminalModification() const)
  AASequence seq1 = AASequence::fromString("DFPIANGER(Amidated)");
  AASequence seq2 = AASequence::fromString("DFPIANGER");
  TEST_EQUAL(seq1.hasCTerminalModification(), true)
  TEST_EQUAL(seq2.hasCTerminalModification(), false)
  seq1.setCTerminalModification("");
  TEST_EQUAL(seq1.hasCTerminalModification(), false)
END_SECTION

START_SECTION(bool isModified() const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  TEST_EQUAL(seq1.isModified(), false);
  AASequence seq2(seq1);
  seq2.setNTerminalModification("MOD:09999");
  TEST_EQUAL(seq2.isModified(), true)

  AASequence seq3(seq1);
  seq3.setCTerminalModification("Amidated");
  TEST_EQUAL(seq3.isModified(), true);

  AASequence seq4 = AASequence::fromString("DFPIANGER(MOD:00177)");
  TEST_EQUAL(seq4.isModified(), true);
END_SECTION

START_SECTION(bool operator<(const AASequence &rhs) const)
  AASequence seq1 = AASequence::fromString("DFPIANGER");
  AASequence seq2 = AASequence::fromString("DFBIANGER");
  TEST_EQUAL(seq2 < seq1, true)
  TEST_EQUAL(seq1 < seq2, false)
  AASequence seq3 = AASequence::fromString("DFPIANGFR");
  TEST_EQUAL(seq3 < seq1, false)

  // shorter residue sequence is smaller than longer one
  TEST_EQUAL(AASequence::fromString("PPP") < AASequence::fromString("AAAA"), true)
  TEST_EQUAL(AASequence::fromString("PM(Oxidation)P") < AASequence::fromString("AAAA"), true)

  // modified is larger than unmodified
  TEST_EQUAL(AASequence::fromString("MMM") < AASequence::fromString("MM(Oxidation)M"), true)
  TEST_EQUAL(AASequence::fromString("ARRR") < AASequence::fromString("ARRR(Label:13C(6))"), true)
  TEST_EQUAL(AASequence::fromString("CNR") < AASequence::fromString("(ICPL:2H(4))CNR"), true)
  TEST_EQUAL(AASequence::fromString("(ICPL:2H(4))CNAR") < AASequence::fromString("(ICPL:13C(6))YCYCY"), true)

  // alphabetic order
  TEST_EQUAL(AASequence::fromString("AAA") < AASequence::fromString("AAM"), true)
  TEST_EQUAL(AASequence::fromString("AAM") < AASequence::fromString("AMA"), true)
  TEST_EQUAL(AASequence::fromString("AMA") < AASequence::fromString("MAA"), true)

  // if N-terminal mods. are the same, check the sequence
  TEST_EQUAL(AASequence::fromString("(ICPL:2H(4))AMA") < AASequence::fromString("(ICPL:2H(4))MAA"), true)
  TEST_EQUAL(AASequence::fromString("(ICPL:2H(4))MAA") < AASequence::fromString("(ICPL:2H(4))AMA"), false)
  // if everything else is the same, check the C-terminal mods.
  TEST_EQUAL(AASequence::fromString("(ICPL:2H(4))AMA(Amidated)") < AASequence::fromString("(ICPL:2H(4))AMA(Label:18O(2))"), true)
  TEST_EQUAL(AASequence::fromString("(ICPL:2H(4))AMA(Label:18O(2))") < AASequence::fromString("(ICPL:2H(4))AMA(Amidated)"), false)

END_SECTION

START_SECTION(bool operator!=(const AASequence& rhs) const)
  AASequence seq1 = AASequence::fromString("(MOD:00051)DFPIANGER");
  AASequence seq2 = AASequence::fromString("DFPIANGER");
  TEST_EQUAL(seq2 != AASequence::fromString("DFPIANGER"), false)
  TEST_EQUAL(seq1 != AASequence::fromString("(MOD:00051)DFPIANGER"), false)

  // test C-terminal mods
  AASequence seq3 = AASequence::fromString("DFPIANGER(MOD:00177)");
  AASequence seq4 = AASequence::fromString("DFPIANGER(Amidated)");
  TEST_EQUAL(seq3 != AASequence::fromString("DFPIANGER"), true)
  TEST_EQUAL(seq3 != AASequence::fromString("DFPIANGER(MOD:00177)"), false)
  TEST_EQUAL(seq4 != AASequence::fromString("DFPIANGER(Amidated)"), false)
  TEST_EQUAL(seq4 != AASequence::fromString("DFPIANGER"), true)

  // test inner mods
  TEST_EQUAL(AASequence::fromString("DFPMIANGER") != AASequence::fromString("DFPM(Oxidation)IANGER"), true)
  TEST_EQUAL(AASequence::fromString("DFPM(Oxidation)IANGER") == AASequence::fromString("DFPM(Oxidation)IANGER"), true)

  AASequence seq5 = AASequence::fromString("DFBIANGER");
  TEST_EQUAL(seq5 != AASequence::fromString("DFPIANGER"), true)
  TEST_EQUAL(seq5 != AASequence::fromString("DFBIANGER"), false)
END_SECTION

START_SECTION(void getAAFrequencies(Map<String, Size>& frequency_table) const)
  AASequence a = AASequence::fromString("THREEAAAWITHYYY");
  Map<String, Size> table;
  a.getAAFrequencies(table);

  TEST_EQUAL(table["T"]==2, true);
  TEST_EQUAL(table["H"]==2, true);
  TEST_EQUAL(table["R"]==1, true);
  TEST_EQUAL(table["E"]==2, true);
  TEST_EQUAL(table["A"]==3, true);
  TEST_EQUAL(table["W"]==1, true);
  TEST_EQUAL(table["I"]==1, true);
  TEST_EQUAL(table["Y"]==3, true);

  TEST_EQUAL(table.size()==8, true);

END_SECTION

START_SECTION([EXTRA] Tag in peptides)
{
  AASequence aa1 = AASequence::fromString("PEPTC[+57.02]IDE"); // 57.021464
  AASequence aa2 = AASequence::fromString("PEPTC(Carbamidomethyl)IDE");
  AASequence aa3 = AASequence::fromString("PEPTC(UniMod:4)IDE");
  AASequence aa4 = AASequence::fromString("PEPTC(Iodoacetamide derivative)IDE");
  AASequence aa5 = AASequence::fromString("PEPTC[160.030654]IDE");
  AASequence aa6 = AASequence::fromString("PEPTX[160.030654]IDE");

  TEST_REAL_SIMILAR(aa1.getMonoWeight(), 959.39066)
  TEST_REAL_SIMILAR(aa2.getMonoWeight(), 959.39066)
  TEST_REAL_SIMILAR(aa3.getMonoWeight(), 959.39066)
  TEST_REAL_SIMILAR(aa4.getMonoWeight(), 959.39066)
  TEST_REAL_SIMILAR(aa5.getMonoWeight(), 959.39066)
  TEST_REAL_SIMILAR(aa6.getMonoWeight(), 959.39066)

  TEST_EQUAL(aa1.size(), 8)
  TEST_EQUAL(aa2.size(), 8)
  TEST_EQUAL(aa3.size(), 8)
  TEST_EQUAL(aa4.size(), 8)
  TEST_EQUAL(aa5.size(), 8)
  TEST_EQUAL(aa6.size(), 8)

  TEST_EQUAL(aa1.isModified(), true)
  TEST_EQUAL(aa2.isModified(), true)
  TEST_EQUAL(aa3.isModified(), true)
  TEST_EQUAL(aa4.isModified(), true)
  TEST_EQUAL(aa5.isModified(), true)
  TEST_EQUAL(aa6.isModified(), true)

  // Test negative mods / losses
  // test without loss
  TEST_REAL_SIMILAR(AASequence::fromString("PEPTMIDE").getMonoWeight(), 930.4004)
  // test with losses
  // known loss from unimod: Homoserine (should actually only happen at c-term but we allow it)
  TEST_REAL_SIMILAR(AASequence::fromString("PEPTM[-30]IDE").getMonoWeight(), 930.4004 - 29.992806)
  // new loss from unimod: Homoserine (should actually only happen at c-term but we allow it)
  TEST_REAL_SIMILAR(AASequence::fromString("PEPTM[-30.4004]IDE").getMonoWeight(), 900.0)
  TEST_EQUAL(AASequence::fromString("PEPTM[-30]IDE").size(), 8)
  TEST_EQUAL(AASequence::fromString("PEPTM[-30]IDE").isModified(), true)
}
END_SECTION

START_SECTION([EXTRA] Arbitrary tag in peptides using square brackets)
{
  // test arbitrary modification
  AASequence aa_original = AASequence::fromString("PEPTIDE");
  TEST_REAL_SIMILAR(aa_original.getMonoWeight(), 799.36001)

  AASequence aa_half = AASequence::fromString("IDE");
  TEST_REAL_SIMILAR(aa_half.getMonoWeight(), 375.1641677975)

  AASequence aa = AASequence::fromString("PEPTX[999]IDE");
  TEST_REAL_SIMILAR(aa.getMonoWeight(), 799.36001 + 999.0)
  TEST_REAL_SIMILAR(aa.getMonoWeight(), aa_original.getMonoWeight() + 999.0)

  {
    AASequence aa = AASequence::fromString("PEPTX[999.0]IDE");
    TEST_REAL_SIMILAR(aa.getMonoWeight(), 799.36001 + 999.0)
    TEST_REAL_SIMILAR(aa.getMonoWeight(), aa_original.getMonoWeight() + 999.0)
  }

  // test arbitrary differences (e.g. it should be possible to encode arbitrary masses and still get the correct weight)
  {
    AASequence test1 = AASequence::fromString("PEPTX[160.030654]IDE");
    TEST_REAL_SIMILAR(test1.getMonoWeight(), aa_original.getMonoWeight() + 160.030654)
    AASequence test2 = AASequence::fromString("PEPTX[160.040654]IDE");
    TEST_REAL_SIMILAR(test2.getMonoWeight(), aa_original.getMonoWeight() + 160.040654)
    AASequence test3 = AASequence::fromString("PEPTX[160.050654]IDE");
    TEST_REAL_SIMILAR(test3.getMonoWeight(), aa_original.getMonoWeight() + 160.050654)
    AASequence test4 = AASequence::fromString("PEPTX[160.130654]IDE");
    TEST_REAL_SIMILAR(test4.getMonoWeight(), aa_original.getMonoWeight() + 160.130654)
    AASequence test5 = AASequence::fromString("PEPTX[160.230654]IDE");
    TEST_REAL_SIMILAR(test5.getMonoWeight(), aa_original.getMonoWeight() + 160.230654)
  }

  // test arbitrary differences (e.g. it should be possible to encode arbitrary masses and still get the correct weight)
  {
    AASequence test1 = AASequence::fromString("PEPTN[160.030654]IDE");
    TEST_REAL_SIMILAR(test1.getMonoWeight(), aa_original.getMonoWeight() + 160.030654)
    AASequence test2 = AASequence::fromString("PEPTN[160.040654]IDE");
    TEST_REAL_SIMILAR(test2.getMonoWeight(), aa_original.getMonoWeight() + 160.040654)
    AASequence test3 = AASequence::fromString("PEPTN[160.050654]IDE");
    TEST_REAL_SIMILAR(test3.getMonoWeight(), aa_original.getMonoWeight() + 160.050654)
    AASequence test4 = AASequence::fromString("PEPTN[160.130654]IDE");
    TEST_REAL_SIMILAR(test4.getMonoWeight(), aa_original.getMonoWeight() + 160.130654)
    AASequence test5 = AASequence::fromString("PEPTN[160.230654]IDE");
    TEST_REAL_SIMILAR(test5.getMonoWeight(), aa_original.getMonoWeight() + 160.230654)
  }

  // test arbitrary differences (e.g. it should be possible to encode arbitrary masses and still get the correct weight)
  {
    AASequence test1 = AASequence::fromString("PEPT[+160.030654]IDE");
    TEST_REAL_SIMILAR(test1.getMonoWeight(), aa_original.getMonoWeight() + 160.030654)

    TEST_REAL_SIMILAR(test1[0].getMonoWeight(), 115.0633292871)
    TEST_EQUAL(test1[3].isModified(), true)
    TEST_REAL_SIMILAR(test1[3].getMonoWeight(), 160.030654 + 119.0582442871) // the weight at this position is the mod + residue itself

    AASequence test2 = AASequence::fromString("PEPT[+160.040654]IDE");
    TEST_REAL_SIMILAR(test2.getMonoWeight(), aa_original.getMonoWeight() + 160.040654)
    AASequence test3 = AASequence::fromString("PEPT[+160.050654]IDE");
    TEST_REAL_SIMILAR(test3.getMonoWeight(), aa_original.getMonoWeight() + 160.050654)
    AASequence test4 = AASequence::fromString("PEPT[+160.130654]IDE");
    TEST_REAL_SIMILAR(test4.getMonoWeight(), aa_original.getMonoWeight() + 160.130654)
    AASequence test5 = AASequence::fromString("PEPT[+160.230654]IDE");
    TEST_REAL_SIMILAR(test5.getMonoWeight(), aa_original.getMonoWeight() + 160.230654)
  }

  // test arbitrary differences when writing them out
  AASequence test6 = AASequence::fromString("PEPTX[1600.230654]IDE");
  TEST_EQUAL(test6.size(), 8)
  TEST_EQUAL(test6.toString(), "PEPTX[1600.230654]IDE")
  TEST_EQUAL(test6.toUniModString(), "PEPTX[1600.230654]IDE")
  TEST_EQUAL(test6.toBracketString(), "PEPTX[1600]IDE")
  TEST_EQUAL(test6.toBracketString(false), "PEPTX[1600.230654]IDE")
  TEST_EQUAL(test6.toUnmodifiedString(), "PEPTXIDE")

  TEST_EQUAL(test6[4].isModified(), true)
  TEST_REAL_SIMILAR(test6[4].getModification()->getMonoMass(), 1600.230654 + 18.0105650638) // because of the H2O loss
  TEST_REAL_SIMILAR(test6[4].getMonoWeight(), 1600.230654 + 18.0105650638) // computed as full (not internal) -> H2O loss
  TEST_REAL_SIMILAR(test6[4].getMonoWeight(Residue::Internal), 1600.230654)
  TEST_REAL_SIMILAR(test6.getMonoWeight(), aa_original.getMonoWeight() + 1600.230654)

  AASequence test7 = AASequence::fromString(test6.toString());
  TEST_EQUAL(test7.size(), 8)
  TEST_EQUAL(test7.toString(), "PEPTX[1600.230654]IDE")
  TEST_EQUAL(test7.toUniModString(), "PEPTX[1600.230654]IDE")
  TEST_EQUAL(test7.toBracketString(), "PEPTX[1600]IDE")
  TEST_EQUAL(test7.toBracketString(false), "PEPTX[1600.230654]IDE")
  TEST_EQUAL(test7.toUnmodifiedString(), "PEPTXIDE")
  TEST_EQUAL(test6, test7) // the peptides should be equal

  // test arbitrary modification on N
  {
    AASequence test_seq = AASequence::fromString("PEPTN[1600.230654]IDE");
    TEST_EQUAL(test_seq.size(), 8)
    TEST_EQUAL(test_seq.toString(), "PEPTN[1600.230654]IDE")
    TEST_EQUAL(test_seq.toUniModString(), "PEPTN[1600.230654]IDE")
    TEST_EQUAL(test_seq.toBracketString(false), "PEPTN[1600.230654]IDE")
    TEST_EQUAL(test_seq.toUnmodifiedString(), "PEPTNIDE")
    TEST_REAL_SIMILAR(test_seq.getMonoWeight(), aa_original.getMonoWeight() + 1600.230654)

    // test that we can re-read the string
    AASequence test_other = AASequence::fromString(test_seq.toString());
    TEST_EQUAL(test_other.size(), 8)
    TEST_EQUAL(test_other.toString(), "PEPTN[1600.230654]IDE")
    TEST_EQUAL(test_other.toUniModString(), "PEPTN[1600.230654]IDE")
    TEST_EQUAL(test_other.toBracketString(false), "PEPTN[1600.230654]IDE")
    TEST_EQUAL(test_other.toUnmodifiedString(), "PEPTNIDE")

    TEST_EQUAL(test_other, test_seq) // the peptides should be equal

    // test that we can re-read the string from BracketString
    test_other = AASequence::fromString(test_seq.toBracketString(false));
    TEST_EQUAL(test_other.size(), 8)
    TEST_EQUAL(test_other.toString(), "PEPTN[1600.230654]IDE")
    TEST_EQUAL(test_other.toUniModString(), "PEPTN[1600.230654]IDE")
    TEST_EQUAL(test_other.toBracketString(false), "PEPTN[1600.230654]IDE")
    TEST_EQUAL(test_other.toUnmodifiedString(), "PEPTNIDE")

    TEST_EQUAL(test_other, test_seq) // the peptides should be equal
  }

  // test N-terminal modification
  {
    AASequence test_seq = AASequence::fromString("[1600.230654]IDE");
    TEST_EQUAL(test_seq.size(), 3)
    TEST_EQUAL(test_seq.toString(), ".[1600.230654]IDE")
    TEST_EQUAL(test_seq.toUniModString(), ".[1600.230654]IDE")
    TEST_EQUAL(test_seq.toBracketString(false), "n[1600.230654]IDE")
    TEST_EQUAL(test_seq.toUnmodifiedString(), "IDE")
    TEST_REAL_SIMILAR(test_seq.getMonoWeight(), aa_half.getMonoWeight() + 1600.230654)

    // test that we can re-read the string
    AASequence test_other = AASequence::fromString(test_seq.toString());
    TEST_EQUAL(test_other.size(), 3)
    TEST_EQUAL(test_other.toString(), ".[1600.230654]IDE")

    TEST_EQUAL(test_seq, test_other) // the peptides should be equal

    // test that we can re-read the string from UniModString
    test_other = AASequence::fromString(test_seq.toUniModString());
    TEST_EQUAL(test_other.size(), 3)
    TEST_EQUAL(test_other.toString(), ".[1600.230654]IDE")

    TEST_EQUAL(test_seq, test_other) // the peptides should be equal

    // test that we can re-read the string from BracketString
    test_other = AASequence::fromString(test_seq.toBracketString(false));
    TEST_EQUAL(test_other.size(), 3)
    TEST_EQUAL(test_other.toString(), ".[1600.230654]IDE")

    TEST_EQUAL(test_seq, test_other) // the peptides should be equal
  }

  // test C-terminal modification
  {
    AASequence test_seq = AASequence::fromString("IDE.[1600.230654]");
    TEST_EQUAL(test_seq.size(), 3)
    TEST_EQUAL(test_seq.toString(), "IDE.[1600.230654]")
    TEST_EQUAL(test_seq.toUniModString(), "IDE.[1600.230654]")
    TEST_EQUAL(test_seq.toBracketString(false), "IDEc[1600.230654]")
    TEST_EQUAL(test_seq.toUnmodifiedString(), "IDE")
    TEST_REAL_SIMILAR(test_seq.getMonoWeight(), aa_half.getMonoWeight() + 1600.230654)

    // test that we can re-read the string
    AASequence test_other = AASequence::fromString(test_seq.toString());
    TEST_EQUAL(test_other.size(), 3)
    TEST_EQUAL(test_other.toString(), "IDE.[1600.230654]")

    TEST_EQUAL(test_seq, test_other) // the peptides should be equal

    // test that we can re-read the UniModString
    test_other = AASequence::fromString(test_seq.toUniModString());
    TEST_EQUAL(test_other.size(), 3)
    TEST_EQUAL(test_other.toString(), "IDE.[1600.230654]")

    TEST_EQUAL(test_seq, test_other) // the peptides should be equal

    // test that we can re-read the string from BracketString
    test_other = AASequence::fromString(test_seq.toBracketString(false));
    TEST_EQUAL(test_other.size(), 3)
    TEST_EQUAL(test_other.toString(), "IDE.[1600.230654]")

    TEST_EQUAL(test_seq, test_other) // the peptides should be equal
  }

  // test from #2320
  {
    AASequence aa_original = AASequence::fromString("DFPANGER");
    TEST_REAL_SIMILAR(aa_original.getMonoWeight(), 904.4038997864)

    AASequence aa_withX_original = AASequence::fromString("DFPANGERX");
    // TEST_REAL_SIMILAR(aa_withX_original.getMonoWeight(), 904.4038997864) // cannot call weight on AASequence with X

    AASequence test_seq = AASequence::fromString("DFPANGERX[113.0840643509]");
    TEST_EQUAL(test_seq.size(), 9)
    TEST_EQUAL(test_seq.toString(), "DFPANGERX[113.0840643509]")
    TEST_EQUAL(test_seq.toUniModString(), "DFPANGERX[113.0840643509]")
    TEST_EQUAL(test_seq.toBracketString(true), "DFPANGERX[113]")
    TEST_EQUAL(test_seq.toBracketString(false), "DFPANGERX[113.0840643509]")
    TEST_EQUAL(test_seq.toUnmodifiedString(), "DFPANGERX")
    TEST_EQUAL(test_seq.hasNTerminalModification(), false)
    TEST_EQUAL(test_seq.hasCTerminalModification(), false)
    TEST_REAL_SIMILAR(test_seq.getMonoWeight(), aa_original.getMonoWeight() + 113.0840643509)

    AASequence test_internal_seq = AASequence::fromString("DFPANGEX[113.0840643509]R");
    TEST_REAL_SIMILAR(test_internal_seq.getMonoWeight(), aa_original.getMonoWeight() + 113.0840643509)

    TEST_EQUAL(test_seq[8].isModified(), true)
    TEST_REAL_SIMILAR(test_seq[8].getMonoWeight(), 113.0840643509 + 18.0105650638) // computed as full (not internal) -> H2O loss
    TEST_REAL_SIMILAR(test_seq[8].getMonoWeight(Residue::Internal), 113.0840643509)

    // test that we can re-read the string
    AASequence test_other = AASequence::fromString(test_seq.toString());
    TEST_EQUAL(test_seq.size(), 9)
    TEST_EQUAL(test_seq.toString(), "DFPANGERX[113.0840643509]")

    TEST_EQUAL(test_seq, test_other) // the peptides should be equal
  }

  // Faulty / nonsense calculations ...
  AASequence test;
  TEST_EXCEPTION(Exception::ParseError, test = AASequence::fromString("PEPTX[+160.230654]IDE"));

  AASequence seq11 = AASequence::fromString("PEPM[147.035405]TIDEK");
  TEST_EQUAL(seq11.isModified(), true);
  TEST_STRING_EQUAL(seq11[3].getModificationName(), "Oxidation");
}
END_SECTION

START_SECTION([EXTRA] Test integer vs float tags)
{
  /// Test absolute masses

  // Test a few modifications with the "correct" accurate mass
  {
  AASequence seq11 = AASequence::fromString("PEPM[147.035405]TIDEK"); // UniMod oxMet is 147.035405
  TEST_EQUAL(seq11.isModified(), true);
  TEST_STRING_EQUAL(seq11[3].getModificationName(), "Oxidation");
  TEST_EQUAL(seq11[3].getModification()->getUniModRecordId(), 35)
  TEST_EQUAL(seq11[3].getModification()->getUniModAccession(), "UniMod:35")

  AASequence seq12 = AASequence::fromString("PEPT[181.014]TIDEK");
  TEST_EQUAL(seq12.isModified(), true);
  TEST_STRING_EQUAL(seq12[3].getModificationName(), "Phospho");
  TEST_EQUAL(seq12[3].getModification()->getUniModRecordId(), 21)
  TEST_EQUAL(seq12[3].getModification()->getUniModAccession(), "UniMod:21")

  AASequence seq13 = AASequence::fromString("PEPY[243.03]TIDEK");
  TEST_EQUAL(seq13.isModified(), true);
  TEST_STRING_EQUAL(seq13[3].getModificationName(), "Phospho");
  TEST_EQUAL(seq13[3].getModification()->getUniModRecordId(), 21)
  TEST_EQUAL(seq13[3].getModification()->getUniModAccession(), "UniMod:21")

  AASequence seq15 = AASequence::fromString("PEPC[160.0306]TIDE");
  TEST_EQUAL(seq15.isModified(), true);
  TEST_STRING_EQUAL(seq15[3].getModificationName(), "Carbamidomethyl");
  TEST_EQUAL(seq15[3].getModification()->getUniModRecordId(), 4)
  TEST_EQUAL(seq15[3].getModification()->getUniModAccession(), "UniMod:4")
  }

  // Test a few modifications with the accurate mass slightly off to match some other modification
  {
  AASequence seq11 = AASequence::fromString("PEPM[147.01]TIDEK");
  TEST_EQUAL(seq11.isModified(), true);
  TEST_STRING_EQUAL(seq11[3].getModificationName(), "Oxidation")
  TEST_EQUAL(seq11[3].getModification()->getUniModRecordId(), 35)

  AASequence seq12 = AASequence::fromString("PEPT[181.004]TIDEK");
  TEST_EQUAL(seq12.isModified(), true);
  TEST_STRING_EQUAL(seq12[3].getModificationName(), "Sulfo");
  TEST_EQUAL(seq12[3].getModification()->getUniModRecordId(), 40)

  AASequence seq13 = AASequence::fromString("PEPY[243.02]TIDEK");
  TEST_EQUAL(seq13.isModified(), true);
  TEST_STRING_EQUAL(seq13[3].getModificationName(), "Sulfo");
  TEST_EQUAL(seq13[3].getModification()->getUniModRecordId(), 40)

  AASequence seq14 = AASequence::fromString("PEPTC[159.035405]IDE");
  TEST_EQUAL(seq14.isModified(), true);
  TEST_STRING_EQUAL(seq14[4].getModificationName(), "Delta:H(4)C(3)O(1)");
  TEST_EQUAL(seq14[4].getModification()->getUniModRecordId(), 206)
  }


  /// Test delta masses

  // Test a few modifications with the "correct" accurate mass
  {
  AASequence seq11 = AASequence::fromString("PEPM[+15.994915]TIDEK"); // UniMod oxMet is 15.994915
  TEST_EQUAL(seq11.isModified(), true);
  TEST_STRING_EQUAL(seq11[3].getModificationName(), "Oxidation");

  AASequence seq12 = AASequence::fromString("PEPT[+79.96632]TIDEK");
  TEST_EQUAL(seq12.isModified(), true);
  TEST_STRING_EQUAL(seq12[3].getModificationName(), "Phospho");

  AASequence seq13 = AASequence::fromString("PEPY[+79.966331]TIDEK");
  TEST_EQUAL(seq13.isModified(), true);
  TEST_STRING_EQUAL(seq13[3].getModificationName(), "Phospho");

  AASequence seq14 = AASequence::fromString("PEPC[+57.02]TIDE");
  TEST_EQUAL(seq14.isModified(), true);
  TEST_STRING_EQUAL(seq14[3].getModificationName(), "Carbamidomethyl");
  }

  // Test a few modifications with the accurate mass slightly off to match some other modification
  {
  // this does not work any more since there is no difference in the oxygen atom
  // AASequence seq11("PEPM[+15.994909]TIDEK"); // PSI-MOD oxMet is 15.994909
  // TEST_EQUAL(seq11.isModified(),true);
  // TEST_STRING_EQUAL(seq11[3].getModification(), "MOD:00719")

  AASequence seq12 = AASequence::fromString("PEPT[+79.957]TIDEK");
  TEST_EQUAL(seq12.isModified(), true);
  TEST_STRING_EQUAL(seq12[3].getModificationName(), "Sulfo")

  AASequence seq13 = AASequence::fromString("PEPY[+79.9568]TIDEK");
  TEST_EQUAL(seq13.isModified(), true);
  TEST_STRING_EQUAL(seq13[3].getModificationName(), "Sulfo")

  AASequence seq14 = AASequence::fromString("PEPTC[+56.026215]IDE");
  TEST_EQUAL(seq14.isModified(), true);
  TEST_STRING_EQUAL(seq14[4].getModificationName(), "Delta:H(4)C(3)O(1)");
  }

}
END_SECTION

START_SECTION([EXTRA] Peptide equivalence)
{

  // Test float mass tag leading to internal Acetylation
  TEST_EQUAL(AASequence::fromString("PEPTC(Acetyl)IDE"), AASequence::fromString("PEPTC[+42.011]IDE"))

  // Test mass tag leading to a "Dimethyl:2H(4)13C(2)" at the N-term (pepXML format)
  TEST_EQUAL(AASequence::fromString("(Dimethyl:2H(4)13C(2))TGSESSQTGTSTTSSR"), AASequence::fromString("n[+34]TGSESSQTGTSTTSSR"))

  // Test absolute mass (34 + 1H) leading to a "Dimethyl:2H(4)13C(2)" at the N-term (pepXML format)
  TEST_EQUAL(AASequence::fromString("(Dimethyl:2H(4)13C(2))TGSESSQTGTSTTSSR"), AASequence::fromString("n[35]TGSESSQTGTSTTSSR"))

  // Test absolute integer mass at the N-term, internal residue and C-terminus (mixing pepXML bracket format and OpenMS round brackets)
  TEST_EQUAL(AASequence::fromString("n(Acetyl)PEPC(Carbamidomethyl)PEPM(Oxidation)PEPRc(Amidated)"), AASequence::fromString("n[43]PEPC(Carbamidomethyl)PEPM[147]PEPRc[16]"));

  // Test float mass tag leading to N-terminal Acetylation
  TEST_EQUAL(AASequence::fromString("(Acetyl)PEPTCIDE"), AASequence::fromString("[+42.011]PEPTCIDE"))

  // Test integer mass tag leading to N-terminal Acetylation
  TEST_EQUAL(AASequence::fromString("(Acetyl)PEPTCIDE"), AASequence::fromString("[+42]PEPTCIDE"))

  // Test Carbamidomethyl
  TEST_EQUAL(AASequence::fromString("PEPTC(UniMod:4)IDE"), AASequence::fromString("PEPTC(Carbamidomethyl)IDE"))
  TEST_EQUAL(AASequence::fromString("PEPTC(UniMod:4)IDE"), AASequence::fromString("PEPTC(Iodoacetamide derivative)IDE"))
  TEST_EQUAL(AASequence::fromString("PEPTC(UniMod:4)IDE"), AASequence::fromString("PEPTC[160.030654]IDE")) // 103.00919 + 57.02
  TEST_EQUAL(AASequence::fromString("PEPTC(UniMod:4)IDE"), AASequence::fromString("PEPTC[+57.02]IDE"))
  TEST_EQUAL(AASequence::fromString("PEPTC(UniMod:4)IDE"), AASequence::fromString("PEPTC[160]IDE")) // 103.00919 + 57.02
  TEST_EQUAL(AASequence::fromString("PEPTC(UniMod:4)IDE"), AASequence::fromString("PEPTC[+57]IDE"))

  // Test Oxidation
  TEST_EQUAL(AASequence::fromString("DFPIAM(UniMod:35)GER"), AASequence::fromString("DFPIAM[+16]GER"))
  TEST_EQUAL(AASequence::fromString("DFPIAM(UniMod:35)GER"), AASequence::fromString("DFPIAM[147]GER"))
  TEST_EQUAL(AASequence::fromString("DFPIAM(UniMod:35)GER"), AASequence::fromString("DFPIAM[+15.99]GER"))
  TEST_EQUAL(AASequence::fromString("DFPIAM(UniMod:35)GER"), AASequence::fromString("DFPIAM[147.035405]GER"))
  TEST_EQUAL(AASequence::fromString("DFPIAM(UniMod:35)GER"), AASequence::fromString("DFPIAM(Oxidation)GER"))

  // Test Phosphorylation
  TEST_EQUAL(AASequence::fromString("PEPT(UniMod:21)TIDEK"), AASequence::fromString("PEPT(Phospho)TIDEK"))
  TEST_EQUAL(AASequence::fromString("PEPT(UniMod:21)TIDEK"), AASequence::fromString("PEPT[181]TIDEK"))
  TEST_EQUAL(AASequence::fromString("PEPT(UniMod:21)TIDEK"), AASequence::fromString("PEPT[+80]TIDEK"))

  TEST_EQUAL(AASequence::fromString("PEPY(UniMod:21)TIDEK"), AASequence::fromString("PEPY(Phospho)TIDEK"))
  TEST_EQUAL(AASequence::fromString("PEPY(UniMod:21)TIDEK"), AASequence::fromString("PEPY[243]TIDEK"))
  TEST_EQUAL(AASequence::fromString("PEPY(UniMod:21)TIDEK"), AASequence::fromString("PEPY[+80]TIDEK"))

  TEST_EQUAL(AASequence::fromString("PEPS(UniMod:21)TIDEK"), AASequence::fromString("PEPS(Phospho)TIDEK"))
  TEST_EQUAL(AASequence::fromString("PEPS(UniMod:21)TIDEK"), AASequence::fromString("PEPS[167]TIDEK"))
  TEST_EQUAL(AASequence::fromString("PEPS(UniMod:21)TIDEK"), AASequence::fromString("PEPS[+80]TIDEK"))

  // Test loss
  TEST_EQUAL(AASequence::fromString("PEPTIDEM(UniMod:10)"), AASequence::fromString("PEPTIDEM(Met->Hse)"))

  TEST_EQUAL(AASequence::fromString("PEPTM(Met->Thr)IDE"), AASequence::fromString("PEPTM[-30]IDE"))
  TEST_EQUAL(AASequence::fromString("PEPTM(Met->Thr)IDE"), AASequence::fromString("PEPTM[101]IDE"))
}
END_SECTION

START_SECTION([EXTRA] Tag in peptides)
{
  String I_weight = String(ResidueDB::getInstance()->getResidue("I")->getMonoWeight(Residue::Internal));
  AASequence aa1 = AASequence::fromString("DFPIANGER");
  AASequence aa2 = AASequence::fromString("DPFX[" + I_weight + "]ANGER");
  AASequence aa3 = AASequence::fromString("X[" + I_weight + "]DFPANGER");
  AASequence aa4 = AASequence::fromString("DFPANGERX[" + I_weight + "]");
  TEST_REAL_SIMILAR(aa1.getMonoWeight(), 1017.487958568)
  TEST_EQUAL(aa2.isModified(), true)
  TEST_EQUAL(aa3.hasNTerminalModification(), false)
  TEST_EQUAL(aa4.hasCTerminalModification(), false)
  TEST_REAL_SIMILAR(aa2.getMonoWeight(), 1017.487958568)
  TEST_REAL_SIMILAR(aa3.getMonoWeight(), 1017.487958568)
  TEST_REAL_SIMILAR(aa4.getMonoWeight(), 1017.487958568)
}
END_SECTION

START_SECTION([EXTRA] testing terminal modifications)
{
  AASequence aaNoMod = AASequence::fromString(".DFPIANGER.");
  AASequence aaNtermMod = AASequence::fromString(".(Dimethyl)DFPIANGER");
  AASequence aaCtermMod = AASequence::fromString("DFPIANGER.(Label:18O(2))");
  TEST_EQUAL(aaNoMod.isModified(), false)
  TEST_EQUAL(aaNtermMod.isModified(), true)
  TEST_EQUAL(aaCtermMod.isModified(), true)
  TEST_EQUAL(aaNoMod.getNTerminalModificationName(), "")
  TEST_EQUAL(aaNtermMod.getNTerminalModificationName(), "Dimethyl")
  TEST_EQUAL(aaCtermMod.getNTerminalModificationName(), "")
  TEST_EQUAL(aaNoMod.getCTerminalModificationName(), "")
  TEST_EQUAL(aaNtermMod.getCTerminalModificationName(), "")
  TEST_EQUAL(aaCtermMod.getCTerminalModificationName(), "Label:18O(2)")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
