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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

///////////////////////////
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModifiedPeptideGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModifiedPeptideGenerator* ptr = nullptr;
ModifiedPeptideGenerator* null_ptr = nullptr;
START_SECTION(ModifiedPeptideGenerator())
{
  ptr = new ModifiedPeptideGenerator();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ModifiedPeptideGenerator())
{
  delete ptr;
}
END_SECTION

START_SECTION((static void applyFixedModifications(const std::vector< ResidueModification >::const_iterator &fixed_mods_begin, const std::vector< ResidueModification >::const_iterator& fixed_mods_end, AASequence& peptide)))
{
  // query modification of interest from ModificationsDB
  StringList modNames;
  modNames << "Carbamidomethyl (C)";
  vector<ResidueModification> fixed_mods;
  for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
  {
    String modification(*mod_it);
    fixed_mods.push_back( ModificationsDB::getInstance()->getModification(modification));
  }

  AASequence seq0 = AASequence::fromString("AAAACAAAA"); // exactly one target site
  AASequence seq1 = AASequence::fromString("AAAAAAAAA"); // no target site
  AASequence seq2 = AASequence::fromString("CCCCCCCCC"); // all target sites
  AASequence seq3 = AASequence::fromString("AAAACAAC(Carbamidomethyl)AAA"); // one of two target sites already modified
  AASequence seq4 = AASequence::fromString("AAAACAAC(Oxidation)AAA"); // one of two target sites already modified

  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods.begin(), fixed_mods.end(), seq0);
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods.begin(), fixed_mods.end(), seq1);
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods.begin(), fixed_mods.end(), seq2);
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods.begin(), fixed_mods.end(), seq3);
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods.begin(), fixed_mods.end(), seq4);

  TEST_EQUAL(seq0.toString(), "AAAAC(Carbamidomethyl)AAAA");
  TEST_EQUAL(seq1.toString(), "AAAAAAAAA");
  TEST_EQUAL(seq2.toString(), "C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)");
  TEST_EQUAL(seq3.toString(), "AAAAC(Carbamidomethyl)AAC(Carbamidomethyl)AAA");
  TEST_EQUAL(seq4.toString(), "AAAAC(Carbamidomethyl)AAC(Oxidation)AAA");

   // test terminal modifications
   modNames.clear();
   modNames << "Carbamyl (N-term)";
   fixed_mods.clear();
   for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
   {
     String modification(*mod_it);
     fixed_mods.push_back(ModificationsDB::getInstance()->getModification(modification));
   }
  seq0 = AASequence::fromString("KAAAAAAAA"); // exactly one target site
  seq1 = AASequence::fromString("K(Carbamyl)AAAAAAAA"); // ambigous case: is mod Carbamyl (K) or (N-Term)?
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods.begin(), fixed_mods.end(), seq0);
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods.begin(), fixed_mods.end(), seq1);
  TEST_EQUAL(seq0.toString(), ".(Carbamyl)KAAAAAAAA");
  TEST_EQUAL(seq1.toString(), ".(Carbamyl)K(Carbamyl)AAAAAAAA");
 
}
END_SECTION

START_SECTION((static void applyVariableModifications(const std::vector< ResidueModification >::const_iterator &var_mods_begin, const std::vector< ResidueModification >::const_iterator &var_mods_end, const AASequence& peptide, Size max_variable_mods_per_peptide, std::vector< AASequence > &all_modified_peptides, bool keep_unmodified=true)))
{
  // query modification of interest from ModificationsDB
  StringList modNames;
  modNames << "Oxidation (M)";
  vector<ResidueModification> fixed_mods;
  for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
  {
    String modification(*mod_it);
    fixed_mods.push_back( ModificationsDB::getInstance()->getModification(modification));
  }

  vector<AASequence> modified_peptides;

  // test behavior if sequence empty
  AASequence seq;
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();

  // test behavior if peptide empty
  seq = AASequence::fromString("");
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 0, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 2, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();

  // test behavior if no target site in sequence
  seq = AASequence::fromString("AAAAAAAAA");
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();

  // test flag to preserve passed peptide
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 1, modified_peptides, true);
  TEST_EQUAL(modified_peptides[0], AASequence::fromString("AAAAAAAAA")); // only the original peptide
  modified_peptides.clear();

  // test behavior if one target site in sequence and different number of maximum variable modifications are choosen
  seq = AASequence::fromString("AAAAMAAAA"); // exactly one target site
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 0, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 1); // one generated peptide
  TEST_EQUAL(modified_peptides[0].toString(), "AAAAM(Oxidation)AAAA");
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 2, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 1); // also only one generated peptide as there is only one target site
  TEST_EQUAL(modified_peptides[0].toString(), "AAAAM(Oxidation)AAAA");
  modified_peptides.clear();
  // test again keeping of the original peptides
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 1, modified_peptides, true);
  TEST_EQUAL(modified_peptides.size(), 2); // original and modified peptide
  TEST_EQUAL(modified_peptides[0].toString(), "AAAAMAAAA");
  TEST_EQUAL(modified_peptides[1].toString(), "AAAAM(Oxidation)AAAA");
  modified_peptides.clear();

  // test behavior if two target sites are in peptide and we need some combinatorics
  // only one additional variable modification
  seq = AASequence::fromString("AAMAAAMAA"); // two target site
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 2); // two modified peptides each modified at a different site
  TEST_EQUAL(modified_peptides[0].toString(), "AAMAAAM(Oxidation)AA");
  TEST_EQUAL(modified_peptides[1].toString(), "AAM(Oxidation)AAAMAA");
  modified_peptides.clear();

  // up to two variable modifications per petide
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 2, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 3); // three modified peptides, 1 modifed at first M, 1 modifed at second M and both M modified
  TEST_EQUAL(modified_peptides[0].toString(), "AAMAAAM(Oxidation)AA");
  TEST_EQUAL(modified_peptides[1].toString(), "AAM(Oxidation)AAAMAA");
  TEST_EQUAL(modified_peptides[2].toString(), "AAM(Oxidation)AAAM(Oxidation)AA");
  modified_peptides.clear();

  // two different modifications with same target site
  modNames.clear();
  modNames << "Glutathione (C)" << "Carbamidomethyl (C)";
  fixed_mods.clear();
  for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
  {
    String modification(*mod_it);
    fixed_mods.push_back( ModificationsDB::getInstance()->getModification(modification));
  }

  seq = AASequence::fromString("ACAACAACA"); // three target sites
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 6); // six modified peptides, as both C modifications can occur at all 3 sites
  TEST_EQUAL(modified_peptides[0].toString(), "ACAACAAC(Glutathione)A");
  TEST_EQUAL(modified_peptides[1].toString(), "ACAACAAC(Carbamidomethyl)A");
  TEST_EQUAL(modified_peptides[2].toString(), "ACAAC(Glutathione)AACA");
  TEST_EQUAL(modified_peptides[3].toString(), "ACAAC(Carbamidomethyl)AACA");
  TEST_EQUAL(modified_peptides[4].toString(), "AC(Glutathione)AACAACA");
  TEST_EQUAL(modified_peptides[5].toString(), "AC(Carbamidomethyl)AACAACA");
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 1, modified_peptides, true);
  TEST_EQUAL(modified_peptides.size(), 7); // same as above but +1 unmodified peptide

  modified_peptides.clear();

  seq = AASequence::fromString("ACAACAACA"); // three target sites and maximum of three occurances of the two modifications Glutathione and Carbamidomethyl
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 3, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 3*3*3-1); // three sites with 3 possibilities (none, Glut., Carb.) each - but we need to subtract (none, none, none). The unmodified peptide

  modified_peptides.clear();
  seq = AASequence::fromString("AAAAC(Carbamidomethyl)AAAA"); // target site already occupied
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 3, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide as target site was already occupied

  // three different modifications
  modNames.clear();
  modNames << "Glutathione (C)" << "Carbamidomethyl (C)" << "Oxidation (M)";
  fixed_mods.clear();
  for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
  {
    String modification(*mod_it);
    fixed_mods.push_back( ModificationsDB::getInstance()->getModification(modification));
  }
  modified_peptides.clear();

  seq = AASequence::fromString("ACMACMACA"); // three target sites (C) and two target sites (M)
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 3, modified_peptides, false);

  // for exactly one mod: 2*3 + 2 = 8 (two possibilities each for each C sites and 1 for each of the two M sites)
  // for exactly two mods: 1 (iff two Ox) + 6 * 2 (iff one Ox at two pos.) + (3 choose 2) * 4 (iff no Ox) = 25
  // for exactly three mod: 6 (iff two Ox) + ((3 choose2) * 4) * 2 (iff one Ox at two pos.) + 8 (iff no Ox) = 38
  TEST_EQUAL(modified_peptides.size(), 8 + 25 + 38);

  // test terminal modifications
  modNames.clear();
  modNames << "Carbamyl (N-term)" << "Oxidation (M)";
  fixed_mods.clear();
  for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
  {
    String modification(*mod_it);
    fixed_mods.push_back(ModificationsDB::getInstance()->getModification(modification));
  }

  modified_peptides.clear();
  seq = AASequence::fromString("KAAAAAAAMA"); // exactly one target site
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 2, modified_peptides, true);
  TEST_EQUAL(modified_peptides.size(), 4);
  TEST_EQUAL(modified_peptides[0].toString(), "KAAAAAAAMA");
  TEST_EQUAL(modified_peptides[1].toString(), "KAAAAAAAM(Oxidation)A");
  TEST_EQUAL(modified_peptides[2].toString(), ".(Carbamyl)KAAAAAAAMA");
  TEST_EQUAL(modified_peptides[3].toString(), ".(Carbamyl)KAAAAAAAM(Oxidation)A");
  
  modified_peptides.clear();
  seq = AASequence::fromString("K(Carbamyl)AAAAAAAMA"); 
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seq, 2, modified_peptides, true);
  TEST_EQUAL(modified_peptides.size(), 4);
  TEST_EQUAL(modified_peptides[0].toString(), "K(Carbamyl)AAAAAAAMA");
  TEST_EQUAL(modified_peptides[1].toString(), "K(Carbamyl)AAAAAAAM(Oxidation)A");
  TEST_EQUAL(modified_peptides[2].toString(), ".(Carbamyl)K(Carbamyl)AAAAAAAMA");
  TEST_EQUAL(modified_peptides[3].toString(), ".(Carbamyl)K(Carbamyl)AAAAAAAM(Oxidation)A");
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



