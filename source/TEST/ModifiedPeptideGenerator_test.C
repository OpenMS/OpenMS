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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

///////////////////////////
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModifiedPeptideGenerator, "$Id ModifiedPeptideGenerator_test.C 11707 2013-08-28 15:35:07Z timosachsenberg $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ModifiedPeptideGenerator* ptr = 0;
ModifiedPeptideGenerator* null_ptr = 0;
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

START_SECTION((static void applyFixedModifications(const std::vector< ResidueModification >::const_iterator &fixed_mods_begin, const std::vector< ResidueModification >::const_iterator &fixed_mods_end, std::vector< AASequence >::iterator digested_peptides_begin, std::vector< AASequence >::iterator digested_peptides_end)))
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

  vector<AASequence> seqs;
  seqs.push_back("AAAACAAAA"); // exactly one target site 
  seqs.push_back("AAAAAAAAA"); // no target site
  seqs.push_back("CCCCCCCCC"); // all target sites
  seqs.push_back("AAAACAAC(Carbamidomethyl)AAA"); // one of two target sites already modified
  seqs.push_back("AAAACAAC(Oxidation)AAA"); // one of two target sites already modified
  
  ModifiedPeptideGenerator::applyFixedModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end());
  TEST_EQUAL(seqs[0].toString(), "AAAAC(Carbamidomethyl)AAAA");  
  TEST_EQUAL(seqs[1].toString(), "AAAAAAAAA");  
  TEST_EQUAL(seqs[2].toString(), "C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)C(Carbamidomethyl)");  
  TEST_EQUAL(seqs[3].toString(), "AAAAC(Carbamidomethyl)AAC(Carbamidomethyl)AAA");  
  TEST_EQUAL(seqs[4].toString(), "AAAAC(Carbamidomethyl)AAC(Oxidation)AAA");  
}
END_SECTION
 
START_SECTION((static void applyVariableModifications(const std::vector< ResidueModification >::const_iterator &var_mods_begin, const std::vector< ResidueModification >::const_iterator &var_mods_end, std::vector< AASequence >::iterator digested_peptides_begin, std::vector< AASequence >::iterator digested_peptides_end, Size max_variable_mods_per_peptide, std::vector< AASequence > &all_modified_peptides, bool keep_unmodified=true)))
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

  vector<AASequence> seqs;
  vector<AASequence> modified_peptides;

  // test behavior if sequence empty
  seqs.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  // test behavior if peptide empty
  seqs.push_back("");
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 0, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 2, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  // test behavior if no target site in sequence
  seqs.clear();
  seqs.push_back("AAAAAAAAA");
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 1, modified_peptides, true);
  TEST_EQUAL(modified_peptides[0], "AAAAAAAAA"); // only the original peptide
  seqs.clear();
  modified_peptides.clear();

  // test behavior if one target site in sequence and different number of maximum variable modifications are choosen
  seqs.push_back("AAAAMAAAA"); // exactly one target site 
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 0, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 0); // no generated peptide
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 1); // one generated peptide
  TEST_EQUAL(modified_peptides[0].toString(), "AAAAM(Oxidation)AAAA");
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 2, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 1); // also only one generated peptide as there is only one target site
  TEST_EQUAL(modified_peptides[0].toString(), "AAAAM(Oxidation)AAAA");
  modified_peptides.clear();
  // test again keeping of the original peptides
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 1, modified_peptides, true);
  TEST_EQUAL(modified_peptides.size(), 2); // original and modified peptide
  TEST_EQUAL(modified_peptides[0].toString(), "AAAAMAAAA");
  TEST_EQUAL(modified_peptides[1].toString(), "AAAAM(Oxidation)AAAA");
  modified_peptides.clear();
  seqs.clear();

  // test behavior if two target sites are in peptide and we need some combinatorics
  // only one additional variable modification 
  seqs.push_back("AAMAAAMAA"); // two target site
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 2); // two modified peptides each modified at a different site
  TEST_EQUAL(modified_peptides[0].toString(), "AAMAAAM(Oxidation)AA");
  TEST_EQUAL(modified_peptides[1].toString(), "AAM(Oxidation)AAAMAA");
  modified_peptides.clear();

  // up to two variable modifications per petide
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 2, modified_peptides, false);
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

  seqs.clear();
  seqs.push_back("ACAACAACA"); // three target sites
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 1, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 6); // six modified peptides, as both C modifications can occur at all 3 sites 
  TEST_EQUAL(modified_peptides[0].toString(), "ACAACAAC(Glutathione)A");
  TEST_EQUAL(modified_peptides[1].toString(), "ACAACAAC(Carbamidomethyl)A");
  TEST_EQUAL(modified_peptides[2].toString(), "ACAAC(Glutathione)AACA");
  TEST_EQUAL(modified_peptides[3].toString(), "ACAAC(Carbamidomethyl)AACA");
  TEST_EQUAL(modified_peptides[4].toString(), "AC(Glutathione)AACAACA");
  TEST_EQUAL(modified_peptides[5].toString(), "AC(Carbamidomethyl)AACAACA");
  modified_peptides.clear();
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 1, modified_peptides, true);
  TEST_EQUAL(modified_peptides.size(), 7); // same as above but +1 unmodified peptide

  seqs.clear();
  modified_peptides.clear();

  seqs.push_back("ACAACAACA"); // three target sites and maximum of three occurances of the two modifications Glutathione and Carbamidomethyl
  ModifiedPeptideGenerator::applyVariableModifications(fixed_mods.begin(), fixed_mods.end(), seqs.begin(), seqs.end(), 3, modified_peptides, false);
  TEST_EQUAL(modified_peptides.size(), 3*3*3-1); // three sites with 3 possibilities (none, Glut., Carb.) each - but we need to subtract (none, none, none). The unmodified peptide  
}
END_SECTION
  
  

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



