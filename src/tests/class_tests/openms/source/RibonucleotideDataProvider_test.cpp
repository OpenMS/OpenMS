// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CHEMISTRY/RibonucleotideDataProvider.h>
#include <OpenMS/CHEMISTRY/ModomicsJSONDataProvider.h>
#include <OpenMS/CHEMISTRY/RibonucleotideTSVDataProvider.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <fstream>

using namespace OpenMS;
using namespace std;

START_TEST(RibonucleotideDataProvider, "$Id$")

/////////////////////////////////////////////////////////////
// InMemoryRibonucleotideDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(InMemoryRibonucleotideDataProvider::loadRibonucleotides())
{
  vector<RibonucleotideEntry> entries;

  RibonucleotideEntry e1;
  e1.ribo = make_unique<Ribonucleotide>();
  e1.ribo->setCode("TestA");
  e1.ribo->setName("Test Adenosine");
  entries.push_back(std::move(e1));

  RibonucleotideEntry e2;
  e2.ribo = make_unique<Ribonucleotide>();
  e2.ribo->setCode("TestC");
  e2.ribo->setName("Test Cytidine");
  entries.push_back(std::move(e2));

  InMemoryRibonucleotideDataProvider provider(std::move(entries));
  auto result = provider.loadRibonucleotides();

  TEST_EQUAL(result.size(), 2)
  TEST_STRING_EQUAL(result[0].ribo->getCode(), "TestA")
  TEST_STRING_EQUAL(result[0].ribo->getName(), "Test Adenosine")
  TEST_STRING_EQUAL(result[1].ribo->getCode(), "TestC")
  TEST_STRING_EQUAL(result[1].ribo->getName(), "Test Cytidine")
  TEST_EQUAL(result[0].isAmbiguous(), false)
  TEST_EQUAL(result[1].isAmbiguous(), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ModomicsJSONDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(ModomicsJSONDataProvider::loadRibonucleotides())
{
  ModomicsJSONDataProvider provider("CHEMISTRY/Modomics.json");
  auto entries = provider.loadRibonucleotides();

  // Modomics JSON should contain many ribonucleotides
  TEST_EQUAL(entries.size() > 100, true)

  // Check that a well-known entry exists: Am = 2'-O-methyladenosine
  bool found_am = false;
  for (const auto& e : entries)
  {
    if (e.ribo->getCode() == "Am")
    {
      found_am = true;
      TEST_STRING_EQUAL(e.ribo->getName(), "2'-O-methyladenosine")
      break;
    }
  }
  TEST_EQUAL(found_am, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// RibonucleotideTSVDataProvider tests
/////////////////////////////////////////////////////////////

START_SECTION(RibonucleotideTSVDataProvider::loadRibonucleotides())
{
  RibonucleotideTSVDataProvider provider("CHEMISTRY/Custom_RNA_modifications.tsv");
  auto entries = provider.loadRibonucleotides();

  // TSV file should contain entries
  TEST_EQUAL(entries.size() > 0, true)

  // Check that the ambiguous entry "msU?" exists
  bool found_msu = false;
  for (const auto& e : entries)
  {
    if (e.ribo->getCode() == "msU?")
    {
      found_msu = true;
      TEST_EQUAL(e.isAmbiguous(), true)
      break;
    }
  }
  TEST_EQUAL(found_msu, true)

  // Explicit sugar formulas must reach both NASequence and the spectrum generator.
  // The extra notes column ensures that the formula is selected by header name.
  std::string filename;
  NEW_TMP_FILE(filename)
  {
    std::ofstream file(filename);
    file << "name\tshort_name\tnew_nomenclature\toriginating_base\trnamods_abbrev\thtml_abbrev\tformula\tmonoisotopic_mass\taverage_mass\talternatives\tnotes\tbaseloss_formula\n"
         << "2'-O-methyladenosine\tOMeA\t\tA\t\t\tC11H15N5O4\t0\t0\t\t\tC6H12O5\n"
         << "LNA adenosine\tLNA-A\t\tA\t\t\tC11H13N5O4\t0\t0\t\t\tC6H10O5\n"
         << "2'-fluorodeoxyadenosine\tdFA\t\tA\t\t\tC10H12FN5O3\t0\t0\t\t\tC5H9FO4\n"
         << "inferred methyladenosine\tAm\t\tA\t\t\tC11H15N5O4\t0\t0\t\t\t\n"
         << "inferred deoxyadenosine\tdA\t\tA\t\t\tC10H13N5O3\t0\t0\t\t\t-\n"
         << "invalid explicit formula\tbadA\t\tA\t\t\tC10H13N5O4\t0\t0\t\t\tinvalid_formula\n"
         << "ambiguous methyladenosine\tmA?\t\tA\t\t\tC11H15N5O4\t0\t0\tmA Am\t\tC5H10O5\n";
  }
  RibonucleotideTSVDataProvider custom_provider(filename);
  auto custom_entries = custom_provider.loadRibonucleotides();
  TEST_EQUAL(custom_entries.size(), 6)
  TEST_TRUE(custom_entries[5].isAmbiguous())
  TEST_STRING_EQUAL(custom_entries[5].alternative_code_1, "mA")
  TEST_STRING_EQUAL(custom_entries[5].alternative_code_2, "Am")
  TEST_EQUAL(custom_entries[3].ribo->getBaselossFormula(), EmpiricalFormula("C6H12O5"))
  TEST_EQUAL(custom_entries[4].ribo->getBaselossFormula(), EmpiricalFormula("C5H10O4"))
  const vector<EmpiricalFormula> sugar_formulas = {
    EmpiricalFormula("C6H12O5"), EmpiricalFormula("C6H10O5"), EmpiricalFormula("C5H9FO4")};
  const vector<EmpiricalFormula> ab_formulas = {
    EmpiricalFormula("C16H20N5O8P"), EmpiricalFormula("C16H18N5O8P"), EmpiricalFormula("C15H17FN5O7P")};
  const auto* db = RibonucleotideDB::getInstance();
  NucleicAcidSpectrumGenerator generator;
  auto parameters = generator.getParameters();
  parameters.setValue("add_b_ions", "false");
  parameters.setValue("add_y_ions", "false");
  parameters.setValue("add_a-B_ions", "true");
  parameters.setValue("add_metainfo", "true");
  generator.setParameters(parameters);
  for (Size i = 0; i < 3; ++i)
  {
    TEST_EQUAL(custom_entries[i].ribo->getBaselossFormula(), sugar_formulas[i])
    for (const std::string& prefix : {std::string("dA"), std::string("dA*")})
    {
      const EmpiricalFormula expected = ab_formulas[i] +
        (prefix == "dA*" ? EmpiricalFormula("SO-1") : EmpiricalFormula());
      // A '*' on the last nucleotide describes the outgoing bond; it must not
      // change a2-B, whose only retained bond belongs to the first nucleotide.
      Ribonucleotide terminal_thiol(*custom_entries[i].ribo);
      terminal_thiol.setCode(terminal_thiol.getCode() + "*");
      for (const Ribonucleotide* terminal : {custom_entries[i].ribo.get(), &terminal_thiol})
      {
        NASequence sequence({db->getRibonucleotide(prefix), terminal, db->getRibonucleotide("dA")}, nullptr, nullptr);
        TEST_EQUAL(sequence.getPrefix(2).getFormula(NASequence::AminusB), expected)
        for (Int charge : {-2, -1, 1, 2})
        {
          const EmpiricalFormula charged_formula = expected + EmpiricalFormula::hydrogen(charge);
          TEST_EQUAL(sequence.getPrefix(2).getFormula(NASequence::AminusB, charge), charged_formula)
          MSSpectrum spectrum;
          generator.getSpectrum(spectrum, sequence, charge, charge);
          TEST_EQUAL(spectrum.size(), 1)
          TEST_STRING_EQUAL(spectrum.getStringDataArrays()[0][0], "a2-B")
          const double expected_mz = (charged_formula.getMonoWeight() - charge * Constants::ELECTRON_MASS_U) / abs(charge);
          TEST_REAL_SIMILAR(spectrum[0].getMZ(), expected_mz)
          map<Int, MSSpectrum> multiple;
          generator.getMultipleSpectra(multiple, sequence, {charge}, charge);
          TEST_EQUAL(multiple.at(charge).size(), 1)
          TEST_REAL_SIMILAR(multiple.at(charge)[0].getMZ(), expected_mz)
        }
      }
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// RibonucleotideTSVDataProvider line-ending regressions
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA] RibonucleotideTSVDataProvider accepts platform-independent line endings)
{
  for (const std::string& newline : {std::string("\n"), std::string("\r\n"), std::string("\r")})
  {
    std::string filename;
    NEW_TMP_FILE(filename)
    {
      // Binary output preserves the requested bytes on Windows too.
      std::ofstream file(filename, std::ios::binary);
      file << "# Custom RNA modifications" << newline
           << "name\tshort_name\tnew_nomenclature\toriginating_base\trnamods_abbrev\thtml_abbrev\tformula\tmonoisotopic_mass\taverage_mass\talternatives\tbaseloss_formula" << newline
           << "LNA adenosine\tLNA-A\t\tA\t\t\tC11H13N5O4\t0\t0\t\tC6H10O5" << newline
           << "inferred methyladenosine\tAm\t\tA\t\t\tC11H15N5O4\t0\t0\t\t" << newline
           << "inferred deoxyadenosine\tdA\t\tA\t\t\tC10H13N5O3\t0\t0\t\t-" << newline
           << "omitted optional fields\tUm\t\tU\t\t\tC10H14N2O6\t0\t0";
    }

    RibonucleotideTSVDataProvider provider(filename);
    auto entries = provider.loadRibonucleotides();
    TEST_EQUAL(entries.size(), 4)
    if (entries.size() == 4)
    {
      TEST_EQUAL(entries[0].ribo->getBaselossFormula(), EmpiricalFormula("C6H10O5"))
      TEST_EQUAL(entries[1].ribo->getBaselossFormula(), EmpiricalFormula("C6H12O5"))
      TEST_EQUAL(entries[2].ribo->getBaselossFormula(), EmpiricalFormula("C5H10O4"))
      TEST_EQUAL(entries[3].ribo->getBaselossFormula(), EmpiricalFormula("C6H12O5"))
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// RibonucleotideDB with InMemoryRibonucleotideDataProvider
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA] RibonucleotideDB with InMemoryRibonucleotideDataProvider)
{
  vector<RibonucleotideEntry> entries;

  RibonucleotideEntry e1;
  e1.ribo = make_unique<Ribonucleotide>();
  e1.ribo->setCode("xA");
  e1.ribo->setName("Test Modified Adenosine");
  entries.push_back(std::move(e1));

  RibonucleotideEntry e2;
  e2.ribo = make_unique<Ribonucleotide>();
  e2.ribo->setCode("xC");
  e2.ribo->setName("Test Modified Cytidine");
  entries.push_back(std::move(e2));

  vector<unique_ptr<RibonucleotideDataProvider>> providers;
  providers.push_back(make_unique<InMemoryRibonucleotideDataProvider>(std::move(entries)));

  RibonucleotideDB db(std::move(providers));

  // verify getRibonucleotide works
  auto ribo_a = db.getRibonucleotide("xA");
  TEST_STRING_EQUAL(ribo_a->getCode(), "xA")
  TEST_STRING_EQUAL(ribo_a->getName(), "Test Modified Adenosine")

  auto ribo_c = db.getRibonucleotide("xC");
  TEST_STRING_EQUAL(ribo_c->getCode(), "xC")
  TEST_STRING_EQUAL(ribo_c->getName(), "Test Modified Cytidine")

  // verify iterator count
  Size count = 0;
  for (auto it = db.begin(); it != db.end(); ++it)
  {
    ++count;
  }
  TEST_EQUAL(count, 2)

  // verify ElementNotFound exception for unknown code
  TEST_EXCEPTION(Exception::ElementNotFound, db.getRibonucleotide("Nonexistent"))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
