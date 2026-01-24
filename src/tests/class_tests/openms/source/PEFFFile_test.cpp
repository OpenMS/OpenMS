// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/PEFFFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <fstream>
#include <vector>

///////////////////////////

START_TEST(PEFFFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PEFFFile* ptr = nullptr;

START_SECTION((PEFFFile()))
  ptr = new PEFFFile();
  TEST_NOT_EQUAL(ptr, nullptr)
END_SECTION

START_SECTION((~PEFFFile()))
  delete ptr;
END_SECTION

START_SECTION((void load(const String& filename, std::vector<PEFFEntry>& entries, std::vector<PEFFDatabaseMetadata>& headers) const))
{
  vector<PEFFEntry> entries;
  vector<PEFFDatabaseMetadata> headers;
  PEFFFile file;

  // Test file not found
  TEST_EXCEPTION(Exception::FileNotFound, file.load("PEFFFile_test_this_file_does_not_exist.peff", entries, headers))

  // Load test file
  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_test.peff"), entries, headers);

  // Check header
  TEST_EQUAL(headers.size(), 1)
  TEST_EQUAL(headers[0].db_name, "TestDB")
  TEST_EQUAL(headers[0].prefix, "test")
  TEST_EQUAL(headers[0].db_description, "Test database for PEFF file parsing")
  TEST_EQUAL(headers[0].is_decoy, false)
  TEST_EQUAL(headers[0].db_sources.size(), 1)
  TEST_EQUAL(headers[0].db_version, "1.0")
  TEST_EQUAL(headers[0].db_date, "20240115")
  TEST_EQUAL(headers[0].number_of_entries, 3)
  TEST_EQUAL(headers[0].sequence_type, PEFFDatabaseMetadata::SequenceType::AA)
  TEST_EQUAL(headers[0].general_comment, "This is a test file")

  // Check entries
  TEST_EQUAL(entries.size(), 3)

  // First entry
  TEST_EQUAL(entries[0].identifier, "test:P12345")
  TEST_EQUAL(entries[0].sequence, "MSEQVENCEOFTHEFIRSTPROTEINWITHSOMEMODIFICATIONSOK")
  TEST_EQUAL(entries[0].protein_names.size(), 1)
  TEST_EQUAL(entries[0].protein_names[0], "Test Protein 1")
  TEST_EQUAL(entries[0].gene_name, "TEST1")
  TEST_EQUAL(entries[0].ncbi_tax_id, 9606)
  TEST_EQUAL(entries[0].taxonomy_name, "Homo sapiens")
  TEST_EQUAL(entries[0].sequence_length, 49)
  TEST_EQUAL(entries[0].sequence_version, "1")
  TEST_EQUAL(entries[0].entry_version, "1")
  TEST_EQUAL(entries[0].protein_existence, 1)

  // Check modifications
  TEST_EQUAL(entries[0].modifications.size(), 2)
  TEST_EQUAL(entries[0].modifications[0].position, 10)
  TEST_EQUAL(entries[0].modifications[0].accession, "UNIMOD:35")
  TEST_EQUAL(entries[0].modifications[0].name, "Oxidation")
  TEST_EQUAL(entries[0].modifications[0].type, PEFFModification::Type::UNIMOD)
  TEST_EQUAL(entries[0].modifications[1].position, 25)
  TEST_EQUAL(entries[0].modifications[1].accession, "MOD:00046")
  TEST_EQUAL(entries[0].modifications[1].name, "Phosphorylation")
  TEST_EQUAL(entries[0].modifications[1].type, PEFFModification::Type::PSI_MOD)

  // Check simple variants
  TEST_EQUAL(entries[0].simple_variants.size(), 1)
  TEST_EQUAL(entries[0].simple_variants[0].position, 5)
  TEST_EQUAL(entries[0].simple_variants[0].variant_aa, 'R')
  TEST_EQUAL(entries[0].simple_variants[0].sources, "dbSNP:rs123456")

  // Check processed regions
  TEST_EQUAL(entries[0].processed_regions.size(), 1)
  TEST_EQUAL(entries[0].processed_regions[0].start_position, 1)
  TEST_EQUAL(entries[0].processed_regions[0].end_position, 20)
  TEST_EQUAL(entries[0].processed_regions[0].type, "PEFF:0001021")
  TEST_EQUAL(entries[0].processed_regions[0].name, "signal peptide")

  // Second entry with multiple protein names and complex variant
  TEST_EQUAL(entries[1].identifier, "test:P12346")
  TEST_EQUAL(entries[1].protein_names.size(), 2)
  TEST_EQUAL(entries[1].protein_names[0], "Test Protein 2")
  TEST_EQUAL(entries[1].protein_names[1], "Alternative Name")
  TEST_EQUAL(entries[1].complex_variants.size(), 1)
  TEST_EQUAL(entries[1].complex_variants[0].start_position, 5)
  TEST_EQUAL(entries[1].complex_variants[0].end_position, 10)
  TEST_EQUAL(entries[1].complex_variants[0].replacement, "LONGER")

  // Third entry (minimal annotations)
  TEST_EQUAL(entries[2].identifier, "test:P12347")
  TEST_EQUAL(entries[2].protein_names.size(), 1)
  TEST_EQUAL(entries[2].protein_names[0], "Simple Protein")
}
END_SECTION

START_SECTION((void store(const String& filename, const std::vector<PEFFEntry>& entries, const PEFFDatabaseMetadata& header) const))
{
  vector<PEFFEntry> entries, entries2;
  vector<PEFFDatabaseMetadata> headers, headers2;
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  PEFFFile file;

  // Load test file
  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_test.peff"), entries, headers);

  // Test unable to create file
  TEST_EXCEPTION(Exception::UnableToCreateFile, file.store("/bla/bluff/blblb/sdfhsdjf/test.peff", entries, headers[0]))

  // Store and reload
  file.store(tmp_filename, entries, headers[0]);
  file.load(tmp_filename, entries2, headers2);

  // Compare entries
  TEST_EQUAL(entries.size(), entries2.size())
  for (Size i = 0; i < entries.size(); ++i)
  {
    TEST_EQUAL(entries[i].identifier, entries2[i].identifier)
    TEST_EQUAL(entries[i].sequence, entries2[i].sequence)
    TEST_EQUAL(entries[i].protein_names.size(), entries2[i].protein_names.size())
    TEST_EQUAL(entries[i].modifications.size(), entries2[i].modifications.size())
    TEST_EQUAL(entries[i].simple_variants.size(), entries2[i].simple_variants.size())
    TEST_EQUAL(entries[i].complex_variants.size(), entries2[i].complex_variants.size())
    TEST_EQUAL(entries[i].processed_regions.size(), entries2[i].processed_regions.size())
  }
}
END_SECTION

START_SECTION((void readStart(const String& filename)))
{
  PEFFFile file;
  TEST_EXCEPTION(Exception::FileNotFound, file.readStart("nonexistent.peff"))

  file.readStart(OPENMS_GET_TEST_DATA_PATH("PEFFFile_test.peff"));
  const auto& headers = file.getHeaders();
  TEST_EQUAL(headers.size(), 1)
  TEST_EQUAL(headers[0].db_name, "TestDB")
}
END_SECTION

START_SECTION((bool readNext(PEFFEntry& entry)))
{
  PEFFFile file;
  PEFFEntry entry;

  file.readStart(OPENMS_GET_TEST_DATA_PATH("PEFFFile_test.peff"));

  TEST_EQUAL(file.readNext(entry), true)
  TEST_EQUAL(entry.identifier, "test:P12345")

  TEST_EQUAL(file.readNext(entry), true)
  TEST_EQUAL(entry.identifier, "test:P12346")

  TEST_EQUAL(file.readNext(entry), true)
  TEST_EQUAL(entry.identifier, "test:P12347")

  TEST_EQUAL(file.readNext(entry), false)
}
END_SECTION

START_SECTION((const std::vector<PEFFDatabaseMetadata>& getHeaders() const))
{
  PEFFFile file;
  file.readStart(OPENMS_GET_TEST_DATA_PATH("PEFFFile_test.peff"));
  const auto& headers = file.getHeaders();
  TEST_EQUAL(headers.size(), 1)
  TEST_EQUAL(headers[0].version, "1.0")
}
END_SECTION

START_SECTION((bool atEnd() const))
{
  PEFFFile file;
  PEFFEntry entry;

  file.readStart(OPENMS_GET_TEST_DATA_PATH("PEFFFile_minimal.peff"));
  TEST_EQUAL(file.atEnd(), false)

  file.readNext(entry);
  TEST_EQUAL(file.atEnd(), true)
}
END_SECTION

START_SECTION((void writeStart(const String& filename, const PEFFDatabaseMetadata& header)))
{
  PEFFFile file;
  PEFFDatabaseMetadata header;
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  file.writeStart(tmp_filename, header);
  file.writeEnd();

  // Verify file was created
  std::ifstream test_file(tmp_filename.c_str());
  TEST_EQUAL(test_file.good(), true)
}
END_SECTION

START_SECTION((void writeNext(const PEFFEntry& entry)))
{
  PEFFFile file;
  PEFFDatabaseMetadata header;
  header.db_name = "WriteTest";
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  PEFFEntry entry;
  entry.identifier = "test:WRITE001";
  entry.sequence = "TESTSEQUENCE";
  entry.protein_names.push_back("Write Test Protein");

  file.writeStart(tmp_filename, header);
  file.writeNext(entry);
  file.writeEnd();

  // Reload and verify
  vector<PEFFEntry> entries;
  vector<PEFFDatabaseMetadata> headers;
  file.load(tmp_filename, entries, headers);

  TEST_EQUAL(entries.size(), 1)
  TEST_EQUAL(entries[0].identifier, "test:WRITE001")
  TEST_EQUAL(entries[0].sequence, "TESTSEQUENCE")
}
END_SECTION

START_SECTION((void writeEnd()))
{
  // Tested in writeStart and writeNext sections
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static bool isPEFFFile(const String& filename)))
{
  TEST_EQUAL(PEFFFile::isPEFFFile(OPENMS_GET_TEST_DATA_PATH("PEFFFile_test.peff")), true)
  TEST_EQUAL(PEFFFile::isPEFFFile(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta")), false)
  TEST_EQUAL(PEFFFile::isPEFFFile("nonexistent.peff"), false)
}
END_SECTION

START_SECTION((static String toProForma(const PEFFEntry& entry)))
{
  PEFFEntry entry;
  entry.sequence = "PEPTIDE";
  String proforma = PEFFFile::toProForma(entry);
  TEST_EQUAL(proforma, "PEPTIDE")
}
END_SECTION

START_SECTION([PEFFEntry] FASTAFile::FASTAEntry toFASTAEntry() const)
{
  PEFFEntry entry;
  entry.identifier = "test:P12345";
  entry.sequence = "SEQUENCE";
  entry.protein_names.push_back("Test Protein");

  FASTAFile::FASTAEntry fasta = entry.toFASTAEntry();
  TEST_EQUAL(fasta.identifier, "test:P12345")
  TEST_EQUAL(fasta.sequence, "SEQUENCE")
  TEST_EQUAL(fasta.description, "Test Protein")
}
END_SECTION

START_SECTION([PEFFEntry] static PEFFEntry fromFASTAEntry(const FASTAFile::FASTAEntry&))
{
  FASTAFile::FASTAEntry fasta("test:P12345", "Test Protein", "SEQUENCE");
  PEFFEntry entry = PEFFEntry::fromFASTAEntry(fasta);

  TEST_EQUAL(entry.identifier, "test:P12345")
  TEST_EQUAL(entry.sequence, "SEQUENCE")
  TEST_EQUAL(entry.protein_names.size(), 1)
  TEST_EQUAL(entry.protein_names[0], "Test Protein")
  TEST_EQUAL(entry.sequence_length, 8)
}
END_SECTION

START_SECTION([PEFFEntry] AASequence getSequence() const)
{
  PEFFEntry entry;
  entry.sequence = "PEPTIDE";
  AASequence seq = entry.getSequence();
  TEST_EQUAL(seq.toString(), "PEPTIDE")
  TEST_EQUAL(seq.size(), 7)
}
END_SECTION

START_SECTION([PEFFEntry] AASequence getModifiedSequence() const)
{
  PEFFEntry entry;
  entry.sequence = "PEPTIDE";
  // Note: This test just verifies the method works; actual modification
  // application would require valid modification names in ModificationsDB
  AASequence seq = entry.getModifiedSequence();
  TEST_EQUAL(seq.size(), 7)
}
END_SECTION

START_SECTION([PEFFEntry] getVariantSequences)
{
  PEFFEntry entry;
  entry.sequence = "PEPTIDE";
  entry.simple_variants.push_back(PEFFVariantSimple(3, 'A', "dbSNP:test"));

  auto variants = entry.getVariantSequences();
  TEST_EQUAL(variants.size(), 1)
  if (!variants.empty())
  {
    TEST_EQUAL(variants[0].second.toString(), "PEATIDE")
    TEST_EQUAL(variants[0].first.hasSubstring("P3A"), true)
  }
}
END_SECTION

START_SECTION([PEFFEntry] AASequence getProcessedSequence(const String& region_type) const)
{
  PEFFEntry entry;
  entry.sequence = "SIGNALPEPTIDERESTOFPROTEIN";
  entry.processed_regions.push_back(PEFFProcessedRegion(1, 13, "PEFF:0001021", "signal peptide"));

  // Get mature protein (after signal peptide)
  AASequence mature = entry.getProcessedSequence("PEFF:0001021");
  TEST_EQUAL(mature.toString(), "RESTOFPROTEIN")
}
END_SECTION

START_SECTION([PEFFModification] PEFFModification())
{
  PEFFModification mod;
  TEST_EQUAL(mod.position, 0)
  TEST_EQUAL(mod.accession, "")
  TEST_EQUAL(mod.name, "")
  TEST_EQUAL(mod.evidence, "")
  TEST_EQUAL(mod.type, PEFFModification::Type::GENERIC)
}
END_SECTION

START_SECTION([PEFFModification] PEFFModification(Size pos, const String& acc, const String& n, const String& ev))
{
  PEFFModification mod(10, "UNIMOD:35", "Oxidation", "experimental");
  TEST_EQUAL(mod.position, 10)
  TEST_EQUAL(mod.accession, "UNIMOD:35")
  TEST_EQUAL(mod.name, "Oxidation")
  TEST_EQUAL(mod.evidence, "experimental")
  TEST_EQUAL(mod.type, PEFFModification::Type::UNIMOD)

  PEFFModification mod2(5, "MOD:00046", "Phospho");
  TEST_EQUAL(mod2.type, PEFFModification::Type::PSI_MOD)
}
END_SECTION

START_SECTION([PEFFVariantSimple] PEFFVariantSimple())
{
  PEFFVariantSimple var;
  TEST_EQUAL(var.position, 0)
  TEST_EQUAL(var.variant_aa, '\0')
  TEST_EQUAL(var.sources, "")
}
END_SECTION

START_SECTION([PEFFVariantComplex] PEFFVariantComplex())
{
  PEFFVariantComplex var;
  TEST_EQUAL(var.start_position, 0)
  TEST_EQUAL(var.end_position, 0)
  TEST_EQUAL(var.replacement, "")
  TEST_EQUAL(var.sources, "")
}
END_SECTION

START_SECTION([PEFFProcessedRegion] PEFFProcessedRegion())
{
  PEFFProcessedRegion reg;
  TEST_EQUAL(reg.start_position, 0)
  TEST_EQUAL(reg.end_position, 0)
  TEST_EQUAL(reg.type, "")
  TEST_EQUAL(reg.name, "")
  TEST_EQUAL(reg.description, "")
}
END_SECTION

START_SECTION([PEFFDatabaseMetadata] PEFFDatabaseMetadata())
{
  PEFFDatabaseMetadata meta;
  TEST_EQUAL(meta.version, "1.0")
  TEST_EQUAL(meta.db_name, "")
  TEST_EQUAL(meta.is_decoy, false)
  TEST_EQUAL(meta.sequence_type, PEFFDatabaseMetadata::SequenceType::AA)
}
END_SECTION

START_SECTION(test_minimal_file)
{
  vector<PEFFEntry> entries;
  vector<PEFFDatabaseMetadata> headers;
  PEFFFile file;

  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_minimal.peff"), entries, headers);

  TEST_EQUAL(entries.size(), 1)
  TEST_EQUAL(entries[0].identifier, "sp:P00001")
  TEST_EQUAL(entries[0].protein_names.size(), 1)
  TEST_EQUAL(entries[0].protein_names[0], "Minimal Test")
  TEST_EQUAL(entries[0].sequence, "SEQUENCE")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
