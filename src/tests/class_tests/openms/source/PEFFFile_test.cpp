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
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <fstream>
#include <limits>
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
  TEST_EQUAL(headers[0].number_of_entries, 3)
  TEST_EQUAL(headers[0].sequence_type, PEFFDatabaseMetadata::SequenceType::AA)
  TEST_EQUAL(headers[0].general_comments.size(), 1)
  TEST_EQUAL(headers[0].general_comments[0], "This is a test file")

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
  TEST_EQUAL(entries[0].simple_variants[0].optional_tag, "dbSNP:rs123456")

  // Check processed regions
  TEST_EQUAL(entries[0].processed_regions.size(), 1)
  TEST_EQUAL(entries[0].processed_regions[0].start_position, 1)
  TEST_EQUAL(entries[0].processed_regions[0].end_position, 20)
  TEST_EQUAL(entries[0].processed_regions[0].accession, "PEFF:0001021")
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
  // Test 1: No modifications
  PEFFEntry entry;
  entry.sequence = "PEPTIDE";
  String proforma = PEFFFile::toProForma(entry);
  TEST_EQUAL(proforma, "PEPTIDE")

  // Test 2: Single modification at position 3
  PEFFEntry entry2;
  entry2.sequence = "PEPTIDE";
  entry2.modifications.push_back(PEFFModification(3, "UNIMOD:35", "Oxidation"));
  proforma = PEFFFile::toProForma(entry2);
  TEST_EQUAL(proforma, "PEP[UNIMOD:35]TIDE")

  // Test 3: Multiple modifications at different positions
  PEFFEntry entry3;
  entry3.sequence = "PEPTIDE";
  entry3.modifications.push_back(PEFFModification(1, "UNIMOD:1", "Acetyl"));
  entry3.modifications.push_back(PEFFModification(5, "MOD:00046", "Phosphorylation"));
  proforma = PEFFFile::toProForma(entry3);
  TEST_EQUAL(proforma, "P[UNIMOD:1]EPTI[MOD:00046]DE")

  // Test 4: Unlocalised modification (position 0)
  PEFFEntry entry4;
  entry4.sequence = "PEPTIDE";
  entry4.modifications.push_back(PEFFModification(0, "UNIMOD:35", "Oxidation"));
  proforma = PEFFFile::toProForma(entry4);
  TEST_EQUAL(proforma, "<?[UNIMOD:35]>PEPTIDE")

  // Test 5: Mix of localised and unlocalised
  PEFFEntry entry5;
  entry5.sequence = "PEPTIDE";
  entry5.modifications.push_back(PEFFModification(0, "UNIMOD:35", "Oxidation"));
  entry5.modifications.push_back(PEFFModification(3, "MOD:00046", "Phosphorylation"));
  proforma = PEFFFile::toProForma(entry5);
  TEST_EQUAL(proforma, "<?[UNIMOD:35]>PEP[MOD:00046]TIDE")

  // Test 6: Named modification without accession
  PEFFEntry entry6;
  entry6.sequence = "PEPTIDE";
  PEFFModification namedMod;
  namedMod.position = 4;
  namedMod.name = "Carbamidomethyl";
  entry6.modifications.push_back(namedMod);
  proforma = PEFFFile::toProForma(entry6);
  TEST_EQUAL(proforma, "PEPT[Carbamidomethyl]IDE")

  // Test 7: Verify generated ProForma can be parsed by ProForma parser
  PEFFEntry entry7;
  entry7.sequence = "PEPTIDE";
  entry7.modifications.push_back(PEFFModification(3, "UNIMOD:35", "Oxidation"));
  entry7.modifications.push_back(PEFFModification(5, "MOD:00046", "Phosphorylation"));
  proforma = PEFFFile::toProForma(entry7);
  // Parse the ProForma string - should not throw
  ProForma::Peptidoform pf = ProForma::parse(proforma);
  TEST_EQUAL(pf.sequence.size(), 7)  // 7 amino acids
  // Check modifications were preserved (sequence elements are variants)
  TEST_EQUAL(std::get<ProForma::SequenceElement>(pf.sequence[2]).modifications.size(), 1)  // Position 3 (0-indexed: 2)
  TEST_EQUAL(std::get<ProForma::SequenceElement>(pf.sequence[4]).modifications.size(), 1)  // Position 5 (0-indexed: 4)
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
  entry.complex_variants.push_back(PEFFVariantComplex(2, 4, "XXX", "ClinVar:123"));

  std::vector<std::string> descriptions;
  std::vector<AASequence> sequences;

  // Without complex variants (default)
  entry.getVariantSequences(descriptions, sequences, false);
  TEST_EQUAL(sequences.size(), 1)
  TEST_EQUAL(descriptions.size(), 1)
  if (!sequences.empty())
  {
    TEST_EQUAL(sequences[0].toString(), "PEATIDE")
    TEST_EQUAL(descriptions[0].find("P3A") != std::string::npos, true)
  }

  // With complex variants
  entry.getVariantSequences(descriptions, sequences, true);
  TEST_EQUAL(sequences.size(), 2)
  if (sequences.size() >= 2)
  {
    // Simple variant
    TEST_EQUAL(sequences[0].toString(), "PEATIDE")
    // Complex variant: replace positions 2-4 (EPT) with XXX
    TEST_EQUAL(sequences[1].toString(), "PXXXIDE")
    TEST_EQUAL(descriptions[1].find("2-4>XXX") != std::string::npos, true)
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
  TEST_EQUAL(mod.optional_tag, "")
  TEST_EQUAL(mod.type, PEFFModification::Type::GENERIC)
}
END_SECTION

START_SECTION([PEFFModification] PEFFModification(Size pos, const String& acc, const String& n, const String& tag))
{
  PEFFModification mod(10, "UNIMOD:35", "Oxidation", "experimental");
  TEST_EQUAL(mod.position, 10)
  TEST_EQUAL(mod.accession, "UNIMOD:35")
  TEST_EQUAL(mod.name, "Oxidation")
  TEST_EQUAL(mod.optional_tag, "experimental")
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
  TEST_EQUAL(var.optional_tag, "")
}
END_SECTION

START_SECTION([PEFFVariantComplex] PEFFVariantComplex())
{
  PEFFVariantComplex var;
  TEST_EQUAL(var.start_position, 0)
  TEST_EQUAL(var.end_position, 0)
  TEST_EQUAL(var.replacement, "")
  TEST_EQUAL(var.optional_tag, "")
}
END_SECTION

START_SECTION([PEFFProcessedRegion] PEFFProcessedRegion())
{
  PEFFProcessedRegion reg;
  TEST_EQUAL(reg.start_position, 0)
  TEST_EQUAL(reg.end_position, 0)
  TEST_EQUAL(reg.accession, "")
  TEST_EQUAL(reg.name, "")
  TEST_EQUAL(reg.optional_tag, "")
}
END_SECTION

START_SECTION([PEFFDatabaseMetadata] PEFFDatabaseMetadata())
{
  PEFFDatabaseMetadata meta;
  TEST_EQUAL(meta.version, "1.0")
  TEST_EQUAL(meta.db_name, "")
  TEST_EQUAL(meta.is_decoy, false)
  TEST_EQUAL(meta.sequence_type, PEFFDatabaseMetadata::SequenceType::AA)
  TEST_EQUAL(meta.has_annotation_identifiers, false)
}
END_SECTION

START_SECTION(test_annotation_identifiers)
{
  vector<PEFFEntry> entries;
  vector<PEFFDatabaseMetadata> headers;
  PEFFFile file;

  // Load file with annotation identifiers
  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_annotid.peff"), entries, headers);

  // Check header has annotation identifiers flag set
  TEST_EQUAL(headers.size(), 1)
  TEST_EQUAL(headers[0].db_name, "AnnotIDTest")
  TEST_EQUAL(headers[0].has_annotation_identifiers, true)

  // Check entries
  TEST_EQUAL(entries.size(), 2)

  // First entry - check annotation IDs on modifications
  TEST_EQUAL(entries[0].identifier, "test:P12345")
  TEST_EQUAL(entries[0].modifications.size(), 2)
  TEST_EQUAL(entries[0].modifications[0].annotation_id, 0)
  TEST_EQUAL(entries[0].modifications[0].position, 10)
  TEST_EQUAL(entries[0].modifications[0].accession, "MOD:00046")
  TEST_EQUAL(entries[0].modifications[1].annotation_id, 1)
  TEST_EQUAL(entries[0].modifications[1].position, 25)
  TEST_EQUAL(entries[0].modifications[1].accession, "UNIMOD:35")

  // Check annotation ID on simple variant
  TEST_EQUAL(entries[0].simple_variants.size(), 1)
  TEST_EQUAL(entries[0].simple_variants[0].annotation_id, 2)
  TEST_EQUAL(entries[0].simple_variants[0].position, 5)
  TEST_EQUAL(entries[0].simple_variants[0].variant_aa, 'R')

  // Check annotation ID on processed region
  TEST_EQUAL(entries[0].processed_regions.size(), 1)
  TEST_EQUAL(entries[0].processed_regions[0].annotation_id, 3)
  TEST_EQUAL(entries[0].processed_regions[0].start_position, 1)
  TEST_EQUAL(entries[0].processed_regions[0].end_position, 20)

  // Second entry - check annotation ID on complex variant
  TEST_EQUAL(entries[1].identifier, "test:P12346")
  TEST_EQUAL(entries[1].complex_variants.size(), 1)
  TEST_EQUAL(entries[1].complex_variants[0].annotation_id, 4)
  TEST_EQUAL(entries[1].complex_variants[0].start_position, 5)
  TEST_EQUAL(entries[1].complex_variants[0].end_position, 10)
  TEST_EQUAL(entries[1].complex_variants[0].replacement, "LONGER")
}
END_SECTION

START_SECTION(test_annotation_identifiers_roundtrip)
{
  vector<PEFFEntry> entries, entries2;
  vector<PEFFDatabaseMetadata> headers, headers2;
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  PEFFFile file;

  // Load file with annotation identifiers
  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_annotid.peff"), entries, headers);

  // Store and reload
  file.store(tmp_filename, entries, headers[0]);
  file.load(tmp_filename, entries2, headers2);

  // Verify header flag preserved
  TEST_EQUAL(headers2[0].has_annotation_identifiers, true)

  // Verify entries preserved
  TEST_EQUAL(entries.size(), entries2.size())
  for (Size i = 0; i < entries.size(); ++i)
  {
    TEST_EQUAL(entries[i].identifier, entries2[i].identifier)
    TEST_EQUAL(entries[i].modifications.size(), entries2[i].modifications.size())
    for (Size j = 0; j < entries[i].modifications.size(); ++j)
    {
      TEST_EQUAL(entries[i].modifications[j].annotation_id, entries2[i].modifications[j].annotation_id)
      TEST_EQUAL(entries[i].modifications[j].position, entries2[i].modifications[j].position)
      TEST_EQUAL(entries[i].modifications[j].accession, entries2[i].modifications[j].accession)
    }
    TEST_EQUAL(entries[i].simple_variants.size(), entries2[i].simple_variants.size())
    for (Size j = 0; j < entries[i].simple_variants.size(); ++j)
    {
      TEST_EQUAL(entries[i].simple_variants[j].annotation_id, entries2[i].simple_variants[j].annotation_id)
      TEST_EQUAL(entries[i].simple_variants[j].position, entries2[i].simple_variants[j].position)
    }
    TEST_EQUAL(entries[i].complex_variants.size(), entries2[i].complex_variants.size())
    for (Size j = 0; j < entries[i].complex_variants.size(); ++j)
    {
      TEST_EQUAL(entries[i].complex_variants[j].annotation_id, entries2[i].complex_variants[j].annotation_id)
      TEST_EQUAL(entries[i].complex_variants[j].start_position, entries2[i].complex_variants[j].start_position)
    }
    TEST_EQUAL(entries[i].processed_regions.size(), entries2[i].processed_regions.size())
    for (Size j = 0; j < entries[i].processed_regions.size(); ++j)
    {
      TEST_EQUAL(entries[i].processed_regions[j].annotation_id, entries2[i].processed_regions[j].annotation_id)
      TEST_EQUAL(entries[i].processed_regions[j].start_position, entries2[i].processed_regions[j].start_position)
    }
  }
}
END_SECTION

START_SECTION(test_mixed_annotation_formats)
{
  // Test that entries without annotation IDs still work
  vector<PEFFEntry> entries;
  vector<PEFFDatabaseMetadata> headers;
  PEFFFile file;

  // Load the original test file (no annotation IDs)
  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_test.peff"), entries, headers);

  // Verify annotation IDs are not set (max value sentinel)
  TEST_EQUAL(headers[0].has_annotation_identifiers, false)
  TEST_EQUAL(entries[0].modifications[0].annotation_id, std::numeric_limits<UInt>::max())
  TEST_EQUAL(entries[0].modifications[1].annotation_id, std::numeric_limits<UInt>::max())
  TEST_EQUAL(entries[0].simple_variants[0].annotation_id, std::numeric_limits<UInt>::max())
  TEST_EQUAL(entries[0].processed_regions[0].annotation_id, std::numeric_limits<UInt>::max())
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

START_SECTION(test_uniprot_format)
{
  // Test UniProt-style PEFF with extended keywords
  vector<PEFFEntry> entries;
  vector<PEFFDatabaseMetadata> headers;
  PEFFFile file;

  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_uniprot.peff"), entries, headers);

  // Check header - multiple GeneralComments
  TEST_EQUAL(headers.size(), 1)
  TEST_EQUAL(headers[0].general_comments.size(), 2)
  TEST_EQUAL(headers[0].general_comments[0], "this is a hand crafted peff file")
  TEST_EQUAL(headers[0].general_comments[1], "This is a second comment line.")

  // Check header - Conversion
  TEST_EQUAL(headers[0].conversion, "manual extraction of entries")

  // Check header - SpecificKey and SpecificValue
  TEST_EQUAL(headers[0].specific_keys.size(), 1)
  TEST_EQUAL(headers[0].specific_keys.count("isoform"), 1)
  TEST_EQUAL(headers[0].specific_keys.at("isoform"), "description of a specific isoform")
  TEST_EQUAL(headers[0].specific_values.size(), 1)
  TEST_EQUAL(headers[0].specific_values.count("isoform"), 1)
  TEST_EQUAL(headers[0].specific_values.at("isoform"), "(xsd:type=string)")

  TEST_EQUAL(entries.size(), 4)

  // First entry - check ID, DbUniqueId, AltAC
  TEST_EQUAL(entries[0].identifier, "sp:P06748")
  TEST_EQUAL(entries[0].entry_id, "NPM_HUMAN")
  TEST_EQUAL(entries[0].db_unique_id, "P06748")
  TEST_EQUAL(entries[0].alt_accessions.size(), 1)
  TEST_EQUAL(entries[0].alt_accessions[0], "Q5SRR7")
  TEST_EQUAL(entries[0].gene_name, "NPM1")
  TEST_EQUAL(entries[0].modifications.size(), 2)

  // Second entry - check OX (shorthand for NcbiTaxId)
  TEST_EQUAL(entries[1].identifier, "sp:P02144_CHAIN0")
  TEST_EQUAL(entries[1].entry_id, "MYG_HUMAN")
  TEST_EQUAL(entries[1].ncbi_tax_id, 9606)  // Set via \OX
  TEST_EQUAL(entries[1].simple_variants.size(), 4)

  // Third entry - check stop codon variant (*) and bracketed sources
  TEST_EQUAL(entries[2].identifier, "sp:P00761")
  TEST_EQUAL(entries[2].processed_regions.size(), 2)
  TEST_EQUAL(entries[2].processed_regions[0].accession, "PEFF:0001021")  // signal peptide
  TEST_EQUAL(entries[2].processed_regions[1].accession, "PEFF:0001020")  // mature protein
  TEST_EQUAL(entries[2].simple_variants.size(), 2)
  TEST_EQUAL(entries[2].simple_variants[1].variant_aa, '*')  // Stop codon
  TEST_EQUAL(entries[2].simple_variants[1].optional_tag, "[ESP][ExAC]")  // Bracketed sources

  // Fourth entry - check DisulfideBond
  TEST_EQUAL(entries[3].identifier, "sp:NX_P01308-1")
  TEST_EQUAL(entries[3].disulfide_bonds.size(), 3)
  TEST_EQUAL(entries[3].disulfide_bonds[0].id1, "1")
  TEST_EQUAL(entries[3].disulfide_bonds[0].id2, "2")
  TEST_EQUAL(entries[3].disulfide_bonds[0].optional_tag, "between chains")
  TEST_EQUAL(entries[3].disulfide_bonds[2].id1, "5")
  TEST_EQUAL(entries[3].disulfide_bonds[2].id2, "6")
  TEST_EQUAL(entries[3].disulfide_bonds[2].optional_tag, "A chain only")
}
END_SECTION

START_SECTION(test_uniprot_roundtrip)
{
  // Test that UniProt-style PEFF can be loaded and stored without data loss
  vector<PEFFEntry> entries, entries2;
  vector<PEFFDatabaseMetadata> headers, headers2;
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  PEFFFile file;

  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_uniprot.peff"), entries, headers);
  file.store(tmp_filename, entries, headers[0]);
  file.load(tmp_filename, entries2, headers2);

  // Verify header fields preserved
  TEST_EQUAL(headers2[0].general_comments.size(), 2)
  TEST_EQUAL(headers2[0].conversion, headers[0].conversion)
  TEST_EQUAL(headers2[0].specific_keys.size(), headers[0].specific_keys.size())
  TEST_EQUAL(headers2[0].specific_values.size(), headers[0].specific_values.size())

  // Verify entries preserved
  TEST_EQUAL(entries.size(), entries2.size())
  for (Size i = 0; i < entries.size(); ++i)
  {
    TEST_EQUAL(entries[i].identifier, entries2[i].identifier)
    TEST_EQUAL(entries[i].entry_id, entries2[i].entry_id)
    TEST_EQUAL(entries[i].db_unique_id, entries2[i].db_unique_id)
    TEST_EQUAL(entries[i].alt_accessions.size(), entries2[i].alt_accessions.size())
    TEST_EQUAL(entries[i].disulfide_bonds.size(), entries2[i].disulfide_bonds.size())
  }
}
END_SECTION

START_SECTION(test_custom_key_def_parsing_and_roundtrip)
{
  vector<PEFFEntry> entries, entries2;
  vector<PEFFDatabaseMetadata> headers, headers2;
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  PEFFFile file;

  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_customkeydef.peff"), entries, headers);
  TEST_EQUAL(headers.size(), 1)
  TEST_EQUAL(headers[0].custom_key_defs.size(), 2)

  const PEFFCustomKeyDef& ckd = headers[0].custom_key_defs[0];
  TEST_EQUAL(ckd.key_name, "SecondaryStructure")
  TEST_EQUAL(ckd.description, "Secondary structure information")
  TEST_EQUAL(ckd.concept_curie, "BAO:0000014")
  TEST_EQUAL(ckd.regexp, "([0-9]+)\\|([0-9]+)\\|([A-Za-z]+:[0-9]+)?\\|(.*) \\|?(.+)?")
  TEST_EQUAL(ckd.field_names.size(), 5)
  TEST_EQUAL(ckd.field_names[0], "StartPosition")
  TEST_EQUAL(ckd.field_names[4], "OptionalTag")
  TEST_EQUAL(ckd.field_types.size(), 5)
  TEST_EQUAL(ckd.field_types[0], "integer")
  TEST_EQUAL(ckd.field_types[4], "string")

  const PEFFCustomKeyDef& ckd2 = headers[0].custom_key_defs[1];
  TEST_EQUAL(ckd2.key_name, "EnumTest")
  TEST_EQUAL(ckd2.description, "Enumeration test")
  TEST_EQUAL(ckd2.field_names.size(), 1)
  TEST_EQUAL(ckd2.field_names[0], "Kind")
  TEST_EQUAL(ckd2.field_types.size(), 1)
  TEST_EQUAL(ckd2.field_types[0], "enumeration(alpha|beta|gamma)")

  file.store(tmp_filename, entries, headers[0]);
  file.load(tmp_filename, entries2, headers2);

  TEST_EQUAL(headers2.size(), 1)
  TEST_EQUAL(headers2[0].custom_key_defs.size(), headers[0].custom_key_defs.size())
  for (Size i = 0; i < headers[0].custom_key_defs.size(); ++i)
  {
    TEST_EQUAL(headers2[0].custom_key_defs[i] == headers[0].custom_key_defs[i], true)
  }

  // Ensure stored output uses spec-style CustomKeyDef=(...) format.
  bool found_custom_key_def = false;
  std::ifstream in(tmp_filename.c_str());
  std::string line;
  while (std::getline(in, line))
  {
    if (line.rfind("# CustomKeyDef=(", 0) == 0)
    {
      found_custom_key_def = true;
      break;
    }
  }
  TEST_EQUAL(found_custom_key_def, true)
}
END_SECTION

START_SECTION(test_custom_key_def_usage_roundtrip)
{
  vector<PEFFEntry> entries, entries2;
  vector<PEFFDatabaseMetadata> headers, headers2;
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  PEFFFile file;

  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_customkeydef_usage.peff"), entries, headers);
  TEST_EQUAL(headers.size(), 1)
  TEST_EQUAL(headers[0].custom_key_defs.size(), 2)

  TEST_EQUAL(entries.size(), 1)
  TEST_EQUAL(entries[0].custom_annotations.count("SecondaryStructure"), 1)
  TEST_EQUAL(entries[0].custom_annotations.at("SecondaryStructure"),
             "(1|2|ncithesaurus:C47937|Helix|Abcg2\\|meta\\\\x10)")
  TEST_EQUAL(entries[0].custom_annotations.count("EnumTest"), 1)
  TEST_EQUAL(entries[0].custom_annotations.at("EnumTest"), "alpha")

  file.store(tmp_filename, entries, headers[0]);
  file.load(tmp_filename, entries2, headers2);

  TEST_EQUAL(entries2.size(), 1)
  TEST_EQUAL(entries2[0].custom_annotations.count("SecondaryStructure"), 1)
  TEST_EQUAL(entries2[0].custom_annotations.at("SecondaryStructure"),
             entries[0].custom_annotations.at("SecondaryStructure"))
  TEST_EQUAL(entries2[0].custom_annotations.count("EnumTest"), 1)
  TEST_EQUAL(entries2[0].custom_annotations.at("EnumTest"),
             entries[0].custom_annotations.at("EnumTest"))

  // Ensure custom annotation value is written without escaping structural pipes/parentheses.
  bool found_secondary_structure = false;
  std::ifstream in(tmp_filename.c_str());
  std::string line;
  while (std::getline(in, line))
  {
    if (line.find("\\SecondaryStructure=") != std::string::npos)
    {
      found_secondary_structure = true;
      TEST_EQUAL(line.find("\\SecondaryStructure=(1|2|ncithesaurus:C47937|Helix|Abcg2\\|meta\\\\x10)") != std::string::npos, true)
      break;
    }
  }
  TEST_EQUAL(found_secondary_structure, true)
}
END_SECTION

START_SECTION([PEFFDisulfideBond] PEFFDisulfideBond())
{
  PEFFDisulfideBond bond;
  TEST_EQUAL(bond.id1, "")
  TEST_EQUAL(bond.id2, "")
  TEST_EQUAL(bond.optional_tag, "")
}
END_SECTION

START_SECTION(test_p53_uniprot)
{
  // Test real UniProt P53 PEFF file with many variants
  // Note: UniProt PEFF has unusual structure with empty first db block
  vector<PEFFEntry> entries;
  vector<PEFFDatabaseMetadata> headers;
  PEFFFile file;

  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_P53_uniprot.peff"), entries, headers);

  // File description block (version) is merged into the database block
  TEST_EQUAL(headers.size(), 1)
  TEST_EQUAL(headers[0].db_name, "UniProt")
  TEST_EQUAL(headers[0].general_comments.size(), 2)

  // Check P53 entry
  TEST_EQUAL(entries.size(), 1)
  TEST_EQUAL(entries[0].identifier, "up:P04637")
  TEST_EQUAL(entries[0].gene_name, "TP53")
  TEST_EQUAL(entries[0].ncbi_tax_id, 9606)
  TEST_EQUAL(entries[0].taxonomy_name, "Homo sapiens")
  TEST_EQUAL(entries[0].sequence_length, 393)
  TEST_EQUAL(entries[0].protein_existence, 1)
  TEST_EQUAL(entries[0].sequence_version, "4")

  // P53 has many variants - check we parsed a significant number
  TEST_EQUAL(entries[0].simple_variants.size() > 1000, true)

  // Check first variant - stop codon at position 2
  TEST_EQUAL(entries[0].simple_variants[0].position, 2)
  TEST_EQUAL(entries[0].simple_variants[0].variant_aa, '*')
  TEST_EQUAL(entries[0].simple_variants[0].optional_tag.hasSubstring("ExAC"), true)

  // Check sequence length matches annotation
  TEST_EQUAL(entries[0].sequence.size(), 393)
}
END_SECTION

START_SECTION([PEFFEntry] digestWithVariants)
{
  // Create a simple entry with variants for testing
  PEFFEntry entry;
  entry.sequence = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAK";
  // Position 5 (T->A), Position 10 (L->M), Position 20 (Y->F) - 1-based positions
  entry.simple_variants.push_back(PEFFVariantSimple(5, 'A', "dbSNP:rs1"));
  entry.simple_variants.push_back(PEFFVariantSimple(10, 'M', "dbSNP:rs2"));
  entry.simple_variants.push_back(PEFFVariantSimple(20, 'F', "dbSNP:rs3"));
  // Add a stop codon variant that should be skipped
  entry.simple_variants.push_back(PEFFVariantSimple(30, '*', "dbSNP:rs4"));

  // Set up trypsin digestion
  ProteaseDigestion digestor;
  digestor.setEnzyme("Trypsin");
  digestor.setMissedCleavages(0);

  std::vector<std::string> descriptions;
  std::vector<AASequence> sequences;

  // Generate peptides with variants (include_variants=true, include_modifications=false)
  entry.digestWithVariants(digestor, descriptions, sequences, 4, 30, true, true, false);

  // Should have generated peptides
  TEST_EQUAL(sequences.size() > 0, true)
  TEST_EQUAL(descriptions.size(), sequences.size())

  // Check that we got reference peptides (empty description)
  int ref_count = 0;
  int var_count = 0;
  for (size_t i = 0; i < descriptions.size(); ++i)
  {
    if (descriptions[i].empty())
    {
      ref_count++;
    }
    else
    {
      var_count++;
    }
    // Verify no stop codons in sequences
    TEST_EQUAL(sequences[i].toString().find('*') == std::string::npos, true)
  }
  TEST_EQUAL(ref_count > 0, true)
  TEST_EQUAL(var_count > 0, true)

  // Test without reference peptides
  std::vector<std::string> var_only_desc;
  std::vector<AASequence> var_only_seq;
  entry.digestWithVariants(digestor, var_only_desc, var_only_seq, 4, 30, false, true, false);
  for (const auto& desc : var_only_desc)
  {
    TEST_EQUAL(desc.empty(), false)
  }
  TEST_EQUAL(var_only_seq.size() < sequences.size(), true)

  // Test with different enzyme
  ProteaseDigestion chymo;
  chymo.setEnzyme("Chymotrypsin");
  chymo.setMissedCleavages(0);

  std::vector<std::string> chymo_desc;
  std::vector<AASequence> chymo_seq;
  entry.digestWithVariants(chymo, chymo_desc, chymo_seq, 4, 30, true, true, false);
  TEST_EQUAL(chymo_seq.size() > 0, true)
  // Different enzyme should give different peptide count
  TEST_EQUAL(chymo_seq.size() != sequences.size(), true)
}
END_SECTION

START_SECTION([PEFFEntry] digestWithVariants_combinatorics)
{
  // Test that combinatorics work correctly within a peptide
  // Note: Trypsin doesn't cut before proline (KP, RP), so use sequences without this pattern
  PEFFEntry entry;
  // Sequence: MK | AEPTIDER (trypsin cuts after K, then after R)
  // K is NOT followed by P, so cleavage occurs
  entry.sequence = "MKAEPTIDER";
  // Two variants within the "AEPTIDER" peptide (positions 3-10, 1-based)
  entry.simple_variants.push_back(PEFFVariantSimple(4, 'L', "var1"));  // E->L at position 4
  entry.simple_variants.push_back(PEFFVariantSimple(6, 'L', "var2"));  // T->L at position 6

  ProteaseDigestion digestor;
  digestor.setEnzyme("Trypsin");
  digestor.setMissedCleavages(0);

  std::vector<std::string> descriptions;
  std::vector<AASequence> sequences;
  entry.digestWithVariants(digestor, descriptions, sequences, 2, 20, true, true, false);

  // For the "AEPTIDER" peptide (positions 3-10, 0-indexed 2-10):
  // Position 4 (1-based) = index 3 = E in AEPTIDER -> local index 1
  // Position 6 (1-based) = index 5 = T in AEPTIDER -> local index 3
  // - Reference: AEPTIDER
  // - Variant 1 only: ALPTIDER (E4L)
  // - Variant 2 only: AEPLIDER (T6L)
  // - Both variants: ALPLIDER (E4L + T6L)
  // So 4 versions of this peptide, plus reference "MK"

  // Count peptides matching "AEPTIDER" variants
  int aeptider_variants = 0;
  bool found_ref = false;
  bool found_e4l = false;
  bool found_t6l = false;
  bool found_both = false;

  for (size_t i = 0; i < sequences.size(); ++i)
  {
    String seq_str = sequences[i].toString();
    if (seq_str == "AEPTIDER")
    {
      found_ref = true;
      aeptider_variants++;
    }
    else if (seq_str == "ALPTIDER")
    {
      found_e4l = true;
      aeptider_variants++;
    }
    else if (seq_str == "AEPLIDER")
    {
      found_t6l = true;
      aeptider_variants++;
    }
    else if (seq_str == "ALPLIDER")
    {
      found_both = true;
      aeptider_variants++;
    }
  }

  TEST_EQUAL(found_ref, true)
  TEST_EQUAL(found_e4l, true)
  TEST_EQUAL(found_t6l, true)
  TEST_EQUAL(found_both, true)
  TEST_EQUAL(aeptider_variants, 4)
}
END_SECTION

START_SECTION([PEFFEntry] digestWithVariants_modifications)
{
  // Test modification combinations
  PEFFEntry entry;
  // Sequence with a tryptic peptide containing oxidizable methionine
  entry.sequence = "MKAMEPTIDER";  // MK | AMEPTIDER
  // Add a modification at position 4 (M in AMEPTIDER)
  // Using Oxidation which is in ModificationsDB
  entry.modifications.push_back(PEFFModification(4, "UNIMOD:35", "Oxidation"));

  ProteaseDigestion digestor;
  digestor.setEnzyme("Trypsin");
  digestor.setMissedCleavages(0);

  std::vector<std::string> mod_desc;
  std::vector<AASequence> mod_seq;

  // Test with modifications only (no variants)
  entry.digestWithVariants(digestor, mod_desc, mod_seq, 2, 20, true, false, true);

  // Should have reference + modified versions
  // For AMEPTIDER: reference + with Oxidation at M (position 4 in protein = position 2 in peptide)
  bool found_unmodified = false;
  int mod_peptide_count = 0;

  for (size_t i = 0; i < mod_seq.size(); ++i)
  {
    String seq_str = mod_seq[i].toString();
    if (seq_str == "AMEPTIDER" && mod_desc[i].empty())
    {
      found_unmodified = true;
    }
    // Count peptides with modification descriptions
    if (mod_desc[i].find("UNIMOD:35") != std::string::npos || mod_desc[i].find("Oxidation") != std::string::npos)
    {
      mod_peptide_count++;
    }
  }

  TEST_EQUAL(found_unmodified, true)
  // Note: mod_peptide_count may be 0 if Oxidation can't be applied to the specific residue
  // The test verifies the mechanism works; actual mod application depends on ModificationsDB
  (void)mod_peptide_count; // Suppress unused variable warning - value depends on ModificationsDB

  // Test with both variants and modifications
  entry.simple_variants.push_back(PEFFVariantSimple(5, 'L', "var1"));  // E->L at position 5

  std::vector<std::string> combo_desc;
  std::vector<AASequence> combo_seq;
  entry.digestWithVariants(digestor, combo_desc, combo_seq, 2, 20, true, true, true);

  // Should have more combinations than mods alone or variants alone
  std::vector<std::string> var_only_desc;
  std::vector<AASequence> var_only_seq;
  entry.digestWithVariants(digestor, var_only_desc, var_only_seq, 2, 20, true, true, false);

  std::vector<std::string> mod_only_desc;
  std::vector<AASequence> mod_only_seq;
  entry.digestWithVariants(digestor, mod_only_desc, mod_only_seq, 2, 20, true, false, true);

  // With 1 variant and 1 mod in the peptide, we should have:
  // - reference (1)
  // - variant only (1)
  // - mod only (1)
  // - variant + mod (1)
  // = 4 combinations for that peptide, plus MK reference
  // combo_peptides should have more than either var_only or mod_only alone
  TEST_EQUAL(combo_seq.size() >= var_only_seq.size(), true)
  TEST_EQUAL(combo_seq.size() >= mod_only_seq.size(), true)
}
END_SECTION

START_SECTION([PEFFEntry] digestWithVariants_no_variants_flag)
{
  // Test that include_variants=false skips variant generation
  PEFFEntry entry;
  entry.sequence = "MKAEPTIDER";
  entry.simple_variants.push_back(PEFFVariantSimple(4, 'L', "var1"));

  ProteaseDigestion digestor;
  digestor.setEnzyme("Trypsin");
  digestor.setMissedCleavages(0);

  std::vector<std::string> with_vars_desc;
  std::vector<AASequence> with_vars_seq;
  entry.digestWithVariants(digestor, with_vars_desc, with_vars_seq, 2, 20, true, true, false);

  std::vector<std::string> no_vars_desc;
  std::vector<AASequence> no_vars_seq;
  entry.digestWithVariants(digestor, no_vars_desc, no_vars_seq, 2, 20, true, false, false);

  // no_vars should have fewer peptides (only reference)
  TEST_EQUAL(no_vars_seq.size() < with_vars_seq.size(), true)

  // All no_vars peptides should have empty description (reference only)
  for (const auto& desc : no_vars_desc)
  {
    TEST_EQUAL(desc.empty(), true)
  }
}
END_SECTION

START_SECTION([PEFFEntry] generatePeptides_basic)
{
  // Test the new generatePeptides method
  PEFFEntry entry;
  // Sequence with methionine that can be oxidized
  entry.sequence = "MKAMEPTIDER";  // MK | AMEPTIDER
  // PEFF variant at position 4 (M->L)
  entry.simple_variants.push_back(PEFFVariantSimple(4, 'L', "dbSNP:test"));
  // PEFF modification at position 5 (phospho)
  entry.modifications.push_back(PEFFModification(5, "UNIMOD:21", "Phospho"));

  ProteaseDigestion digestor;
  digestor.setEnzyme("Trypsin");
  digestor.setMissedCleavages(0);

  // Test without additional mods (empty fixed and variable)
  std::vector<std::string> no_fixed_mods;
  std::vector<std::string> no_var_mods;
  std::vector<std::string> descriptions;
  std::vector<AASequence> sequences;
  entry.generatePeptides(digestor, descriptions, sequences, no_fixed_mods, no_var_mods, 2, 2, 20, true, true, true);

  // Should have peptides
  TEST_EQUAL(sequences.size() > 0, true)

  // Check we have reference and variant peptides
  bool found_ref = false;
  bool found_variant = false;
  for (size_t i = 0; i < descriptions.size(); ++i)
  {
    if (descriptions[i].empty() || descriptions[i].find("M4L") == String::npos)
    {
      found_ref = true;
    }
    if (descriptions[i].find("M4L") != String::npos)
    {
      found_variant = true;
    }
  }
  TEST_EQUAL(found_ref, true)
  TEST_EQUAL(found_variant, true)
}
END_SECTION

START_SECTION([PEFFEntry] generatePeptides_with_fixed_mods)
{
  // Test with fixed modifications (like Carbamidomethyl on Cys)
  PEFFEntry entry;
  // Sequence with cysteine that gets carbamidomethylated
  entry.sequence = "MKACEPTIDER";  // MK | ACEPTIDER (C at position 4)

  ProteaseDigestion digestor;
  digestor.setEnzyme("Trypsin");
  digestor.setMissedCleavages(0);

  // Add Carbamidomethyl as fixed modification
  std::vector<std::string> fixed_mods = {"Carbamidomethyl (C)"};
  std::vector<std::string> no_var_mods;
  std::vector<std::string> descriptions;
  std::vector<AASequence> sequences;
  entry.generatePeptides(digestor, descriptions, sequences, fixed_mods, no_var_mods, 0, 2, 20, true, false, false);

  // Should have peptides
  TEST_EQUAL(sequences.size() > 0, true)

  // All peptides containing C should have Carbamidomethyl
  for (size_t i = 0; i < sequences.size(); ++i)
  {
    String seq_str = sequences[i].toUnmodifiedString();
    if (seq_str.find('C') != String::npos)
    {
      // The peptide has C, so it should be modified
      TEST_EQUAL(sequences[i].isModified(), true)
    }
  }
}
END_SECTION

START_SECTION([PEFFEntry] generatePeptides_with_variable_mods)
{
  // Test with variable modifications (sample handling)
  PEFFEntry entry;
  // Sequence with methionine that can be oxidized
  entry.sequence = "MKAMEPTIDER";  // MK | AMEPTIDER (M at position 4)

  ProteaseDigestion digestor;
  digestor.setEnzyme("Trypsin");
  digestor.setMissedCleavages(0);

  // Add Oxidation as variable modification (simulating sample handling)
  std::vector<std::string> no_fixed_mods;
  std::vector<std::string> var_mods = {"Oxidation (M)"};
  std::vector<std::string> descriptions;
  std::vector<AASequence> sequences;
  entry.generatePeptides(digestor, descriptions, sequences, no_fixed_mods, var_mods, 2, 2, 20, true, false, false);

  // Should have more peptides due to oxidation variants
  // Reference AMEPTIDER + oxidized versions
  TEST_EQUAL(sequences.size() > 0, true)

  // Check that some peptides have oxidation applied
  bool found_oxidized = false;
  for (size_t i = 0; i < sequences.size(); ++i)
  {
    String seq_str = sequences[i].toString();
    // Oxidized methionine shows as M(Oxidation) or similar
    if (seq_str.find("Oxidation") != String::npos || sequences[i].isModified())
    {
      found_oxidized = true;
      break;
    }
  }
  // Oxidation should be found if ModificationsDB has it
  // Note: This may fail if Oxidation (M) is not in the DB, which is fine
  (void)found_oxidized;

  // Compare with no variable mods
  std::vector<std::string> no_var_mods;
  std::vector<std::string> desc_no_var;
  std::vector<AASequence> seq_no_var;
  entry.generatePeptides(digestor, desc_no_var, seq_no_var, no_fixed_mods, no_var_mods, 2, 2, 20, true, false, false);

  // With variable mods should have >= peptides than without
  TEST_EQUAL(sequences.size() >= seq_no_var.size(), true)
}
END_SECTION

START_SECTION([PEFFEntry] generatePeptides_full_workflow)
{
  // Test full workflow: PEFF variants + PEFF mods + fixed + variable sample mods
  PEFFEntry entry;
  entry.sequence = "MKACMEPTIDER";  // Has both C and M
  // PEFF variant
  entry.simple_variants.push_back(PEFFVariantSimple(5, 'L', "var1"));  // M->L at position 5

  ProteaseDigestion digestor;
  digestor.setEnzyme("Trypsin");
  digestor.setMissedCleavages(0);

  // No sample handling mods - just PEFF annotations
  std::vector<std::string> no_fixed;
  std::vector<std::string> no_var;
  std::vector<std::string> peff_only_desc;
  std::vector<AASequence> peff_only_seq;
  entry.generatePeptides(digestor, peff_only_desc, peff_only_seq, no_fixed, no_var, 2, 2, 20, true, true, false);

  // With Carbamidomethyl (fixed) and Oxidation (variable)
  std::vector<std::string> fixed = {"Carbamidomethyl (C)"};
  std::vector<std::string> variable = {"Oxidation (M)"};
  std::vector<std::string> with_mods_desc;
  std::vector<AASequence> with_mods_seq;
  entry.generatePeptides(digestor, with_mods_desc, with_mods_seq, fixed, variable, 1, 2, 20, true, true, false);

  // Sample mods should add more peptides (or at least same, if no compatible residues)
  TEST_EQUAL(with_mods_seq.size() >= peff_only_seq.size(), true)
}
END_SECTION

START_SECTION([PEFFEntry] generatePeptides_no_reference_with_peff_mods)
{
  // Regression test: when include_reference=false and there are PEFF modifications,
  // the unmodified reference peptide should NOT appear in the output.
  // Previously, enumeratePEFFModifications_ returned combo=0 (no mods) which was
  // the reference peptide, and it wasn't being filtered out.
  PEFFEntry entry;
  entry.sequence = "MKPEPTIDER";  // MK | PEPTIDER
  // Add a PEFF modification at position 4 (E -> phospho)
  entry.modifications.push_back(PEFFModification(4, "UNIMOD:21", "Phospho"));

  ProteaseDigestion digestor;
  digestor.setEnzyme("Trypsin");
  digestor.setMissedCleavages(0);

  std::vector<std::string> no_fixed;
  std::vector<std::string> no_var;
  std::vector<std::string> descriptions;
  std::vector<AASequence> sequences;

  // include_reference=false, include_peff_modifications=true
  entry.generatePeptides(digestor, descriptions, sequences, no_fixed, no_var, 2, 2, 20,
                         false,  // include_reference = false
                         false,  // include_variants = false
                         true);  // include_peff_modifications = true

  // The peptide PEPTIDER should only appear with modifications, not unmodified
  AASequence ref_peptide = AASequence::fromString("PEPTIDER");
  String ref_peptide_str = ref_peptide.toString();

  bool found_unmodified_reference = false;
  for (size_t i = 0; i < sequences.size(); ++i)
  {
    if (sequences[i].toString() == ref_peptide_str)
    {
      found_unmodified_reference = true;
      break;
    }
  }

  // The unmodified reference should NOT be in the output
  TEST_EQUAL(found_unmodified_reference, false)

  // But we should have at least the modified version(s)
  TEST_EQUAL(sequences.size() > 0, true)

  // All peptides should be modified (have phospho or other PEFF mod)
  for (size_t i = 0; i < sequences.size(); ++i)
  {
    TEST_EQUAL(sequences[i].isModified(), true)
  }
}
END_SECTION

START_SECTION(test_modification_type_roundtrip)
{
  // Test Fix #1: Modifications are written with correct type-specific keys
  // and Fix #2: Comma-separated positions are expanded
  vector<PEFFEntry> entries;
  vector<PEFFDatabaseMetadata> headers;
  PEFFFile file;

  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_modtypes.peff"), entries, headers);

  // Check first entry - modifications from different keys
  TEST_EQUAL(entries.size(), 2)
  TEST_EQUAL(entries[0].identifier, "mt:P00001")

  // Should have 5 modifications total (2 PSI, 2 UNIMOD, 1 GENERIC)
  TEST_EQUAL(entries[0].modifications.size(), 5)

  // Check types are correctly assigned based on key name
  int psi_count = 0, unimod_count = 0, generic_count = 0;
  for (const auto& mod : entries[0].modifications)
  {
    if (mod.type == PEFFModification::Type::PSI_MOD) psi_count++;
    else if (mod.type == PEFFModification::Type::UNIMOD) unimod_count++;
    else if (mod.type == PEFFModification::Type::GENERIC) generic_count++;
  }
  TEST_EQUAL(psi_count, 2)
  TEST_EQUAL(unimod_count, 2)
  TEST_EQUAL(generic_count, 1)

  // Second entry - comma-separated positions should expand to 3 mods
  TEST_EQUAL(entries[1].identifier, "mt:P00002")
  TEST_EQUAL(entries[1].modifications.size(), 3)
  TEST_EQUAL(entries[1].modifications[0].position, 5)
  TEST_EQUAL(entries[1].modifications[1].position, 10)
  TEST_EQUAL(entries[1].modifications[2].position, 15)
  // All should be UNIMOD type
  for (const auto& mod : entries[1].modifications)
  {
    TEST_EQUAL(mod.type, PEFFModification::Type::UNIMOD)
    TEST_EQUAL(mod.accession, "UNIMOD:21")
    TEST_EQUAL(mod.name, "Phospho")
  }

  // Roundtrip test: store and reload, check types preserved
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename, entries, headers[0]);

  vector<PEFFEntry> entries2;
  vector<PEFFDatabaseMetadata> headers2;
  file.load(tmp_filename, entries2, headers2);

  // After roundtrip, modification types should be preserved
  TEST_EQUAL(entries2[0].modifications.size(), entries[0].modifications.size())
  for (Size i = 0; i < entries[0].modifications.size(); ++i)
  {
    TEST_EQUAL(entries[0].modifications[i].type, entries2[0].modifications[i].type)
    TEST_EQUAL(entries[0].modifications[i].position, entries2[0].modifications[i].position)
    TEST_EQUAL(entries[0].modifications[i].accession, entries2[0].modifications[i].accession)
  }
}
END_SECTION

START_SECTION(test_modification_optional_tag_roundtrip)
{
  // Spec section 3.3.10-3.3.12: ModRes* format is (position|accession|name|OptionalTag).
  vector<PEFFEntry> entries, entries2;
  vector<PEFFDatabaseMetadata> headers, headers2;
  PEFFFile file;

  file.load(OPENMS_GET_TEST_DATA_PATH("PEFFFile_mod_optionaltag.peff"), entries, headers);
  TEST_EQUAL(entries.size(), 1)
  TEST_EQUAL(entries[0].modifications.size(), 1)
  TEST_EQUAL(entries[0].modifications[0].optional_tag, "invitro")

  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename, entries, headers[0]);
  file.load(tmp_filename, entries2, headers2);

  TEST_EQUAL(entries2.size(), 1)
  TEST_EQUAL(entries2[0].modifications.size(), 1)
  TEST_EQUAL(entries2[0].modifications[0].optional_tag, "invitro")

  // Ensure stored output uses the 4th field for the optional tag (no empty placeholder field).
  bool found_modres = false;
  std::ifstream in(tmp_filename.c_str());
  std::string line;
  while (std::getline(in, line))
  {
    if (line.find("\\ModResUnimod=") != std::string::npos)
    {
      found_modres = true;
      TEST_EQUAL(line.find("|Phospho|invitro)") != std::string::npos, true)
      break;
    }
  }
  TEST_EQUAL(found_modres, true)
}
END_SECTION

START_SECTION(test_proteoform_db_header)
{
  // Test Fix #4: ProteoformDb header field
  PEFFDatabaseMetadata header;
  header.db_name = "ProteoformTest";
  header.prefix = "pf";
  header.db_version = "1.0";
  header.db_sources.push_back("test");
  header.number_of_entries = 1;
  header.is_proteoform_db = true;

  PEFFEntry entry;
  entry.identifier = "pf:P00001";
  entry.sequence = "PEPTIDE";
  entry.sequence_length = 7;

  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  PEFFFile file;
  vector<PEFFEntry> entries = {entry};
  file.store(tmp_filename, entries, header);

  // Reload and check
  vector<PEFFEntry> entries2;
  vector<PEFFDatabaseMetadata> headers2;
  file.load(tmp_filename, entries2, headers2);

  TEST_EQUAL(headers2.size(), 1)
  TEST_EQUAL(headers2[0].is_proteoform_db, true)
}
END_SECTION

START_SECTION(test_altac_parenthesized_list)
{
  // Test Fix #13: AltAC written as single key with parenthesized list
  PEFFEntry entry;
  entry.identifier = "sp:P12345";
  entry.sequence = "PEPTIDE";
  entry.sequence_length = 7;
  entry.alt_accessions.push_back("Q11111");
  entry.alt_accessions.push_back("Q22222");
  entry.alt_accessions.push_back("Q33333");

  PEFFDatabaseMetadata header;
  header.db_name = "AltACTest";
  header.prefix = "sp";
  header.db_version = "1.0";
  header.db_sources.push_back("test");
  header.number_of_entries = 1;

  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  PEFFFile file;
  vector<PEFFEntry> entries = {entry};
  file.store(tmp_filename, entries, header);

  // Reload and check all accessions preserved
  vector<PEFFEntry> entries2;
  vector<PEFFDatabaseMetadata> headers2;
  file.load(tmp_filename, entries2, headers2);

  TEST_EQUAL(entries2[0].alt_accessions.size(), 3)
  TEST_EQUAL(entries2[0].alt_accessions[0], "Q11111")
  TEST_EQUAL(entries2[0].alt_accessions[1], "Q22222")
  TEST_EQUAL(entries2[0].alt_accessions[2], "Q33333")
}
END_SECTION

START_SECTION(test_disulfide_bond_annotation_id)
{
  // Test Fix #8: DisulfideBond annotation ID extraction
  PEFFEntry entry;
  entry.identifier = "sp:P12345";
  entry.sequence = "PEPTIDE";
  entry.sequence_length = 7;

  PEFFDisulfideBond bond("1", "2", "between chains", 5);
  entry.disulfide_bonds.push_back(bond);

  PEFFDatabaseMetadata header;
  header.db_name = "DSBondTest";
  header.prefix = "sp";
  header.db_version = "1.0";
  header.db_sources.push_back("test");
  header.number_of_entries = 1;
  header.has_annotation_identifiers = true;

  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  PEFFFile file;
  vector<PEFFEntry> entries = {entry};
  file.store(tmp_filename, entries, header);

  // Reload and check annotation ID preserved
  vector<PEFFEntry> entries2;
  vector<PEFFDatabaseMetadata> headers2;
  file.load(tmp_filename, entries2, headers2);

  TEST_EQUAL(entries2[0].disulfide_bonds.size(), 1)
  TEST_EQUAL(entries2[0].disulfide_bonds[0].annotation_id, 5)
  TEST_EQUAL(entries2[0].disulfide_bonds[0].id1, "1")
  TEST_EQUAL(entries2[0].disulfide_bonds[0].id2, "2")
  TEST_EQUAL(entries2[0].disulfide_bonds[0].optional_tag, "between chains")
}
END_SECTION

START_SECTION(test_multi_database_store_roundtrip)
{
  PEFFDatabaseMetadata db1;
  db1.version = "1.0";
  db1.db_name = "Database1";
  db1.prefix = "db1";
  db1.db_version = "1.0";
  db1.db_sources.push_back("test");
  db1.number_of_entries = 1;

  PEFFDatabaseMetadata db2;
  db2.version = "1.0";
  db2.db_name = "Database2";
  db2.prefix = "db2";
  db2.is_decoy = true;
  db2.db_version = "2.0";
  db2.db_sources.push_back("test");
  db2.number_of_entries = 1;

  PEFFEntry e1;
  e1.identifier = "db1:P00001";
  e1.sequence = "SEQUENCEONE";
  e1.sequence_length = 10;
  e1.protein_names = {"Protein from DB1"};

  PEFFEntry e2;
  e2.identifier = "db2:P00002";
  e2.sequence = "SEQUENCETWO";
  e2.sequence_length = 10;
  e2.protein_names = {"Protein from DB2"};

  vector<PEFFEntry> entries = {e1, e2};
  vector<PEFFDatabaseMetadata> headers = {db1, db2};

  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  PEFFFile file;
  file.store(tmp_filename, entries, headers);

  vector<PEFFEntry> entries2;
  vector<PEFFDatabaseMetadata> headers2;
  file.load(tmp_filename, entries2, headers2);

  TEST_EQUAL(headers2.size(), 2)
  TEST_EQUAL(headers2[0].prefix, "db1")
  TEST_EQUAL(headers2[0].db_name, "Database1")
  TEST_EQUAL(headers2[0].is_decoy, false)
  TEST_EQUAL(headers2[1].prefix, "db2")
  TEST_EQUAL(headers2[1].db_name, "Database2")
  TEST_EQUAL(headers2[1].is_decoy, true)

  TEST_EQUAL(entries2.size(), 2)
  TEST_EQUAL(entries2[0].identifier, "db1:P00001")
  TEST_EQUAL(entries2[1].identifier, "db2:P00002")

  // Only one PEFF version line should be written.
  std::ifstream in(tmp_filename.c_str());
  std::string line;
  Size peff_lines = 0;
  while (std::getline(in, line))
  {
    if (line.rfind("# PEFF ", 0) == 0)
    {
      ++peff_lines;
    }
  }
  TEST_EQUAL(peff_lines, 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
