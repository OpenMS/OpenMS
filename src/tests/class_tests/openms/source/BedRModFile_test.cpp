// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sam Wein $
// $Authors: Sam Wein $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/BedRModFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <optional>
#include <tuple>

using namespace OpenMS;
using namespace std;
using ID = IdentificationData;

namespace
{
  /// An oligonucleotide run as NucleicAcidSearchEngine creates it: one input file, one database and its sequences,
  /// a hyperscore and (unless left out) a PSM-level q-value
  struct RnaRun
  {
    explicit RnaRun(bool with_qvalue = true) :
      run(data.addRun("NASE", ID::MoleculeKind::OLIGONUCLEOTIDE))
    {
      ID::ScoreDefinition hyperscore;
      hyperscore.name = "hyperscore";
      run.setPrimaryScore(run.addScore(hyperscore));
      if (with_qvalue)
      {
        ID::ScoreDefinition qvalue;
        qvalue.name = "PSM-level q-value";
        qvalue.higher_better = false;
        run.addScore(qvalue);
      }
      ID::Database fasta;
      fasta.path = "test.fasta";
      database = run.addDatabase(fasta);
      ID::SourceFile file;
      file.path = "test.mzML";
      source = run.addSource(file);
    }
    RnaRun(const RnaRun&) = delete;

    /// Add an entry of the database (a parent sequence)
    void addSequence(const std::string& accession, const std::string& sequence, bool decoy = false)
    {
      ID::DatabaseSequence entry(database, accession, decoy ? ID::TargetDecoy::DECOY : ID::TargetDecoy::TARGET);
      entry.sequence = sequence;
      sequences.push_back(entry);
      run.setDatabaseSequences(sequences);
    }

    /// Add a spectrum of the input file
    ID::QueryId addObservation(const std::string& data_id, double rt = 100.0, double mz = 500.0)
    {
      ID::Observation observation;
      observation.data_id = data_id;
      observation.rt = rt;
      observation.mz = mz;
      return run.addIdentification(source, observation);
    }

    /// Add a match of @p sequence to @p query that occurs at @p evidence (accession, start, end), with a hyperscore of
    /// 100 and the PSM-level q-value @p qvalue (if the run has one)
    void addMatch(ID::QueryId query, const std::string& sequence,
                  const std::vector<std::tuple<std::string, UInt32, UInt32>>& evidence,
                  std::optional<double> qvalue = 0.01, Int charge = 2)
    {
      ID::MatchData match;
      match.encoding = ID::Encoding::NA_SEQUENCE;
      match.representation = sequence;
      match.charge = charge;
      for (const auto& [accession, start, end] : evidence)
      {
        match.sequence_evidence.push_back({database, accession, start, end});
      }
      std::vector<std::optional<double>> scores {100.0};
      if (run.getScoreDefinitions().size() > 1) scores.push_back(qvalue);
      run.addMatch(query, match, scores);
    }

    ID data;
    ID::Run& run;
    ID::DatabaseId database;
    ID::SourceId source;
    std::vector<ID::DatabaseSequence> sequences;
  };
}

///////////////////////////

START_TEST(BedRModFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BedRModFile* ptr = nullptr;
BedRModFile* nullPointer = nullptr;

START_SECTION((BedRModFile()))
{
  ptr = new BedRModFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~BedRModFile()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void store(const std::string& out_file, const IdentificationData& id_data, const std::string& chebi_mapping_file, const std::map<std::string, std::string>& reference_metadata)))
{
  // Create test identification data with a parent sequence (RNA transcript)
  RnaRun rna;
  rna.addSequence("test_rna", "AUG[m5C]AUGC");

  // A sequence with both unmodified and modified bases (AUG[m5C]) at positions 0-3 of the parent, in one spectrum,
  // with a hyperscore and a PSM-level q-value
  rna.addMatch(rna.addObservation("spectrum_1"), "AUG[m5C]", {{"test_rna", 0, 3}});
  const IdentificationData& id_data = rna.data;

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Count data lines (non-header, non-empty)
  int data_lines = 0;
  bool found_mod_names_header = false;

  for (const auto& line : lines)
  {
    if (line.empty())
    {
      continue;
    }
    if (line.starts_with("#modification_names="))
    {
      found_mod_names_header = (line == "#modification_names=m5C:m5C:C");
    }
    if (line[0] == '#')
    {
      continue;
    }
    data_lines++;
  }

  // Only the modified residue produces a row.
  TEST_EQUAL(data_lines, 1)
  TEST_TRUE(found_mod_names_header) // Should list m5C in modification_names header
}
END_SECTION

START_SECTION(([EXTRA] distinct observations determine modification frequency and identity))
{
  RnaRun rna;
  rna.addSequence("rna", "ACGA");
  ID::Observation intense;
  intense.data_id = "modified";
  intense.setMetaValue("precursor_intensity", 1e12);
  const auto modified = rna.run.addIdentification(rna.source, intense);
  const auto unmodified = rna.addObservation("unmodified");
  const auto add_match = [&](const std::string& sequence, ID::QueryId observation, UInt32 start, double score = 0.01)
  {
    rna.addMatch(observation, sequence, {{"rna", start, start + 2}}, score);
  };
  const IdentificationData& data = rna.data;
  // Both hypotheses from one spectrum put m5C at the same parent position.
  add_match("A[m5C]G", modified, 0, 0.02);
  add_match("[m5C]GA", modified, 1);
  add_match("ACG", unmodified, 0);

  std::string output_file;
  NEW_TMP_FILE(output_file)
  BedRModFile().store(output_file, data);
  TextFile output(output_file);
  TEST_TRUE(std::find(output.begin(), output.end(), "#modification_names=m5C:m5C:C") != output.end())
  Size row_count = 0;
  for (const auto& line : output)
  {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> fields;
    StringUtils::split(line, '\t', fields);
    TEST_EQUAL(fields.size(), 14)
    if (fields.size() == 14)
    {
      TEST_EQUAL(fields[1], "1")
      TEST_EQUAL(fields[2], "2")
      TEST_EQUAL(fields[3], "m5C")
      TEST_EQUAL(fields[4], "0.01")
      TEST_EQUAL(fields[9], "2") // observations, independent of precursor intensity
      TEST_EQUAL(fields[10], "50")
      TEST_EQUAL(fields[12], "1,2")
      TEST_EQUAL(fields[13], "3,4")
    }
    ++row_count;
  }
  TEST_EQUAL(row_count, 1)

  // A supplied ChEBI mapping uses the same identifier in row and dictionary.
  std::string mapping_file;
  NEW_TMP_FILE(mapping_file)
  {
    std::ofstream mapping(mapping_file);
    mapping << "mod,chebi_id\nm5C,20607\n";
  }
  BedRModFile().store(output_file, data, mapping_file);
  TextFile mapped(output_file);
  TEST_TRUE(std::find(mapped.begin(), mapped.end(), "#modification_names=20607:m5C:C") != mapped.end())
  for (const auto& line : mapped)
  {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> fields;
    StringUtils::split(line, '\t', fields);
    TEST_EQUAL(fields.size(), 14)
    if (fields.size() == 14)
    {
      TEST_EQUAL(fields[3], "20607")
    }
  }
}
END_SECTION

START_SECTION(([EXTRA] terminal modifications excluded))
{
  // Test that terminal modifications (5' and 3') are excluded
  RnaRun rna;
  rna.addSequence("test_rna_terminal", "AUGCAUGC");

  // Only the internal modification produces a row, not either terminal group.
  rna.addMatch(rna.addObservation("spectrum_1"), "pA[Um]p", {{"test_rna_terminal", 0, 1}});
  const IdentificationData& id_data = rna.data;

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Count data lines for the internal Um modification.
  int data_lines = 0;
  for (const auto& line : lines)
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    data_lines++;
  }

  TEST_EQUAL(data_lines, 1)
}
END_SECTION

START_SECTION(([EXTRA] matches without q-value are skipped))
{
  // Test that matches without a valid q-value score are excluded from output
  // A run with only a hyperscore: deliberately no q-value score, to simulate "FDR not run"
  RnaRun rna(false);
  rna.addSequence("test_rna_no_fdr", "AUGC");
  rna.addMatch(rna.addObservation("spectrum_1"), "AU[Gm]", {{"test_rna_no_fdr", 0, 2}});
  const IdentificationData& id_data = rna.data;

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Count data lines
  int data_lines = 0;
  for (const auto& line : lines)
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    data_lines++;
  }

  // Should have 0 data lines since no q-value was provided
  TEST_EQUAL(data_lines, 0)
}
END_SECTION

START_SECTION(([EXTRA] sentinel q-values do not produce rows or affect frequency))
{
  RnaRun rna;
  rna.addSequence("test_rna_sentinel", "ACG");
  rna.addMatch(rna.addObservation("valid", 100.0), "A[m5C]G", {{"test_rna_sentinel", 0, 2}}, 0.01);
  rna.addMatch(rna.addObservation("sentinel", 101.0), "ACG", {{"test_rna_sentinel", 0, 2}}, -1.0);
  const IdentificationData& id_data = rna.data;

  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile().store(test_file, id_data);

  TextFile output;
  output.load(test_file);
  int data_lines = 0;
  for (const auto& line : output)
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    vector<std::string> fields;
    StringUtils::split(line, '\t', fields);
    TEST_TRUE(fields.size() >= 11)
    if (fields.size() >= 11)
    {
      TEST_EQUAL(fields[4], "0.01")
      TEST_EQUAL(fields[10], "100")
    }
    ++data_lines;
  }
  TEST_EQUAL(data_lines, 1)
}
END_SECTION

START_SECTION(([EXTRA] target_mapping_count reflects target occurrences))
{
  // Test that target_mapping_count column shows count of target sequence occurrences
  RnaRun rna;

  // Add three parent sequences - 2 targets and 1 decoy
  rna.addSequence("target_rna_1", "AUGCAUGC");
  rna.addSequence("target_rna_2", "CAUGCAUG");
  rna.addSequence("DECOY_rna", "GCAUGCAU", true);

  // A sequence that maps to multiple locations of all three parents
  rna.addMatch(rna.addObservation("spectrum_1"), "AU[Gm]",
               {{"target_rna_1", 0, 2},  // target 1
                {"target_rna_1", 4, 6},  // target 1 (second occurrence)
                {"target_rna_2", 1, 3},  // target 2
                {"DECOY_rna", 2, 4}});   // decoy (should not count)
  const IdentificationData& id_data = rna.data;

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Check that target_mapping_count (column 12) equals 3 for all data lines
  // (2 occurrences in target 1 + 1 occurrence in target 2 = 3 total target occurrences)
  // The decoy occurrence should NOT be counted
  int data_lines_checked = 0;
  for (const auto& line : lines)
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    TEST_FALSE(line.starts_with("DECOY_"))
    
    // Parse the line (tab-separated)
    vector<std::string> fields;
    StringUtils::split(line, '\t', fields);
    
    // Check that we have enough fields
    TEST_TRUE(fields.size() >= 12)
    
    if (fields.size() >= 12)
    {
      // Column 12 (0-indexed: 11) is unique_mapping/target_mapping_count
      Int mapping_count = StringUtils::toInt32(fields[11]);
      
      // Should be 3 (only target occurrences)
      TEST_EQUAL(mapping_count, 3)
      data_lines_checked++;
    }
  }

  // Should have checked some data lines
  TEST_TRUE(data_lines_checked > 0)
}
END_SECTION

START_SECTION(([EXTRA] target_mapping_count equals 1 for unique mapping))
{
  // Test that target_mapping_count is 1 for a sequence with single target occurrence
  // One target parent sequence and a single match to it
  RnaRun rna;
  rna.addSequence("unique_target", "AUGC");
  rna.addMatch(rna.addObservation("spectrum_1"), "AU[Gm]", {{"unique_target", 0, 2}});
  const IdentificationData& id_data = rna.data;

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Check that target_mapping_count equals 1
  int data_lines_checked = 0;
  for (const auto& line : lines)
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    
    vector<std::string> fields;
    StringUtils::split(line, '\t', fields);
    
    TEST_TRUE(fields.size() >= 12)
    
    if (fields.size() >= 12)
    {
      Int mapping_count = StringUtils::toInt32(fields[11]);
      TEST_EQUAL(mapping_count, 1)
      data_lines_checked++;
    }
  }

  TEST_TRUE(data_lines_checked > 0)
}
END_SECTION

START_SECTION(([EXTRA] decoy-only mappings do not produce modification rows))
{
  // Test that sequences mapping only to decoys have target_mapping_count = 0
  // Only decoy parent sequences, and a match to both
  RnaRun rna;
  rna.addSequence("DECOY_rna_1", "AUGC", true);
  rna.addSequence("DECOY_rna_2", "CAUGC", true);
  rna.addMatch(rna.addObservation("spectrum_1"), "AU[Gm]", {{"DECOY_rna_1", 0, 2}, {"DECOY_rna_2", 1, 3}});
  const IdentificationData& id_data = rna.data;

  // Store to file
  std::string test_file;
  NEW_TMP_FILE(test_file);
  BedRModFile file;
  file.store(test_file, id_data);

  // Read and verify the output
  TextFile output;
  output.load(test_file);
  vector<std::string> lines(output.begin(), output.end());

  // Decoy parents are not genomic modification calls.
  int data_lines_checked = 0;
  for (const auto& line : lines)
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    
    vector<std::string> fields;
    StringUtils::split(line, '\t', fields);
    
    TEST_TRUE(fields.size() >= 12)
    
    if (fields.size() >= 12)
    {
      Int mapping_count = StringUtils::toInt32(fields[11]);
      TEST_EQUAL(mapping_count, 0)
      data_lines_checked++;
    }
  }

  TEST_EQUAL(data_lines_checked, 0)
}
END_SECTION


START_SECTION(([EXTRA] reference metadata is explicit and defaults to unknown))
{
  std::string output;
  NEW_TMP_FILE(output)
  IdentificationData data;
  BedRModFile().store(output, data);
  TextFile unknown(output);
  TEST_EQUAL(*std::next(unknown.begin()), "#organism=NA")
  TEST_TRUE(std::find(unknown.begin(), unknown.end(), "#basecalling=NA") != unknown.end())
  TEST_TRUE(std::find(unknown.begin(), unknown.end(), "#sequencing_platform=NA") != unknown.end())
  BedRModFile().store(output, data, "", {{"organism", "10090"}, {"assembly", "GRCm39"}});
  TextFile explicit_reference(output);
  TEST_EQUAL(*std::next(explicit_reference.begin()), "#organism=10090")
  TEST_TRUE(std::find(explicit_reference.begin(), explicit_reference.end(), "#assembly=GRCm39") != explicit_reference.end())
  TEST_EXCEPTION(Exception::InvalidParameter, BedRModFile().store(output, data, "", {{"organism", "bad\nheader"}}))
  BedRModFile().store(output, data, "", {{"sequencing_platform", "custom-MS"}, {"basecalling", "custom-producer"}});
  TextFile producer(output);
  TEST_TRUE(std::find(producer.begin(), producer.end(), "#sequencing_platform=custom-MS") != producer.end())
  TEST_TRUE(std::find(producer.begin(), producer.end(), "#basecalling=custom-producer") != producer.end())
  TEST_EXCEPTION(Exception::InvalidParameter, BedRModFile().store(output, data, "", {{"basecalling", "bad\nheader"}}))
}
END_SECTION

START_SECTION(([EXTRA] buffered output errors are reported))
{
  // Linux's full device opens successfully but fails on write/flush.
  if (File::exists("/dev/full"))
  {
    TEST_EXCEPTION(Exception::FileNotWritable, BedRModFile().store("/dev/full", IdentificationData()))
  }
  else
  {
    NOT_TESTABLE
  }
}
END_SECTION

START_SECTION(([EXTRA] custom modification names require an honest mapping without truncating output))
{
  RnaRun rna;
  rna.addSequence("rna", "AG");
  rna.addMatch(rna.addObservation("scan"), "[mA?]G", {{"rna", 0, 1}});
  const IdentificationData& data = rna.data;
  std::string output_file, mapping_file;
  NEW_TMP_FILE(output_file)
  NEW_TMP_FILE(mapping_file)
  {
    std::ofstream output(output_file);
    output << "preserve this file\n";
  }
  TEST_EXCEPTION(Exception::InvalidParameter, BedRModFile().store(output_file, data))
  TextFile unchanged(output_file);
  TEST_EQUAL(std::distance(unchanged.begin(), unchanged.end()), 1)
  TEST_EQUAL(*unchanged.begin(), "preserve this file")
  {
    std::ofstream mapping(mapping_file);
    mapping << "mod,modomics_short_name\nmA?,m5C\n";
  }
  TEST_EXCEPTION(Exception::InvalidParameter, BedRModFile().store(output_file, data, mapping_file))
  TextFile still_unchanged(output_file);
  TEST_EQUAL(*still_unchanged.begin(), "preserve this file")
  // An explicit user mapping can identify a compatible concrete modification.
  {
    std::ofstream mapping(mapping_file);
    mapping << "mod,modomics_short_name\nmA?,m6A\n";
  }
  BedRModFile().store(output_file, data, mapping_file);
  TextFile mapped(output_file);
  TEST_TRUE(std::find(mapped.begin(), mapped.end(), "#modification_names=m6A:m6A:A") != mapped.end())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
