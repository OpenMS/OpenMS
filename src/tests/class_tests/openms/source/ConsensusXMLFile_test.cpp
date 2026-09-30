// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/TestFileValidation.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
///////////////////////////

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <fstream>
#include <sstream>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

using namespace OpenMS;
using namespace std;


DRange<1> makeRange(double a, double b)
{
  DPosition<1> pa(a), pb(b);
  return DRange<1>(pa, pb);
}

namespace
{
  // registers a tool-defined modification; ModificationsDB is process-wide, so every section uses its own name
  const ResidueModification* defineMod4b(const std::string& id, char origin, const std::string& formula)
  {
    ResidueModification d;
    d.setId(id);
    d.setOrigin(origin);
    d.setTermSpecificity(ResidueModification::ANYWHERE);
    d.setFullId();
    d.setDiffFormula(EmpiricalFormula(formula));
    d.setDiffMonoMass(EmpiricalFormula(formula).getMonoWeight());
    return ModificationsDB::getInstance()->registerDefinition(d);
  }

  // a definition record for a name that is NOT registered in this process
  std::string freshRecord4b(const std::string& id, char origin, const std::string& formula)
  {
    ResidueModification d;
    d.setId(id);
    d.setOrigin(origin);
    d.setTermSpecificity(ResidueModification::ANYWHERE);
    d.setFullId();
    d.setDiffFormula(EmpiricalFormula(formula));
    d.setDiffMonoMass(EmpiricalFormula(formula).getMonoWeight());
    return d.toDefinitionString();
  }

  std::string slurp4b(const std::string& path)
  {
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }

  bool fileContains4b(const std::string& path, const std::string& needle)
  {
    return slurp4b(path).find(needle) != std::string::npos;
  }

  // first occurrence only; returns false when @p from is absent
  bool replaceInFile4b(const std::string& path, const std::string& from, const std::string& to)
  {
    std::string s = slurp4b(path);
    const std::size_t pos = s.find(from);
    if (pos == std::string::npos) return false;
    s.replace(pos, from.size(), to);
    std::ofstream out(path);
    out << s;
    return true;
  }
}

START_TEST(ConsensusXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusXMLFile * ptr = nullptr;
ConsensusXMLFile* nullPointer = nullptr;
START_SECTION((ConsensusXMLFile()))
ptr = new ConsensusXMLFile();
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~ConsensusXMLFile()))
delete ptr;
END_SECTION

TOLERANCE_ABSOLUTE(0.01)

START_SECTION(const PeakFileOptions& getOptions() const)
ConsensusXMLFile file;
TEST_EQUAL(file.getOptions().hasMSLevels(), false)
END_SECTION

START_SECTION(PeakFileOptions& getOptions())
ConsensusXMLFile file;
file.getOptions().addMSLevel(1);
TEST_EQUAL(file.getOptions().hasMSLevels(), true);
END_SECTION

START_SECTION((void load(const std::string &filename, ConsensusMap & map)))
ConsensusMap map;
ConsensusXMLFile file;
file.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), map);

//test DocumentIdentifier addition
TEST_STRING_EQUAL(map.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"));
TEST_STRING_EQUAL(FileTypes::typeToName(map.getLoadedFileType()), "consensusXML");

//meta data
TEST_EQUAL(map.getIdentifier(), "lsid")
TEST_EQUAL(map.getExperimentType() == "label-free", true)
TEST_EQUAL(map.getMetaValue("name1") == DataValue("value1"), true)
TEST_EQUAL(map.getMetaValue("name2") == DataValue(2), true)
//file descriptions
TEST_EQUAL(map.getColumnHeaders()[0].filename == "data/MapAlignmentFeatureMap1.xml", true)
TEST_EQUAL(map.getColumnHeaders()[0].label, "label")
TEST_EQUAL(map.getColumnHeaders()[0].size, 144)
TEST_EQUAL(map.getColumnHeaders()[0].getMetaValue("name3") == DataValue("value3"), true)
TEST_EQUAL(map.getColumnHeaders()[0].getMetaValue("name4") == DataValue(4), true)
TEST_STRING_EQUAL(map.getColumnHeaders()[1].filename, "data/MapAlignmentFeatureMap2.xml")
TEST_EQUAL(map.getColumnHeaders()[1].label, "")
TEST_EQUAL(map.getColumnHeaders()[1].size, 0)
TEST_EQUAL(map.getColumnHeaders()[1].getMetaValue("name5") == DataValue("value5"), true)
TEST_EQUAL(map.getColumnHeaders()[1].getMetaValue("name6") == DataValue(6.0), true)
//data processing
TEST_EQUAL(map.getDataProcessing().size(), 2)
TEST_STRING_EQUAL(map.getDataProcessing()[0].getSoftware().getName(), "Software1")
TEST_STRING_EQUAL(map.getDataProcessing()[0].getSoftware().getVersion(), "0.91a")
TEST_EQUAL(map.getDataProcessing()[0].getProcessingActions().size(), 1)
TEST_EQUAL(map.getDataProcessing()[0].getProcessingActions().count(DataProcessing::DEISOTOPING), 1)
TEST_STRING_EQUAL(map.getDataProcessing()[0].getMetaValue("name"), "dataProcessing")
TEST_STRING_EQUAL(map.getDataProcessing()[1].getSoftware().getName(), "Software2")
TEST_STRING_EQUAL(map.getDataProcessing()[1].getSoftware().getVersion(), "0.92a")
TEST_EQUAL(map.getDataProcessing()[1].getProcessingActions().size(), 2)
TEST_EQUAL(map.getDataProcessing()[1].getProcessingActions().count(DataProcessing::SMOOTHING), 1)
TEST_EQUAL(map.getDataProcessing()[1].getProcessingActions().count(DataProcessing::BASELINE_REDUCTION), 1)
//protein identifications
TEST_EQUAL(map.getProteinIdentifications().size(), 2)
TEST_EQUAL(map.getProteinIdentifications()[0].getHits().size(), 2)
TEST_EQUAL(map.getProteinIdentifications()[0].getHits()[0].getSequence(), "ABCDEFG")
TEST_EQUAL(map.getProteinIdentifications()[0].getHits()[1].getSequence(), "HIJKLMN")
TEST_EQUAL(map.getProteinIdentifications()[1].getHits().size(), 1)
TEST_EQUAL(map.getProteinIdentifications()[1].getHits()[0].getSequence(), "OPQREST")
//peptide identifications
TEST_EQUAL(map[0].getPeptideIdentifications().size(), 2)
TEST_EQUAL(map[0].getPeptideIdentifications()[0].getHits().size(), 1)
TEST_EQUAL(map[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("A"))
TEST_EQUAL(map[0].getPeptideIdentifications()[1].getHits().size(), 2)
TEST_EQUAL(map[0].getPeptideIdentifications()[1].getHits()[0].getSequence(), AASequence::fromString("C"))
TEST_EQUAL(map[0].getPeptideIdentifications()[1].getHits()[1].getSequence(), AASequence::fromString("D"))
TEST_EQUAL(map[1].getPeptideIdentifications().size(), 1)
TEST_EQUAL(map[1].getPeptideIdentifications()[0].getHits().size(), 1)
TEST_EQUAL(map[1].getPeptideIdentifications()[0].getHits()[0].getSequence(),AASequence::fromString( "E"))
//unassigned peptide identifications
TEST_EQUAL(map.getUnassignedPeptideIdentifications().size(), 2)
TEST_EQUAL(map.getUnassignedPeptideIdentifications()[0].getHits().size(), 1)
TEST_EQUAL(map.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("F"))
TEST_EQUAL(map.getUnassignedPeptideIdentifications()[1].getHits().size(), 2)
TEST_EQUAL(map.getUnassignedPeptideIdentifications()[1].getHits()[0].getSequence(), AASequence::fromString("G"))
TEST_EQUAL(map.getUnassignedPeptideIdentifications()[1].getHits()[1].getSequence(), AASequence::fromString("H"))

//features
TEST_EQUAL(map.size(), 6)
ConsensusFeature cons_feature = map[0];
TEST_REAL_SIMILAR(cons_feature.getRT(), 1273.27)
TEST_REAL_SIMILAR(cons_feature.getMZ(), 904.47)
TEST_REAL_SIMILAR(cons_feature.getIntensity(), 3.12539e+07)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().minPosition()[0], 1273.27)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().maxPosition()[0], 1273.27)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().minPosition()[1], 904.47)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().maxPosition()[1], 904.47)
TEST_REAL_SIMILAR(cons_feature.getIntensityRange().minPosition()[0], 3.12539e+07)
TEST_REAL_SIMILAR(cons_feature.getIntensityRange().maxPosition()[0], 3.12539e+07)
TEST_REAL_SIMILAR(cons_feature.getQuality(), 1.1)
TEST_EQUAL(cons_feature.getMetaValue("peptide_id") == DataValue("RefSeq:NC_1234"), true)
ConsensusFeature::HandleSetType::const_iterator it = cons_feature.begin();
TEST_REAL_SIMILAR(it->getIntensity(), 3.12539e+07)

cons_feature = map[5];
TEST_REAL_SIMILAR(cons_feature.getRT(), 1194.82)
TEST_REAL_SIMILAR(cons_feature.getMZ(), 777.101)
TEST_REAL_SIMILAR(cons_feature.getIntensity(), 1.78215e+07)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().minPosition()[0], 1194.82)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().maxPosition()[0], 1194.82)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().minPosition()[1], 777.101)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().maxPosition()[1], 777.101)
TEST_REAL_SIMILAR(cons_feature.getIntensityRange().minPosition()[0], 1.78215e+07)
TEST_REAL_SIMILAR(cons_feature.getIntensityRange().maxPosition()[0], 1.78215e+07)
TEST_REAL_SIMILAR(cons_feature.getQuality(), 0.0)
it = cons_feature.begin();
TEST_REAL_SIMILAR(it->getIntensity(), 1.78215e+07)
++ it;
TEST_REAL_SIMILAR(it->getIntensity(), 1.78215e+07)

// test meta values:
TEST_EQUAL(map[0].getMetaValue("myIntList") == ListUtils::create<Int>("1,10,12"), true);
TEST_EQUAL(map[0].getMetaValue("myDoubleList") == ListUtils::create<double>("1.111,10.999,12.45"), true);
std::cout << "list: " << map[0].getMetaValue("myStringList") << "\n";
TEST_EQUAL(map[0].getMetaValue("myStringList") == ListUtils::create<std::string>("myABC1,Stuff,12"), true);
TEST_EQUAL(map[4].getMetaValue("myDoubleList") == ListUtils::create<double>("6.442"), true);

//PeakFileOptions tests

file.getOptions().setRTRange(makeRange(815, 818));
file.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_2_options.consensusXML"), map);
TEST_EQUAL(map.size(), 1)
TEST_REAL_SIMILAR(map[0].getRT(), 817.266)

file.getOptions() = PeakFileOptions();
file.getOptions().setMZRange(makeRange(944, 945));
file.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_2_options.consensusXML"), map);
TEST_EQUAL(map.size(), 1)
TEST_REAL_SIMILAR(map[0].getMZ(), 944.96)

file.getOptions() = PeakFileOptions();
file.getOptions().setIntensityRange(makeRange(15000, 24000));
file.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_2_options.consensusXML"), map);
TEST_EQUAL(map.size(), 1)
TEST_REAL_SIMILAR(map[0].getIntensity(), 23000.238)

END_SECTION

START_SECTION((void store(const std::string &filename, const ConsensusMap &consensus_map)))
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  ConsensusMap map;
  ConsensusXMLFile f;

  f.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), map);

  // make protIDs non-unique
  map.getProteinIdentifications().push_back(map.getProteinIdentifications()[0]);
  TEST_EXCEPTION(Exception::ParseError, f.store(tmp_filename, map))
  map.getProteinIdentifications().pop_back(); // undo

  f.store(tmp_filename, map);
  WHITELIST("?xml-stylesheet")
  TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), tmp_filename)


END_SECTION

START_SECTION([EXTRA](bool isValid(const std::string &filename)))
  ConsensusXMLFile f;
  TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), std::cerr), true);
  TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_2_options.consensusXML"), std::cerr), true);

  //test if written empty file
  // - this is invalid, so it is not tested :)

  //test if written full file is valid
  ConsensusMap m;
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  f.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), m);
  f.store(tmp_filename, m);
  TEST_EQUAL(f.isValid(tmp_filename, std::cerr), true);
END_SECTION

START_SECTION([EXTRA] Compressed file writing - gzip round-trip)
  ConsensusXMLFile f;
  ConsensusMap map;
  f.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), map);

  // Store as gzip-compressed file
  std::string tmp_gz;
  NEW_TMP_FILE(tmp_gz);
  tmp_gz += ".gz";
  f.store(tmp_gz, map);

  // Load back from compressed file
  ConsensusMap map_gz;
  f.load(tmp_gz, map_gz);

  // Verify round-trip integrity
  TEST_EQUAL(map_gz.size(), map.size())
  TEST_EQUAL(map_gz.getColumnHeaders().size(), map.getColumnHeaders().size())
END_SECTION

START_SECTION([EXTRA] Compressed file writing - bzip2 round-trip)
  ConsensusXMLFile f;
  ConsensusMap map;
  f.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), map);

  // Store as bzip2-compressed file
  std::string tmp_bz2;
  NEW_TMP_FILE(tmp_bz2);
  tmp_bz2 += ".bz2";
  f.store(tmp_bz2, map);

  // Load back from compressed file
  ConsensusMap map_bz2;
  f.load(tmp_bz2, map_bz2);

  // Verify round-trip integrity
  TEST_EQUAL(map_bz2.size(), map.size())
  TEST_EQUAL(map_bz2.getColumnHeaders().size(), map.getColumnHeaders().size())
END_SECTION

START_SECTION([EXTRA] Protein group quantities round-trip)
{
  ConsensusXMLFile f;
  ConsensusMap map;
  f.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), map);

  // attach quantities the way PeptideAndProteinQuant::annotateQuantificationsToProteins() does
  ProteinIdentification& prot = map.getProteinIdentifications()[0];
  ProteinIdentification::ProteinGroup group;
  for (const ProteinHit& hit : prot.getHits()) { group.accessions.push_back(hit.getAccession()); }
  TEST_NOT_EQUAL(group.accessions.size(), 0)
  group.probability = 0.75;

  group.getFloatDataArrays().resize(4);
  group.getStringDataArrays().resize(2);
  group.getIntegerDataArrays().resize(4);
  group.getFloatDataArrays()[0].setName("psm_count");
  group.getFloatDataArrays()[1].setName("distinct_peptides");
  group.getFloatDataArrays()[2].setName("file_channel_level_abundance");
  group.getFloatDataArrays()[2].assign({10.0f, 20.0f, 30.0f, 40.0f});
  group.getFloatDataArrays()[3].setName("fraction_group_level_abundance");
  group.getFloatDataArrays()[3].assign({1.5f, 2.5f});
  group.getStringDataArrays()[0].setName("file_channel_level_filename");
  group.getStringDataArrays()[0].assign({"fileA", "fileA", "fileB", "fileB"});
  group.getStringDataArrays()[1].setName("file_level_filename");
  group.getIntegerDataArrays()[0].setName("file_channel_level_channel");
  group.getIntegerDataArrays()[0].assign({1, 2, 1, 2});
  group.getIntegerDataArrays()[1].setName("file_level_psm_count");
  group.getIntegerDataArrays()[2].setName("fraction_group_level_fraction_group");
  group.getIntegerDataArrays()[2].assign({1, 2});
  group.getIntegerDataArrays()[3].setName("fraction_group_level_label");
  group.getIntegerDataArrays()[3].assign({1, 1});

  prot.insertIndistinguishableProteins(group);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  f.store(tmp_filename, map);

  ConsensusMap loaded;
  f.load(tmp_filename, loaded);

  TEST_EQUAL(loaded.getProteinIdentifications()[0].getIndistinguishableProteins().size(), 1)
  const ProteinIdentification::ProteinGroup& rt = loaded.getProteinIdentifications()[0].getIndistinguishableProteins()[0];

  // The assay-only layout is restored without inventing a legacy sample array.
  TEST_EQUAL(rt.getFloatDataArrays().size(), 4)
  TEST_EQUAL(rt.getStringDataArrays().size(), 2)
  TEST_EQUAL(rt.getIntegerDataArrays().size(), 4)
  TEST_EQUAL(rt.getFloatDataArrays()[2].getName(), "file_channel_level_abundance")
  TEST_EQUAL(rt.getFloatDataArrays()[3].getName(), "fraction_group_level_abundance")
  TEST_EQUAL(rt.getStringDataArrays()[0].getName(), "file_channel_level_filename")
  TEST_EQUAL(rt.getIntegerDataArrays()[0].getName(), "file_channel_level_channel")

  TEST_EQUAL(rt.getFloatDataArrays()[2].size(), 4)
  TEST_REAL_SIMILAR(rt.getFloatDataArrays()[2][0], 10.0)
  TEST_REAL_SIMILAR(rt.getFloatDataArrays()[2][3], 40.0)
  TEST_EQUAL(rt.getFloatDataArrays()[3].size(), 2)
  TEST_REAL_SIMILAR(rt.getFloatDataArrays()[3][0], 1.5)
  TEST_REAL_SIMILAR(rt.getFloatDataArrays()[3][1], 2.5)
  TEST_EQUAL(rt.getStringDataArrays()[0].size(), 4)
  TEST_EQUAL(rt.getStringDataArrays()[0][0], "fileA")
  TEST_EQUAL(rt.getStringDataArrays()[0][3], "fileB")
  TEST_EQUAL(rt.getIntegerDataArrays()[0].size(), 4)
  TEST_EQUAL(rt.getIntegerDataArrays()[0][0], 1)
  TEST_EQUAL(rt.getIntegerDataArrays()[0][3], 2)

  // All-zero count arrays are not written and no longer borrow a length from sample abundances.
  TEST_TRUE(rt.getFloatDataArrays()[0].empty())
  TEST_TRUE(rt.getFloatDataArrays()[1].empty())

  // no leftovers of the encoding on the ProteinIdentification
  std::vector<std::string> keys;
  loaded.getProteinIdentifications()[0].getKeys(keys);
  for (const std::string& key : keys)
  {
    TEST_FALSE(StringUtils::hasSubstring(key, "_abundances"))
    TEST_FALSE(StringUtils::hasSubstring(key, "_quantified_proteins"))
  }

  // and the file is still schema-valid
  TEST_TRUE(f.isValid(tmp_filename, std::cerr))

  // storing again must be stable
  std::string tmp_filename2;
  NEW_TMP_FILE(tmp_filename2);
  f.store(tmp_filename2, loaded);
  ConsensusMap loaded2;
  f.load(tmp_filename2, loaded2);
  TEST_TRUE(loaded2.getProteinIdentifications()[0].getIndistinguishableProteins()[0] == rt)
}
END_SECTION

START_SECTION([EXTRA] Quantities whose owner no longer matches the group are discarded)
{
  // Each quantified group records which proteins its quantities were computed for, so that a tool
  // which filters or renumbers groups without knowing about the quantity params cannot cause them to
  // be silently adopted by a different protein.
  //
  // The recorded value must be ACCESSIONS. The obvious alternative - the "PH_<n>" references the group
  // itself is stored with - does not work: those are positional indices into the protein hit list,
  // reassigned on every store, so dropping a leading hit shifts the PH_ numbers and the group indices
  // by the same amount and a stale entry still compares equal.
  ConsensusXMLFile f;
  ConsensusMap map;
  f.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), map);

  ProteinIdentification& prot = map.getProteinIdentifications()[0];
  TEST_EQUAL(prot.getHits().size(), 2)
  const std::string quantified_accession = prot.getHits()[0].getAccession();

  ProteinIdentification::ProteinGroup g;
  g.accessions.push_back(quantified_accession);
  g.probability = 0.9;
  g.getFloatDataArrays().resize(1);
  g.getFloatDataArrays()[0].setName("abundances");
  g.getFloatDataArrays()[0].assign({11.0f, 22.0f});
  prot.insertIndistinguishableProteins(g);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  f.store(tmp_filename, map);

  // the ownership record holds the accession, not a positional PH_ reference
  TextFile stored;
  stored.load(tmp_filename);
  std::string guard_line;
  for (const std::string& line : stored)
  {
    if (StringUtils::hasSubstring(line, "_quantified_proteins")) { guard_line = line; }
  }
  TEST_EQUAL(guard_line.empty(), false)
  TEST_EQUAL(StringUtils::hasSubstring(guard_line, quantified_accession), true)
  TEST_EQUAL(StringUtils::hasSubstring(guard_line, "PH_"), false)

  // rewrite it to name a different protein, as a renumbering would effectively do
  TextFile rewritten;
  for (const std::string& line : stored)
  {
    std::string out = line;
    if (StringUtils::hasSubstring(out, "_quantified_proteins"))
    {
      StringUtils::substitute(out, quantified_accession, "SOME_OTHER_PROTEIN");
    }
    rewritten.addLine(out);
  }
  std::string tampered;
  NEW_TMP_FILE(tampered);
  rewritten.store(tampered);

  ConsensusMap loaded;
  f.load(tampered, loaded);
  const auto& groups = loaded.getProteinIdentifications()[0].getIndistinguishableProteins();
  TEST_EQUAL(groups.size(), 1)
  TEST_EQUAL(groups[0].accessions[0], quantified_accession)

  // the quantities belong to a different protein now, so they must be dropped, not handed over
  bool has_abundances = false;
  for (const auto& fda : groups[0].getFloatDataArrays())
  {
    if (fda.getName() == "abundances" && !fda.empty()) { has_abundances = true; }
  }
  TEST_EQUAL(has_abundances, false)

  // and nothing of the discarded block may linger as a stray UserParam
  std::vector<std::string> keys;
  loaded.getProteinIdentifications()[0].getKeys(keys);
  for (const std::string& key : keys)
  {
    TEST_EQUAL(StringUtils::hasSubstring(key, "indistinguishable_proteins_0_"), false)
  }
}
END_SECTION

START_SECTION([EXTRA] Legacy all-zero sample abundances keep their length across a round-trip)
{
  // Preserve the shape of old consensusXML annotations even though sample abundances no longer
  // identify a quantified group. This compatibility path must not invent assay arrays.
  ConsensusXMLFile f;
  ConsensusMap map;
  f.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), map);

  ProteinIdentification& prot = map.getProteinIdentifications()[0];
  ProteinIdentification::ProteinGroup group;
  for (const ProteinHit& hit : prot.getHits()) { group.accessions.push_back(hit.getAccession()); }
  group.probability = 0.5;
  group.getFloatDataArrays().resize(4);
  group.getStringDataArrays().resize(2);
  group.getIntegerDataArrays().resize(2);
  group.getFloatDataArrays()[0].setName("abundances");
  group.getFloatDataArrays()[0].assign(3, 0.0f); // quantified everywhere as zero
  group.getFloatDataArrays()[1].setName("psm_count");
  group.getFloatDataArrays()[1].resize(3);
  group.getFloatDataArrays()[2].setName("distinct_peptides");
  group.getFloatDataArrays()[2].resize(3);
  group.getFloatDataArrays()[3].setName("file_channel_level_abundance");
  group.getFloatDataArrays()[3].assign(3, 0.0f);
  group.getStringDataArrays()[0].setName("file_channel_level_filename");
  group.getStringDataArrays()[0].assign({"fA", "fB", "fC"});
  group.getStringDataArrays()[1].setName("file_level_filename");
  group.getIntegerDataArrays()[0].setName("file_channel_level_channel");
  group.getIntegerDataArrays()[0].assign({1, 1, 1});
  group.getIntegerDataArrays()[1].setName("file_level_psm_count");
  prot.insertIndistinguishableProteins(group);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  f.store(tmp_filename, map);
  ConsensusMap loaded;
  f.load(tmp_filename, loaded);

  const ProteinIdentification::ProteinGroup& rt = loaded.getProteinIdentifications()[0].getIndistinguishableProteins().back();
  TEST_EQUAL(rt.getFloatDataArrays()[0].getName(), "abundances")
  TEST_EQUAL(rt.getFloatDataArrays()[0].size(), 3)
  TEST_EQUAL(rt.getFloatDataArrays()[1].size(), 3)
  TEST_EQUAL(rt.getFloatDataArrays()[2].size(), 3)
  TEST_EQUAL(rt.getFloatDataArrays()[3].size(), 3)
  TEST_REAL_SIMILAR(rt.getFloatDataArrays()[0][0], 0.0)
}
END_SECTION

START_SECTION([EXTRA] Protein groups without quantities are written unchanged)
{
  // Groups that carry no data arrays - everything outside ProteomicsLFQ/IsobaricWorkflow - must not
  // gain any output, otherwise every reference file with protein groups would change.
  ConsensusXMLFile f;
  ConsensusMap map;
  f.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), map);

  ProteinIdentification& prot = map.getProteinIdentifications()[0];
  ProteinIdentification::ProteinGroup group;
  for (const ProteinHit& hit : prot.getHits()) { group.accessions.push_back(hit.getAccession()); }
  group.probability = 0.5;
  prot.insertIndistinguishableProteins(group);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  f.store(tmp_filename, map);

  // nothing beyond the plain "indistinguishable_proteins_0" entry may be emitted for this group
  // (the fixture carries unrelated list-typed UserParams, so only the group prefix can be checked)
  TextFile stored;
  stored.load(tmp_filename);
  Size quant_lines = 0;
  for (const std::string& line : stored)
  {
    if (StringUtils::hasSubstring(line, "indistinguishable_proteins_0_")) { ++quant_lines; }
  }
  TEST_EQUAL(quant_lines, 0)

  ConsensusMap loaded;
  f.load(tmp_filename, loaded);
  TEST_EQUAL(loaded.getProteinIdentifications()[0].getIndistinguishableProteins().size(), 1)
  TEST_EQUAL(loaded.getProteinIdentifications()[0].getIndistinguishableProteins()[0].getFloatDataArrays().empty(), true)
}
END_SECTION

START_SECTION([EXTRA] store/load - tool-defined modifications on assigned and unassigned identifications travel with their definitions)
{
  TEST_TRUE(defineMod4b("TestCXML:Assigned", 'K', "C2H2O") != nullptr)
  TEST_TRUE(defineMod4b("TestCXML:Unassigned", 'R', "CH2") != nullptr)
  ConsensusMap map;
  map.ensureUniqueId();
  map.getColumnHeaders()[0].filename = "file0.mzML";
  map.getColumnHeaders()[0].size = 1;
  ProteinIdentification prot;
  prot.setIdentifier("run4b");
  prot.setDateTime(DateTime::now());
  map.getProteinIdentifications().push_back(prot);

  ConsensusFeature f;
  f.setRT(100.0);
  f.setMZ(500.0);
  f.setIntensity(1000.0);
  f.ensureUniqueId();
  f.insert(FeatureHandle(0, f));
  PeptideIdentification pa;
  pa.setIdentifier("run4b");
  PeptideHit ha;
  ha.setSequence(AASequence::fromString("PEPK(TestCXML:Assigned)IDE"));
  pa.insertHit(ha);
  f.getPeptideIdentifications().push_back(pa);
  map.push_back(f);

  PeptideIdentification pu;
  pu.setIdentifier("run4b");
  PeptideHit hu;
  hu.setSequence(AASequence::fromString("PEPR(TestCXML:Unassigned)IDE"));
  pu.insertHit(hu);
  map.getUnassignedPeptideIdentifications().push_back(pu);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename)
  ConsensusXMLFile().store(tmp_filename, map);
  TEST_TRUE(fileContains4b(tmp_filename, "name=\"modification_definitions\""))
  TEST_TRUE(fileContains4b(tmp_filename, "1|TestCXML:Assigned|TestCXML:Assigned (K)|"))
  TEST_TRUE(fileContains4b(tmp_filename, "1|TestCXML:Unassigned|TestCXML:Unassigned (R)|"))

  ConsensusMap in;
  ConsensusXMLFile().load(tmp_filename, in);
  TEST_EQUAL(in.size(), 1)
  TEST_EQUAL(in.getUnassignedPeptideIdentifications().size(), 1)
  if (in.size() == 1 && !in[0].getPeptideIdentifications().empty() && !in[0].getPeptideIdentifications()[0].getHits().empty())
  {
    TEST_EQUAL(in[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPK(TestCXML:Assigned)IDE")
  }
  if (in.getUnassignedPeptideIdentifications().size() == 1 && !in.getUnassignedPeptideIdentifications()[0].getHits().empty())
  {
    TEST_EQUAL(in.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPR(TestCXML:Unassigned)IDE")
  }
}
END_SECTION

START_SECTION([EXTRA] load - definitions are registered before the sequences are parsed)
{
  const ModificationsDB* db = ModificationsDB::getInstance();
  TEST_FALSE(db->hasDefinedModification("TestCXML:Fresh"))
  ConsensusMap map;
  map.ensureUniqueId();
  map.getColumnHeaders()[0].filename = "file0.mzML";
  map.getColumnHeaders()[0].size = 1;
  ProteinIdentification prot;
  prot.setIdentifier("run4b_fresh");
  prot.setDateTime(DateTime::now());
  ProteinIdentification::SearchParameters sp;
  sp.setMetaValue(Constants::UserParam::MODIFICATION_DEFINITIONS, freshRecord4b("TestCXML:Fresh", 'K', "C2H2O"));
  prot.setSearchParameters(sp);
  map.getProteinIdentifications().push_back(prot);
  ConsensusFeature f;
  f.setRT(100.0);
  f.setMZ(500.0);
  f.setIntensity(1000.0);
  f.ensureUniqueId();
  f.insert(FeatureHandle(0, f));
  PeptideIdentification pep;
  pep.setIdentifier("run4b_fresh");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTKIDE"));
  pep.insertHit(hit);
  f.getPeptideIdentifications().push_back(pep);
  map.push_back(f);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename)
  ConsensusXMLFile().store(tmp_filename, map);
  TEST_FALSE(db->hasDefinedModification("TestCXML:Fresh")) // storing registers nothing
  // a hit using the not-yet-registered name, as a file from another process would carry it
  TEST_TRUE(replaceInFile4b(tmp_filename, "sequence=\"PEPTKIDE\"", "sequence=\"PEPTK(TestCXML:Fresh)IDE\""))

  ConsensusMap in;
  ConsensusXMLFile().load(tmp_filename, in);
  TEST_TRUE(db->hasDefinedModification("TestCXML:Fresh"))
  if (in.size() == 1 && !in[0].getPeptideIdentifications().empty() && !in[0].getPeptideIdentifications()[0].getHits().empty())
  {
    TEST_EQUAL(in[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPTK(TestCXML:Fresh)IDE")
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// check the temporary files written above against their XML schema (types without a validator are skipped)
VALIDATE_TMP_FILES

END_TEST
