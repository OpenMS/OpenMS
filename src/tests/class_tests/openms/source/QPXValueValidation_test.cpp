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
#include <OpenMS/FORMAT/QPXValueValidation.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <arrow/api.h>

using namespace OpenMS;
using namespace std;

namespace
{
  /// Swap one column of @p table for @p values, keeping the schema. The way every section below
  /// turns a valid table into the one invalid shape it is about.
  std::shared_ptr<arrow::Table> replaceColumn(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name,
    const std::shared_ptr<arrow::Array>& values)
  {
    std::vector<std::shared_ptr<arrow::ChunkedArray>> columns = table->columns();
    columns.at(static_cast<size_t>(table->schema()->GetFieldIndex(name))) =
      std::make_shared<arrow::ChunkedArray>(values);
    return arrow::Table::Make(table->schema(), std::move(columns));
  }

  std::shared_ptr<arrow::Array> stringArray(const std::vector<std::string>& values)
  {
    arrow::StringBuilder b;
    for (const auto& v : values) { (void)b.Append(v); }
    return b.Finish().ValueOrDie();
  }

  std::shared_ptr<arrow::Array> floatArray(const std::vector<float>& values)
  {
    arrow::FloatBuilder b;
    for (float v : values) { (void)b.Append(v); }
    return b.Finish().ValueOrDie();
  }

  /// One identification to put in a PSM table. An empty @p spectrum_reference is how a real
  /// producer (protein inference, an ID converter, a transferred ID) yields an empty scan list.
  struct PSMSpec
  {
    std::string sequence;
    Int charge;
    std::string spectrum_reference;
  };

  /// Build the psm view from @p specs, all belonging to one run whose path is set.
  std::shared_ptr<arrow::Table> makePSMTableFrom(const std::vector<PSMSpec>& specs)
  {
    ProteinIdentification prot;
    prot.setIdentifier("run1");
    prot.setScoreType("q-value");
    prot.setPrimaryMSRunPath({"/data/run1.mzML"});
    ProteinHit ph;
    ph.setAccession("PROT_A");
    prot.setHits({ph});

    PeptideIdentificationList peps;
    double rt = 100.0;
    for (const auto& spec : specs)
    {
      PeptideIdentification pep;
      pep.setIdentifier("run1");
      pep.setRT(rt);
      pep.setMZ(500.0 + rt);
      rt += 1.0;
      pep.setScoreType("q-value");
      if (!spec.spectrum_reference.empty()) { pep.setSpectrumReference(spec.spectrum_reference); }
      PeptideHit hit;
      hit.setSequence(AASequence::fromString(spec.sequence));
      hit.setCharge(spec.charge);
      hit.setScore(0.01);
      PeptideEvidence ev;
      ev.setProteinAccession("PROT_A");
      hit.setPeptideEvidences({ev});
      pep.setHits({hit});
      peps.push_back(pep);
    }
    return QPXFile::exportPSMsToQPXArrow({prot}, peps);
  }

  /// A minimal but complete PSM table: @p n identifications of one run, distinct scans.
  std::shared_ptr<arrow::Table> makePSMTable(size_t n = 2)
  {
    std::vector<PSMSpec> specs;
    for (size_t i = 0; i < n; ++i)
    {
      specs.push_back({i % 2 == 0 ? "PEPTIDER" : "DFPIANGER", 2,
                       "controllerType=0 controllerNumber=1 scan=" + StringUtils::toStr(1000 + i)});
    }
    return makePSMTableFrom(specs);
  }

  /// A label-free consensus map with @p n unmapped features in one run, all charge 2.
  ConsensusMap makeUnmappedFeatureMap(const std::vector<std::pair<double, double>>& rt_mz)
  {
    ConsensusMap cmap;
    cmap.setExperimentType("label-free");
    ConsensusMap::ColumnHeader h;
    h.filename = "run1.mzML";
    h.label = "";
    cmap.getColumnHeaders()[0] = h;

    ProteinIdentification prot;
    prot.setIdentifier("run1");
    prot.setScoreType("q-value");
    prot.setPrimaryMSRunPath({"/data/run1.mzML"});
    cmap.setProteinIdentifications({prot});

    for (const auto& [rt, mz] : rt_mz)
    {
      ConsensusFeature cf;
      cf.setRT(rt);
      cf.setMZ(mz);
      cf.setCharge(2);
      cf.setIntensity(1000.0f);
      BaseFeature bf;
      bf.setRT(rt);
      bf.setMZ(mz);
      bf.setCharge(2);
      bf.setIntensity(1000.0f);
      cf.insert(0, bf);
      cmap.push_back(cf); // no identification: an unmapped feature
    }
    return cmap;
  }
} // namespace

START_TEST(QPXValueValidation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

QPXValueValidation* ptr = nullptr;
QPXValueValidation* null_ptr = nullptr;

START_SECTION((QPXValueValidation(View view)))
{
  ptr = new QPXValueValidation(QPXValueValidation::View::PSM);
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((~QPXValueValidation()))
{
  delete ptr;
}
END_SECTION

START_SECTION(([QPXValueValidation::Result] std::string toString() const))
{
  QPXValueValidation::Result result;
  TEST_TRUE(result.valid)
  TEST_EQUAL(result.toString(), "")

  result.valid = false;
  result.errors.push_back("first");
  TEST_EQUAL(result.toString(), "first")
  result.errors.push_back("second");
  // One log line per refusal, so the diagnostics have to be joined rather than dropped.
  TEST_TRUE(result.toString().find("first") != std::string::npos)
  TEST_TRUE(result.toString().find("second") != std::string::npos)
}
END_SECTION

START_SECTION((Result validate(const std::shared_ptr<arrow::Table>& table)))
{
  auto table = makePSMTable();
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 2)

  QPXValueValidation validator(QPXValueValidation::View::PSM);
  const auto result = validator.validate(table);
  TEST_TRUE(result.valid)
  TEST_EQUAL(result.errors.size(), 0)

  // A null table is a caller error, not a valid empty export.
  QPXValueValidation null_validator(QPXValueValidation::View::PSM);
  TEST_FALSE(null_validator.validate(nullptr).valid)

  // A 0-row table is legal: the streaming writers emit one for empty input.
  QPXValueValidation empty_validator(QPXValueValidation::View::PSM);
  TEST_TRUE(empty_validator.validate(makePSMTable(0)).valid)
}
END_SECTION

START_SECTION(([EXTRA] psm view: identity-bearing key components must carry a value))
{
  auto table = makePSMTable();

  // scan: derived from the spectrum reference. Protein-inference output, some ID converters and
  // transferred IDs carry none, and QPXFile deliberately emits the row with an empty scan list -
  // so build that shape rather than patching the column, which is a list<int32>, not a string.
  {
    auto no_scan = makePSMTableFrom({{"PEPTIDER", 2, ""},
                                     {"DFPIANGER", 2, "controllerType=0 controllerNumber=1 scan=1001"}});
    QPXValueValidation validator(QPXValueValidation::View::PSM);
    const auto result = validator.validate(no_scan);
    TEST_FALSE(result.valid)
    TEST_TRUE(result.toString().find("scan") != std::string::npos)
  }

  // run_file_name: the key the psm view joins to the feature and pg views on.
  {
    auto no_run = replaceColumn(table, QPXPSMSchema::RUN_FILE_NAME, stringArray({"", "run1"}));
    QPXValueValidation validator(QPXValueValidation::View::PSM);
    const auto result = validator.validate(no_run);
    TEST_FALSE(result.valid)
    TEST_TRUE(result.toString().find("run_file_name") != std::string::npos)
  }

  // Blank is not a value either - it would join to nothing just as an empty string does.
  {
    auto blank_run = replaceColumn(table, QPXPSMSchema::RUN_FILE_NAME, stringArray({"   ", "run1"}));
    QPXValueValidation validator(QPXValueValidation::View::PSM);
    TEST_FALSE(validator.validate(blank_run).valid)
  }
}
END_SECTION

START_SECTION(([EXTRA] psm view: the primary key must be unique))
{
  // The shape IDMapper produces: one MS2 identification copied onto two consensus features, so
  // the same (peptidoform, charge, run, scan) is exported twice.
  const std::string ref = "controllerType=0 controllerNumber=1 scan=1000";
  auto duplicate = makePSMTableFrom({{"PEPTIDER", 2, ref}, {"PEPTIDER", 2, ref}});
  TEST_EQUAL(duplicate->num_rows(), 2)

  QPXValueValidation validator(QPXValueValidation::View::PSM);
  const auto result = validator.validate(duplicate);
  TEST_FALSE(result.valid)
  TEST_TRUE(result.toString().find("primary key") != std::string::npos)
}
END_SECTION

START_SECTION(([EXTRA] feature view: unmapped features are separated by observed_mz))
{
  // Two unmapped features of the same charge in one run, at retention times that differ but
  // round to the same float32 (one ULP is 244 us at 3000 s). Without observed_mz in the key
  // they would collide and the whole export would be refused - and unmapped features are the
  // common case for ProteomicsLFQ with seeds, not an exotic one.
  TEST_EQUAL(static_cast<float>(3000.0), static_cast<float>(3000.0001))
  auto cmap = makeUnmappedFeatureMap({{3000.0, 400.0}, {3000.0001, 500.0}});
  auto table = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 2)

  QPXValueValidation validator(QPXValueValidation::View::FEATURE);
  TEST_TRUE(validator.validate(table).valid)

  // Agreeing on the m/z as well makes the two rows genuinely indistinguishable: still refused.
  auto same_mz = replaceColumn(table, QPXFeatureSchema::OBSERVED_MZ, floatArray({400.0f, 400.0f}));
  same_mz = replaceColumn(same_mz, QPXFeatureSchema::RT, floatArray({3000.0f, 3000.0f}));
  QPXValueValidation duplicate_validator(QPXValueValidation::View::FEATURE);
  const auto duplicate_result = duplicate_validator.validate(same_mz);
  TEST_FALSE(duplicate_result.valid)
  TEST_TRUE(duplicate_result.toString().find("primary key") != std::string::npos)
}
END_SECTION

START_SECTION(([EXTRA] feature view: non-finite key coordinates are refused))
{
  auto cmap = makeUnmappedFeatureMap({{3000.0, 400.0}});
  auto table = ConsensusMapArrowExport::exportToArrow(cmap);
  TEST_NOT_EQUAL(table, nullptr)

  {
    auto nan_rt = replaceColumn(table, QPXFeatureSchema::RT,
                                floatArray({std::numeric_limits<float>::quiet_NaN()}));
    QPXValueValidation validator(QPXValueValidation::View::FEATURE);
    const auto result = validator.validate(nan_rt);
    TEST_FALSE(result.valid)
    TEST_TRUE(result.toString().find("rt") != std::string::npos)
  }
  {
    auto inf_mz = replaceColumn(table, QPXFeatureSchema::OBSERVED_MZ,
                                floatArray({std::numeric_limits<float>::infinity()}));
    QPXValueValidation validator(QPXValueValidation::View::FEATURE);
    const auto result = validator.validate(inf_mz);
    TEST_FALSE(result.valid)
    TEST_TRUE(result.toString().find("observed_mz") != std::string::npos)
  }
}
END_SECTION

START_SECTION(([EXTRA] a validator spans batches, so duplicates across a batch boundary are found))
{
  // This is why the streaming writers keep one instance for the whole file rather than one per
  // batch: a key repeated in a later batch is still a duplicate in the finished Parquet file.
  auto table = makePSMTable();
  QPXValueValidation validator(QPXValueValidation::View::PSM);
  TEST_TRUE(validator.validate(table).valid)
  const auto second = validator.validate(table); // same rows again
  TEST_FALSE(second.valid)
  TEST_TRUE(second.toString().find("primary key") != std::string::npos)
}
END_SECTION

START_SECTION((void reset()))
{
  auto table = makePSMTable();
  QPXValueValidation validator(QPXValueValidation::View::PSM);
  TEST_TRUE(validator.validate(table).valid)
  TEST_FALSE(validator.validate(table).valid)

  // reset() starts a different output file, so the same keys are legal again.
  validator.reset();
  TEST_TRUE(validator.validate(table).valid)
}
END_SECTION

START_SECTION(([EXTRA] a rejected batch does not advance the accumulated key state))
{
  // Documented behaviour: state advances only when the whole table is valid, so a caller can
  // fix a refused batch and submit it again without first discarding the keys of the batches
  // that already passed.
  auto good = makePSMTable();
  auto bad = replaceColumn(good, QPXPSMSchema::RUN_FILE_NAME, stringArray({"", "run1"}));

  QPXValueValidation validator(QPXValueValidation::View::PSM);
  TEST_FALSE(validator.validate(bad).valid)
  // The valid row of the refused batch must not have been remembered.
  TEST_TRUE(validator.validate(good).valid)
}
END_SECTION

START_SECTION(([EXTRA] a table that does not carry the schema of its view is refused))
{
  auto psm_table = makePSMTable();
  // The pg view's contract cannot be checked against a psm table; a mismatch must be reported,
  // not silently accepted by a validator that finds none of its columns.
  QPXValueValidation validator(QPXValueValidation::View::PROTEIN_GROUP);
  const auto result = validator.validate(psm_table);
  TEST_FALSE(result.valid)
  TEST_TRUE(result.errors.size() > 0)
}
END_SECTION

START_SECTION(([EXTRA] a multi-chunk column is validated as one logical column))
{
  // Tables handed to validate() are not always single-chunk (a read-modify-write, or a
  // concatenation of batches), and a duplicate that straddles two chunks is still a duplicate.
  auto table = makePSMTable(1);
  auto doubled = arrow::ConcatenateTables({table, table}).ValueOrDie();
  TEST_EQUAL(doubled->num_rows(), 2)
  TEST_TRUE(doubled->column(0)->num_chunks() > 1)

  QPXValueValidation validator(QPXValueValidation::View::PSM);
  const auto result = validator.validate(doubled);
  TEST_FALSE(result.valid)
  TEST_TRUE(result.toString().find("primary key") != std::string::npos)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
