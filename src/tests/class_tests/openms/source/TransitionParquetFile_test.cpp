// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryIDNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

#include <map>
#include <set>
#include <vector>

using namespace OpenMS;
using namespace std;

namespace
{
  template <typename BuilderT, typename ValueT>
  void appendOk_(BuilderT& builder, const ValueT& value)
  {
    TEST_EQUAL(builder.Append(value).ok(), true)
  }

  template <typename BuilderT>
  std::shared_ptr<arrow::Array> finishArray_(BuilderT& builder)
  {
    auto result = builder.Finish();
    TEST_EQUAL(result.ok(), true)
    if (result.ok())
    {
      return result.ValueOrDie();
    }
    return nullptr;
  }

  ::arrow::Status writeParquetTable_(const std::shared_ptr<arrow::Table>& table, const std::string& filename)
  {
    auto outfile_result = arrow::io::FileOutputStream::Open(std::string(filename));
    if (!outfile_result.ok())
    {
      return outfile_result.status();
    }
    auto outfile = outfile_result.ValueOrDie();
    return parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 1024);
  }

  std::string joinProteinAccessions_(const std::vector<std::string>& accessions)
  {
    std::string joined;
    for (Size i = 0; i < accessions.size(); ++i)
    {
      if (i > 0) joined += ";";
      joined += accessions[i];
    }
    return std::string(joined);
  }
}

START_TEST(TransitionParquetFile, "$Id$")

TransitionParquetFile* ptr = nullptr;
TransitionParquetFile* nullPointer = nullptr;

START_SECTION(TransitionParquetFile())
{
  ptr = new TransitionParquetFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~TransitionParquetFile())
{
  delete ptr;
}
END_SECTION

START_SECTION(void convertParquetToTargetedExperiment(const std::string& oswpq_dir, OpenSwath::LightTargetedExperiment& targeted_exp) const)
{
  const std::string input_file = OPENMS_GET_TEST_DATA_PATH("MRMAssay_detectingTransistionCompound_input.TraML");
  TraMLFile traml;
  TargetedExperiment targeted_exp;
  traml.load(input_file, targeted_exp);

  OpenSwath::LightTargetedExperiment light_exp;
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);

  Size compound_count = std::min<Size>(2, light_exp.compounds.size());
  TEST_EQUAL(compound_count, 2)

  std::map<std::string, int64_t> compound_to_precursor;
  int64_t precursor_id = 1;
  for (Size i = 0; i < compound_count; ++i)
  {
    compound_to_precursor[light_exp.compounds[i].id] = precursor_id++;
  }

  std::vector<OpenSwath::LightTransition> transitions;
  std::map<std::string, double> precursor_mz;
  for (const auto& transition : light_exp.transitions)
  {
    if (compound_to_precursor.find(transition.peptide_ref) != compound_to_precursor.end())
    {
      transitions.push_back(transition);
      if (precursor_mz.find(transition.peptide_ref) == precursor_mz.end())
      {
        precursor_mz[transition.peptide_ref] = transition.precursor_mz;
      }
    }
  }
  TEST_EQUAL(transitions.size() > 0, true)

  File::TempDir tmp_dir;
  const std::string base_dir = tmp_dir.getPath() + "/test.oswpq";
  const std::string library_dir = base_dir + "/library";
  File::makeDir(base_dir);
  File::makeDir(library_dir);

  // Build precursors table
  arrow::Int64Builder precursor_id_builder;
  arrow::DoubleBuilder precursor_mz_builder;
  arrow::Int32Builder precursor_charge_builder;
  arrow::DoubleBuilder library_rt_builder;
  arrow::DoubleBuilder drift_time_builder;
  arrow::BooleanBuilder decoy_builder;
  arrow::StringBuilder traml_id_builder;
  arrow::StringBuilder modified_sequence_builder;
  arrow::StringBuilder unmodified_sequence_builder;
  arrow::StringBuilder protein_accessions_builder;

  for (Size i = 0; i < compound_count; ++i)
  {
    const auto& compound = light_exp.compounds[i];
    const std::string compound_id = compound.id;
    const int64_t id = compound_to_precursor[compound_id];

    appendOk_(precursor_id_builder, id);
    appendOk_(precursor_mz_builder, precursor_mz[compound_id]);
    appendOk_(precursor_charge_builder, compound.charge);
    appendOk_(library_rt_builder, compound.rt);
    appendOk_(drift_time_builder, compound.drift_time);
    appendOk_(decoy_builder, false);
    appendOk_(traml_id_builder, std::string(compound.id));
    appendOk_(modified_sequence_builder, std::string(compound.sequence));

    std::string unmodified_sequence;
    if (!compound.sequence.empty())
    {
      try
      {
        unmodified_sequence = AASequence::fromString(compound.sequence).toUnmodifiedString();
      }
      catch (Exception::InvalidValue&)
      {
        unmodified_sequence = "";
      }
    }
    appendOk_(unmodified_sequence_builder, std::string(unmodified_sequence));
    appendOk_(protein_accessions_builder, joinProteinAccessions_(compound.protein_refs));
  }

  auto precursor_id_array = finishArray_(precursor_id_builder);
  auto precursor_mz_array = finishArray_(precursor_mz_builder);
  auto precursor_charge_array = finishArray_(precursor_charge_builder);
  auto library_rt_array = finishArray_(library_rt_builder);
  auto drift_time_array = finishArray_(drift_time_builder);
  auto decoy_array = finishArray_(decoy_builder);
  auto traml_id_array = finishArray_(traml_id_builder);
  auto modified_sequence_array = finishArray_(modified_sequence_builder);
  auto unmodified_sequence_array = finishArray_(unmodified_sequence_builder);
  auto protein_accessions_array = finishArray_(protein_accessions_builder);

  auto precursor_schema = arrow::schema({
    arrow::field("precursor_id", arrow::int64()),
    arrow::field("precursor_mz", arrow::float64()),
    arrow::field("charge", arrow::int32()),
    arrow::field("library_rt", arrow::float64()),
    arrow::field("library_drift_time", arrow::float64()),
    arrow::field("decoy", arrow::boolean()),
    arrow::field("traml_id", arrow::utf8()),
    arrow::field("modified_sequence", arrow::utf8()),
    arrow::field("unmodified_sequence", arrow::utf8()),
    arrow::field("protein_accessions", arrow::utf8())
  });

  auto precursors_table = arrow::Table::Make(
    precursor_schema,
    {precursor_id_array, precursor_mz_array, precursor_charge_array, library_rt_array,
     drift_time_array, decoy_array, traml_id_array, modified_sequence_array,
     unmodified_sequence_array, protein_accessions_array});

  TEST_EQUAL(writeParquetTable_(precursors_table, library_dir + "/precursors.parquet").ok(), true)

  // Build transitions table
  arrow::Int64Builder transition_id_builder;
  arrow::Int64Builder transition_precursor_id_builder;
  arrow::StringBuilder transition_traml_id_builder;
  arrow::DoubleBuilder product_mz_builder;
  arrow::Int32Builder fragment_charge_builder;
  arrow::StringBuilder fragment_type_builder;
  arrow::StringBuilder annotation_builder;
  arrow::Int32Builder ordinal_builder;
  arrow::BooleanBuilder detecting_builder;
  arrow::BooleanBuilder identifying_builder;
  arrow::BooleanBuilder quantifying_builder;
  arrow::DoubleBuilder transition_intensity_builder;
  arrow::BooleanBuilder transition_decoy_builder;

  int64_t transition_id = 1;
  for (const auto& transition : transitions)
  {
    const int64_t precursor_ref = compound_to_precursor[transition.peptide_ref];
    appendOk_(transition_id_builder, transition_id++);
    appendOk_(transition_precursor_id_builder, precursor_ref);
    appendOk_(transition_traml_id_builder, std::string(transition.transition_name));
    appendOk_(product_mz_builder, transition.product_mz);
    appendOk_(fragment_charge_builder, transition.fragment_charge);
    appendOk_(fragment_type_builder, std::string(transition.getFragmentType()));
    appendOk_(annotation_builder, std::string(transition.getAnnotation()));
    appendOk_(ordinal_builder, transition.fragment_nr);
    appendOk_(detecting_builder, transition.isDetectingTransition());
    appendOk_(identifying_builder, transition.isIdentifyingTransition());
    appendOk_(quantifying_builder, transition.isQuantifyingTransition());
    appendOk_(transition_intensity_builder, transition.library_intensity);
    appendOk_(transition_decoy_builder, transition.getDecoy());
  }

  auto transition_id_array = finishArray_(transition_id_builder);
  auto transition_precursor_id_array = finishArray_(transition_precursor_id_builder);
  auto transition_traml_id_array = finishArray_(transition_traml_id_builder);
  auto product_mz_array = finishArray_(product_mz_builder);
  auto fragment_charge_array = finishArray_(fragment_charge_builder);
  auto fragment_type_array = finishArray_(fragment_type_builder);
  auto annotation_array = finishArray_(annotation_builder);
  auto ordinal_array = finishArray_(ordinal_builder);
  auto detecting_array = finishArray_(detecting_builder);
  auto identifying_array = finishArray_(identifying_builder);
  auto quantifying_array = finishArray_(quantifying_builder);
  auto transition_intensity_array = finishArray_(transition_intensity_builder);
  auto transition_decoy_array = finishArray_(transition_decoy_builder);

  auto transition_schema = arrow::schema({
    arrow::field("transition_id", arrow::int64()),
    arrow::field("precursor_id", arrow::int64()),
    arrow::field("traml_id", arrow::utf8()),
    arrow::field("product_mz", arrow::float64()),
    arrow::field("charge", arrow::int32()),
    arrow::field("type", arrow::utf8()),
    arrow::field("annotation", arrow::utf8()),
    arrow::field("ordinal", arrow::int32()),
    arrow::field("detecting", arrow::boolean()),
    arrow::field("identifying", arrow::boolean()),
    arrow::field("quantifying", arrow::boolean()),
    arrow::field("library_intensity", arrow::float64()),
    arrow::field("decoy", arrow::boolean())
  });

  auto transitions_table = arrow::Table::Make(
    transition_schema,
    {transition_id_array, transition_precursor_id_array, transition_traml_id_array,
     product_mz_array, fragment_charge_array, fragment_type_array, annotation_array,
     ordinal_array, detecting_array, identifying_array, quantifying_array,
     transition_intensity_array, transition_decoy_array});

  TEST_EQUAL(writeParquetTable_(transitions_table, library_dir + "/transitions.parquet").ok(), true)

  TransitionParquetFile reader;
  OpenSwath::LightTargetedExperiment out_exp;
  OpenSwathLibraryIDNormalizer::SourceIDMapping source_ids;
  reader.convertParquetToTargetedExperiment(base_dir, out_exp, &source_ids);

  TEST_EQUAL(out_exp.compounds.size(), compound_count)
  TEST_EQUAL(out_exp.transitions.size(), transitions.size())

  // OSWPQ persistent numeric IDs define operational identity on read. traml_id
  // remains source/provenance metadata and must not replace those foreign keys.
  std::set<std::string> expected_precursor_ids;
  for (const auto& entry : compound_to_precursor)
  {
    expected_precursor_ids.insert(std::to_string(entry.second));
  }

  for (const auto& entry : compound_to_precursor)
  {
    const std::string canonical_id = std::to_string(entry.second);
    TEST_EQUAL(source_ids.precursor_source_to_canonical.at(entry.first), canonical_id)
    TEST_EQUAL(source_ids.precursor_canonical_to_source.at(canonical_id), entry.first)
  }

  std::set<std::string> compound_refs;
  for (const auto& compound : out_exp.compounds)
  {
    TEST_EQUAL(expected_precursor_ids.count(compound.id), 1)
    compound_refs.insert(compound.id);
  }
  TEST_EQUAL(compound_refs.size(), expected_precursor_ids.size())

  for (Size i = 0; i < out_exp.transitions.size(); ++i)
  {
    TEST_EQUAL(out_exp.transitions[i].transition_name, std::to_string(i + 1))
    TEST_EQUAL(out_exp.transitions[i].peptide_ref,
               std::to_string(compound_to_precursor[transitions[i].peptide_ref]))
    TEST_EQUAL(compound_refs.count(out_exp.transitions[i].peptide_ref), 1)
  }

  // Reader-boundary regression: duplicate persistent precursor IDs must be
  // rejected before unordered_map materialization can silently discard one row.
  arrow::Int64Builder duplicate_precursor_id_builder;
  for (Size i = 0; i < compound_count; ++i)
  {
    appendOk_(duplicate_precursor_id_builder, static_cast<int64_t>(1));
  }
  auto duplicate_precursor_id_array = finishArray_(duplicate_precursor_id_builder);
  auto duplicate_precursors_table = arrow::Table::Make(
    precursor_schema,
    {duplicate_precursor_id_array, precursor_mz_array, precursor_charge_array, library_rt_array,
     drift_time_array, decoy_array, traml_id_array, modified_sequence_array,
     unmodified_sequence_array, protein_accessions_array});

  // Keep all transition foreign keys valid so the duplicated precursor ID is the
  // only invariant violated by this fixture.
  arrow::Int64Builder duplicate_transition_precursor_id_builder;
  for (Size i = 0; i < transitions.size(); ++i)
  {
    appendOk_(duplicate_transition_precursor_id_builder, static_cast<int64_t>(1));
  }
  auto duplicate_transition_precursor_id_array = finishArray_(duplicate_transition_precursor_id_builder);
  auto duplicate_transitions_table = arrow::Table::Make(
    transition_schema,
    {transition_id_array, duplicate_transition_precursor_id_array, transition_traml_id_array,
     product_mz_array, fragment_charge_array, fragment_type_array, annotation_array,
     ordinal_array, detecting_array, identifying_array, quantifying_array,
     transition_intensity_array, transition_decoy_array});

  const std::string duplicate_dir = tmp_dir.getPath() + "/duplicate.oswpq";
  const std::string duplicate_library_dir = duplicate_dir + "/library";
  File::makeDir(duplicate_dir);
  File::makeDir(duplicate_library_dir);
  TEST_EQUAL(writeParquetTable_(duplicate_precursors_table, duplicate_library_dir + "/precursors.parquet").ok(), true)
  TEST_EQUAL(writeParquetTable_(duplicate_transitions_table, duplicate_library_dir + "/transitions.parquet").ok(), true)

  OpenSwath::LightTargetedExperiment duplicate_out;
  TEST_EXCEPTION(Exception::InvalidValue,
                 reader.convertParquetToTargetedExperiment(duplicate_dir, duplicate_out))
}
END_SECTION

START_SECTION(void convertLightTargetedExperimentToParquet(const std::string& oswpq_path, const OpenSwath::LightTargetedExperiment& targeted_exp) const)
{
  // --- Build a reference LightTargetedExperiment from a TraML file ---
  const std::string input_file = OPENMS_GET_TEST_DATA_PATH("MRMAssay_detectingTransistionCompound_input.TraML");
  TraMLFile traml;
  TargetedExperiment targeted_exp;
  traml.load(input_file, targeted_exp);

  OpenSwath::LightTargetedExperiment light_exp;
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);
  TEST_EQUAL(light_exp.compounds.size() > 0, true)
  TEST_EQUAL(light_exp.transitions.size() > 0, true)

  // --- Write to a temporary .oswpq directory ---
  File::TempDir tmp_dir;
  const std::string out_dir = tmp_dir.getPath() + "/roundtrip.oswpq";
  File::makeDir(out_dir);

  const auto source_ids = OpenSwathLibraryIDNormalizer::normalizeSourceIDs(light_exp);

  TransitionParquetFile writer;
  writer.convertLightTargetedExperimentToParquet(out_dir, light_exp, &source_ids);

  // Verify library files exist
  TEST_EQUAL(File::exists(out_dir + "/library/precursors.parquet"), true)
  TEST_EQUAL(File::exists(out_dir + "/library/transitions.parquet"), true)
  TEST_EQUAL(File::exists(out_dir + "/library/metadata.json"), true)

  // --- Read back and compare ---
  TransitionParquetFile reader;
  OpenSwath::LightTargetedExperiment roundtrip_exp;
  OpenSwathLibraryIDNormalizer::SourceIDMapping roundtrip_source_ids;
  reader.convertParquetToTargetedExperiment(out_dir, roundtrip_exp, &roundtrip_source_ids);

  TEST_EQUAL(roundtrip_source_ids.precursor_source_to_canonical.size(), source_ids.precursor_source_to_canonical.size())
  for (const auto& [source_id, canonical_id] : source_ids.precursor_source_to_canonical)
  {
    TEST_EQUAL(roundtrip_source_ids.precursor_source_to_canonical.at(source_id), canonical_id)
  }
  TEST_EQUAL(roundtrip_source_ids.transition_canonical_to_source.size(), source_ids.transition_canonical_to_source.size())
  for (const auto& [canonical_id, source_id] : source_ids.transition_canonical_to_source)
  {
    TEST_EQUAL(roundtrip_source_ids.transition_canonical_to_source.at(canonical_id), source_id)
  }

  // The persisted OSWPQ IDs must remain valid canonical operational IDs after
  // the writer -> reader round trip.
  OpenSwathLibraryIDNormalizer::validateCanonicalIDs(roundtrip_exp);

  TEST_EQUAL(roundtrip_exp.compounds.size(), light_exp.compounds.size())
  TEST_EQUAL(roundtrip_exp.transitions.size(), light_exp.transitions.size())
  TEST_EQUAL(roundtrip_exp.proteins.size(), light_exp.proteins.size())

  // Verify each transition references a valid compound
  std::set<std::string> roundtrip_compound_ids;
  for (const auto& compound : roundtrip_exp.compounds)
  {
    roundtrip_compound_ids.insert(compound.id);
  }
  for (const auto& transition : roundtrip_exp.transitions)
  {
    TEST_EQUAL(roundtrip_compound_ids.count(transition.peptide_ref) > 0, true)
  }

  // Source normalization assigned canonical transition IDs in row order. The
  // writer preserves them exactly and the reader must not restore transition traml_id.
  for (Size i = 0; i < roundtrip_exp.transitions.size(); ++i)
  {
    TEST_EQUAL(roundtrip_exp.transitions[i].transition_name, std::to_string(i))
    TEST_REAL_SIMILAR(roundtrip_exp.transitions[i].product_mz, light_exp.transitions[i].product_mz)
  }

  // A canonical-looking but malformed experiment must be rejected instead of being
  // silently renumbered by the compatibility writer.
  TEST_EQUAL(light_exp.transitions.size() > 1, true)
  OpenSwath::LightTargetedExperiment malformed = light_exp;
  malformed.transitions[1].transition_name = malformed.transitions[0].transition_name;
  const std::string malformed_out = tmp_dir.getPath() + "/malformed.oswpq";
  TEST_EXCEPTION(Exception::InvalidValue,
                 writer.convertLightTargetedExperimentToParquet(malformed_out, malformed))
}
END_SECTION

END_TEST
