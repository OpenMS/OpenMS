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
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

#include <map>
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

  ::arrow::Status writeParquetTable_(const std::shared_ptr<arrow::Table>& table, const String& filename)
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
    String joined;
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

START_SECTION(void convertParquetToTargetedExperiment(const String& oswpq_dir, OpenSwath::LightTargetedExperiment& targeted_exp) const)
{
  const String input_file = OPENMS_GET_TEST_DATA_PATH("MRMAssay_detectingTransistionCompound_input.TraML");
  TraMLFile traml;
  TargetedExperiment targeted_exp;
  traml.load(input_file, targeted_exp);

  OpenSwath::LightTargetedExperiment light_exp;
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);

  Size compound_count = std::min<Size>(2, light_exp.compounds.size());
  TEST_EQUAL(compound_count > 0, true)

  std::map<String, int64_t> compound_to_precursor;
  int64_t precursor_id = 1;
  for (Size i = 0; i < compound_count; ++i)
  {
    compound_to_precursor[light_exp.compounds[i].id] = precursor_id++;
  }

  std::vector<OpenSwath::LightTransition> transitions;
  std::map<String, double> precursor_mz;
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
  const String base_dir = tmp_dir.getPath() + "/test.oswpq";
  const String library_dir = base_dir + "/library";
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
    const String compound_id = compound.id;
    const int64_t id = compound_to_precursor[compound_id];

    appendOk_(precursor_id_builder, id);
    appendOk_(precursor_mz_builder, precursor_mz[compound_id]);
    appendOk_(precursor_charge_builder, compound.charge);
    appendOk_(library_rt_builder, compound.rt);
    appendOk_(drift_time_builder, compound.drift_time);
    appendOk_(decoy_builder, false);
    appendOk_(traml_id_builder, std::string(compound.id));
    appendOk_(modified_sequence_builder, std::string(compound.sequence));

    String unmodified_sequence;
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
  reader.convertParquetToTargetedExperiment(base_dir, out_exp);

  TEST_EQUAL(out_exp.compounds.size(), compound_count)
  TEST_EQUAL(out_exp.transitions.size(), transitions.size())

  std::map<String, int> compound_refs;
  for (const auto& compound : out_exp.compounds)
  {
    compound_refs[compound.id] = 1;
  }
  for (const auto& transition : out_exp.transitions)
  {
    TEST_EQUAL(compound_refs.find(transition.peptide_ref) != compound_refs.end(), true)
  }
}
END_SECTION

START_SECTION(void convertLightTargetedExperimentToParquet(const String& oswpq_path, const OpenSwath::LightTargetedExperiment& targeted_exp) const)
{
  // --- Build a reference LightTargetedExperiment from a TraML file ---
  const String input_file = OPENMS_GET_TEST_DATA_PATH("MRMAssay_detectingTransistionCompound_input.TraML");
  TraMLFile traml;
  TargetedExperiment targeted_exp;
  traml.load(input_file, targeted_exp);

  OpenSwath::LightTargetedExperiment light_exp;
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);
  TEST_EQUAL(light_exp.compounds.size() > 0, true)
  TEST_EQUAL(light_exp.transitions.size() > 0, true)

  // --- Write to a temporary .oswpq directory ---
  File::TempDir tmp_dir;
  const String out_dir = tmp_dir.getPath() + "/roundtrip.oswpq";
  File::makeDir(out_dir);

  TransitionParquetFile writer;
  writer.convertLightTargetedExperimentToParquet(out_dir, light_exp);

  // Verify library files exist
  TEST_EQUAL(File::exists(out_dir + "/library/precursors.parquet"), true)
  TEST_EQUAL(File::exists(out_dir + "/library/transitions.parquet"), true)
  TEST_EQUAL(File::exists(out_dir + "/library/metadata.json"), true)

  // --- Read back and compare ---
  TransitionParquetFile reader;
  OpenSwath::LightTargetedExperiment roundtrip_exp;
  reader.convertParquetToTargetedExperiment(out_dir, roundtrip_exp);

  TEST_EQUAL(roundtrip_exp.compounds.size(), light_exp.compounds.size())
  TEST_EQUAL(roundtrip_exp.transitions.size(), light_exp.transitions.size())
  TEST_EQUAL(roundtrip_exp.proteins.size(), light_exp.proteins.size())

  // Verify each transition references a valid compound
  std::set<String> roundtrip_compound_ids;
  for (const auto& compound : roundtrip_exp.compounds)
  {
    roundtrip_compound_ids.insert(compound.id);
  }
  for (const auto& transition : roundtrip_exp.transitions)
  {
    TEST_EQUAL(roundtrip_compound_ids.count(transition.peptide_ref) > 0, true)
  }

  // Verify transition product_mz values are preserved
  // Build a map from transition name -> product_mz for the original
  std::map<String, double> orig_transition_mz;
  for (const auto& tr : light_exp.transitions)
  {
    orig_transition_mz[tr.transition_name] = tr.product_mz;
  }
  for (const auto& tr : roundtrip_exp.transitions)
  {
    auto it = orig_transition_mz.find(tr.transition_name);
    if (it != orig_transition_mz.end())
    {
      TEST_REAL_SIMILAR(tr.product_mz, it->second)
    }
  }
}
END_SECTION

END_TEST
