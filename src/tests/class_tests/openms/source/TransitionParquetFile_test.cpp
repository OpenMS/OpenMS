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

#ifdef WITH_PARQUET
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#endif

#include <map>
#include <vector>

using namespace OpenMS;
using namespace std;

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

#ifdef WITH_PARQUET
namespace
{
  ::arrow::Status writeParquetTable_(const std::shared_ptr<arrow::Table>& table, const String& filename)
  {
    ARROW_ASSIGN_OR_RAISE(auto outfile, arrow::io::FileOutputStream::Open(filename.toStdString()));
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
    return joined;
  }
}
#endif

START_SECTION(void convertParquetToTargetedExperiment(const String& pqp_parquet_dir, OpenSwath::LightTargetedExperiment& targeted_exp) const)
{
#ifdef WITH_PARQUET
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
  const String base_dir = tmp_dir.getPath() + "/test.pqp_parquet";
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

    precursor_id_builder.Append(id);
    precursor_mz_builder.Append(precursor_mz[compound_id]);
    precursor_charge_builder.Append(compound.charge);
    library_rt_builder.Append(compound.rt);
    drift_time_builder.Append(compound.drift_time);
    decoy_builder.Append(false);
    traml_id_builder.Append(compound.id);
    modified_sequence_builder.Append(compound.sequence);

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
    unmodified_sequence_builder.Append(unmodified_sequence);
    protein_accessions_builder.Append(joinProteinAccessions_(compound.protein_refs));
  }

  std::shared_ptr<arrow::Array> precursor_id_array;
  std::shared_ptr<arrow::Array> precursor_mz_array;
  std::shared_ptr<arrow::Array> precursor_charge_array;
  std::shared_ptr<arrow::Array> library_rt_array;
  std::shared_ptr<arrow::Array> drift_time_array;
  std::shared_ptr<arrow::Array> decoy_array;
  std::shared_ptr<arrow::Array> traml_id_array;
  std::shared_ptr<arrow::Array> modified_sequence_array;
  std::shared_ptr<arrow::Array> unmodified_sequence_array;
  std::shared_ptr<arrow::Array> protein_accessions_array;

  ARROW_ASSIGN_OR_RAISE(precursor_id_array, precursor_id_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(precursor_mz_array, precursor_mz_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(precursor_charge_array, precursor_charge_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(library_rt_array, library_rt_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(drift_time_array, drift_time_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(decoy_array, decoy_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(traml_id_array, traml_id_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(modified_sequence_array, modified_sequence_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(unmodified_sequence_array, unmodified_sequence_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(protein_accessions_array, protein_accessions_builder.Finish());

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
    transition_id_builder.Append(transition_id++);
    transition_precursor_id_builder.Append(precursor_ref);
    transition_traml_id_builder.Append(transition.transition_name);
    product_mz_builder.Append(transition.product_mz);
    fragment_charge_builder.Append(transition.fragment_charge);
    fragment_type_builder.Append(transition.getFragmentType());
    annotation_builder.Append(transition.getAnnotation());
    ordinal_builder.Append(transition.fragment_nr);
    detecting_builder.Append(transition.isDetectingTransition());
    identifying_builder.Append(transition.isIdentifyingTransition());
    quantifying_builder.Append(transition.isQuantifyingTransition());
    transition_intensity_builder.Append(transition.library_intensity);
    transition_decoy_builder.Append(transition.getDecoy());
  }

  std::shared_ptr<arrow::Array> transition_id_array;
  std::shared_ptr<arrow::Array> transition_precursor_id_array;
  std::shared_ptr<arrow::Array> transition_traml_id_array;
  std::shared_ptr<arrow::Array> product_mz_array;
  std::shared_ptr<arrow::Array> fragment_charge_array;
  std::shared_ptr<arrow::Array> fragment_type_array;
  std::shared_ptr<arrow::Array> annotation_array;
  std::shared_ptr<arrow::Array> ordinal_array;
  std::shared_ptr<arrow::Array> detecting_array;
  std::shared_ptr<arrow::Array> identifying_array;
  std::shared_ptr<arrow::Array> quantifying_array;
  std::shared_ptr<arrow::Array> transition_intensity_array;
  std::shared_ptr<arrow::Array> transition_decoy_array;

  ARROW_ASSIGN_OR_RAISE(transition_id_array, transition_id_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(transition_precursor_id_array, transition_precursor_id_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(transition_traml_id_array, transition_traml_id_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(product_mz_array, product_mz_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(fragment_charge_array, fragment_charge_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(fragment_type_array, fragment_type_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(annotation_array, annotation_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(ordinal_array, ordinal_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(detecting_array, detecting_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(identifying_array, identifying_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(quantifying_array, quantifying_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(transition_intensity_array, transition_intensity_builder.Finish());
  ARROW_ASSIGN_OR_RAISE(transition_decoy_array, transition_decoy_builder.Finish());

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
#else
  NOT_TESTABLE
#endif
}
END_SECTION

END_TEST
