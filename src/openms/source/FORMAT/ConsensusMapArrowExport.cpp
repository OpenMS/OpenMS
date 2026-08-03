// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>

#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/QPXValueValidation.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/Scores.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/METADATA/IdentifierMSRunMapper.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/METADATA/SpectrumNativeIDParser.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <algorithm>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

namespace // anonymous
{
  /// Try to extract a scan number from a native ID string.
  /// Uses SpectrumNativeIDParser to auto-detect the native ID format.
  /// Returns -1 on failure.
  Int extractScan(const std::string& native_id)
  {
    std::string regex_str = SpectrumNativeIDParser::getRegExFromNativeID(native_id);
    if (regex_str.empty())
    {
      return -1;
    }
    boost::regex scan_regex(regex_str);
    return SpectrumNativeIDParser::extractScanNumber(native_id, scan_regex, true);
  }

  /// Resolve a MetaInfo index to its registered name, caching results in @p cache.
  /// MetaInfoRegistry::getName() takes an `omp critical` lock, so caching by index pays
  /// that lock at most once per unique key for the whole range instead of once per
  /// metavalue per feature (as getKeys()/getMetaValue() would). Registered names are
  /// non-empty, so an empty cache slot reliably means "not yet resolved".
  inline const std::string& cachedMetaName_(UInt index, std::vector<std::string>& cache)
  {
    if (index >= cache.size()) { cache.resize(static_cast<size_t>(index) + 1); }
    std::string& name = cache[index];
    if (name.empty()) { name = MetaInfoInterface::metaRegistry().getName(index); }
    return name;
  }
} // anonymous namespace


namespace // anonymous (feature range builder + shared per-map lookups)
{

/// Per-map, read-only lookups derived once from the ConsensusMap's protein
/// identifications. Shared by the whole-table and streaming feature export paths
/// (built once, reused across all batches).
struct FeatureIdLookups
{
  std::unordered_map<std::string, double> pg_qvalue;               ///< accession -> protein-group q-value
  std::unordered_map<std::string, std::vector<std::string>> gg_accessions; ///< accession -> gene accessions
  std::unordered_map<std::string, std::vector<std::string>> gg_names;      ///< accession -> gene names
  /// Resolves a PeptideIdentification to its own source run via id_merge_index, for
  /// id_run_file_name. Read-only after construction, so it is safe to share across the
  /// streaming path's OpenMP workers.
  IdentifierMSRunMapper run_mapper;

  /// Per column header (map index): the QPX run name and intensity label.
  /// Both are expensive to derive and depend only on the header, of which there are at most
  /// a few dozen -- while the export loop touches them once per FeatureHandle per row.
  /// File::stemName() scans the 73-entry file-type table, and reading the header's
  /// channel_name by string takes the MetaInfoRegistry omp-critical lock, which would
  /// serialize the parallel batch build (see the note on pre-resolved metavalue indices).
  struct RunColumnInfo
  {
    std::string source;     ///< ColumnHeader::filename, the run's identity
    std::string run_name;   ///< stemmed run_file_name, the emitted value
    std::string label;      ///< QPX intensities[].label
  };
  std::unordered_map<UInt64, RunColumnInfo> run_columns;

  /// Whether the map is isobaric, deciding where a row's coordinates come from.
  bool is_isobaric = false;
  /// Set when a column header names a channel QPX cannot label; the export must not proceed,
  /// since the label is a join key.
  bool has_unlabelable_channel = false;
};

/// Build the per-map protein-group q-value and gene-info lookups from the
/// ConsensusMap's ProteinIdentifications.
FeatureIdLookups buildFeatureIdLookups(const ConsensusMap& cmap)
{
  FeatureIdLookups lut;

  // Deliberately OUTSIDE the try below: that catch exists to tolerate a duplicate PATH LIST,
  // for which the forward map is still complete. A duplicate IDENTIFIER is the opposite -- the
  // forward map assigns with operator[], so only the last run's paths survive and a feature of
  // the first would be stamped with a file it never came from, silently, since a single row
  // gives nothing away. The pg exporter refuses the same input.
  {
    std::set<std::string> seen_identifiers;
    for (const auto& prot_id : cmap.getProteinIdentifications())
    {
      if (!seen_identifiers.insert(prot_id.getIdentifier()).second)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "ConsensusMapArrowExport: several identification runs share the identifier '"
          + prot_id.getIdentifier() + "'. Their MS run paths cannot be told apart, so a feature "
          "would be attributed to the wrong run file.", prot_id.getIdentifier());
      }
    }
  }

  // The mapper's constructor throws when two identifiers share a path list; that guards only
  // the reverse map, which is never used here, and the forward map is fully populated before
  // the throw. Degrade to a warning rather than failing an otherwise-exportable map.
  try
  {
    lut.run_mapper.create(cmap.getProteinIdentifications());
  }
  catch (const Exception::InvalidValue& e)
  {
    OPENMS_LOG_WARN << "ConsensusMapArrowExport: ambiguous MS run paths across identification "
                       "runs (" << e.getMessage() << "); id_run_file_name resolution continues "
                       "on the identifier -> path mapping." << std::endl;
  }
  // A map is isobaric if it says so OR if its headers carry channel annotations. Neither
  // signal alone is reliable: IsobaricWorkflow annotates channels without ever calling
  // setExperimentType, while other producers declare a labelled experiment type without
  // writing channel_name. Getting this wrong writes reporter-ion m/z into observed_mz, so
  // accept either. (The sibling pg exporter reaches the same conclusion via
  // ColumnHeader::getLabelAsUInt(getExperimentType()).)
  lut.is_isobaric = (cmap.getExperimentType() != "label-free");

  // Same policy as the PSM exporter: two source paths sharing a stem are indistinguishable
  // as a join key, but same-named files in different directories are a legitimate layout, so
  // warn rather than fail. Rows stay separate -- grouping keys on the source path, not the
  // stem -- so this reports an ambiguous key, never a merged one.
  {
    std::vector<std::string> header_paths;
    header_paths.reserve(cmap.getColumnHeaders().size());
    for (const auto& [map_index, header] : cmap.getColumnHeaders()) { header_paths.push_back(header.filename); }
    (void)ArrowIOHelpers::qpxWarnOnRunNameCollisions("ConsensusMapArrowExport", header_paths);
  }
  const auto qpx_labels = ArrowIOHelpers::qpxIntensityLabels(cmap);
  for (const auto& [map_index, header] : cmap.getColumnHeaders())
  {
    const std::string channel_name = header.metaValueExists("channel_name")
                                   ? header.getMetaValue("channel_name").toString() : "";
    if (!channel_name.empty()) { lut.is_isobaric = true; }

    std::string label = qpx_labels.at(map_index);
    if (label.empty())
    {
      // qpxIntensityLabels() documents that callers must handle this: the channel's method
      // could not be identified, so any token written here would be a guessed join key.
      lut.has_unlabelable_channel = true;
      OPENMS_LOG_ERROR << "ConsensusMapArrowExport: column header '" << header.label
                       << "' names channel '" << channel_name << "' whose quantitation method "
                          "cannot be identified, so no QPX intensity label can be derived for it."
                       << std::endl;
    }
    lut.run_columns[map_index] = FeatureIdLookups::RunColumnInfo{
      header.filename, ArrowIOHelpers::qpxRunFileName(header.filename), std::move(label)};
  }

  for (const auto& prot_id : cmap.getProteinIdentifications())
  {
    // Gene info from ProteinHit metavalues
    for (const auto& ph : prot_id.getHits())
    {
      const std::string& acc = ph.getAccession();
      if (ph.metaValueExists("gene_accession"))
      {
        lut.gg_accessions[acc].push_back(ph.getMetaValue("gene_accession").toString());
      }
      if (ph.metaValueExists("gene_name"))
      {
        lut.gg_names[acc].push_back(ph.getMetaValue("gene_name").toString());
      }
    }

    for (const auto& pg : prot_id.getProteinGroups())
    {
      for (const auto& acc : pg.accessions)
      {
        lut.pg_qvalue[acc] = pg.probability;
      }
    }
    for (const auto& ig : prot_id.getIndistinguishableProteins())
    {
      for (const auto& acc : ig.accessions)
      {
        if (!lut.pg_qvalue.contains(acc))
        {
          lut.pg_qvalue[acc] = ig.probability;
        }
      }
    }
  }
  return lut;
}

/// Select the identity of a consensus feature: the best-scoring hit across ALL of its
/// peptide identifications, together with the index of the identification holding it.
///
/// Deliberately local rather than an IDFilter addition. IDFilter::getBestHit() reports only
/// the hit, and its internal "best" iterator is the score-type *reference* (pinned to the
/// first non-empty identification), not the winner -- while the QPX feature view needs the
/// winner's parent to resolve id_run_file_name. Exposing that would mean new public API in a
/// widely-included header for a single caller, and it could not replace the existing
/// per-identification helpers (MapAlignmentAlgorithmIdentification::getBestScoringHit is
/// protected, works within one identification and returns a pointer; MzTab aggregates scores
/// with different tie-breaking). Score comparison itself is not duplicated -- it delegates to
/// PeptideIdentification::getScoreComparator, the same primitive those helpers use.
///
/// Comparison uses ONE orientation, taken from the first identification that has a hit. Using
/// each identification's own orientation would compare in different directions within a single
/// selection and make the winner depend on identification order.
///
/// Identifications that disagree on score type or orientation are a pipeline defect, not a
/// case to recover from: scores of different types are not comparable, harmonizing them is an
/// explicit step (IDScoreSwitcher / IDScoreSwitcherAlgorithm, which ProteomicsLFQ, Epifany and
/// FalseDiscoveryRate all run), and IDFilter::getBestHit throws outright on the same condition.
/// This warns via @p inconsistent_scores and proceeds with the reference orientation rather
/// than inventing a recovery, which would only turn a broken input into an arbitrary answer.
///
/// @param[in] pep_ids The feature's peptide identifications
/// @param[out] best_id_index Index into @p pep_ids of the identification holding the winner
/// @param[out] inconsistent_scores Set when the identifications disagree on score type or orientation
/// @return Pointer to the winning hit, or nullptr when no identification has one
const PeptideHit* selectFeatureIdentity(
  const PeptideIdentificationList& pep_ids,
  Size& best_id_index,
  bool& inconsistent_scores)
{
  const PeptideHit* best_hit = nullptr;
  const PeptideIdentification* ref_id = nullptr;   // fixes the score type and orientation
  std::function<bool(const PeptideHit&, const PeptideHit&)> better;
  best_id_index = 0;
  inconsistent_scores = false;

  for (Size i = 0; i < pep_ids.size(); ++i)
  {
    const auto& pid = pep_ids[i];
    if (pid.getHits().empty()) { continue; }       // hitless ids must not veto later ones

    if (ref_id == nullptr)
    {
      ref_id = &pid;
      better = PeptideIdentification::getScoreComparator(pid.isHigherScoreBetter());
    }
    else if (ref_id->getScoreType() != pid.getScoreType()
          || ref_id->isHigherScoreBetter() != pid.isHigherScoreBetter())
    {
      inconsistent_scores = true;                  // reported once per export by the caller
    }

    for (const auto& hit : pid.getHits())
    {
      if (best_hit == nullptr || better(hit, *best_hit))
      {
        best_hit = &hit;
        best_id_index = i;
      }
    }
  }
  return best_hit;
}

/// Build Parquet WriterProperties from the OpenMS config (shared by the one-shot
/// and streaming feature export paths).
std::shared_ptr<parquet::WriterProperties> makeFeatureWriterProperties(const ParquetWriteConfig& config)
{
  auto builder = parquet::WriterProperties::Builder();
  switch (config.compression)
  {
    case ParquetWriteConfig::Compression::NONE:   builder.compression(arrow::Compression::UNCOMPRESSED); break;
    case ParquetWriteConfig::Compression::SNAPPY: builder.compression(arrow::Compression::SNAPPY); break;
    case ParquetWriteConfig::Compression::GZIP:   builder.compression(arrow::Compression::GZIP); builder.compression_level(config.compression_level); break;
    case ParquetWriteConfig::Compression::LZ4:    builder.compression(arrow::Compression::LZ4); break;
    case ParquetWriteConfig::Compression::ZSTD:   builder.compression(arrow::Compression::ZSTD); builder.compression_level(config.compression_level); break;
  }
  builder.data_pagesize(config.data_page_size);
  if (config.write_statistics) { builder.enable_statistics(); }
  else                         { builder.disable_statistics(); }
  return builder.build();
}

/// Build a QPXFeatureSchema Arrow table for the half-open feature range
/// [range_begin, range_end) of @p cmap. @p id_lut holds the prebuilt per-map
/// protein-group/gene lookups. Each row is independent (no cross-row state), so
/// any contiguous sub-range yields a self-contained, schema-valid table; this is
/// what lets exportToParquetStreaming() flush one batch at a time. Shared by
/// exportToArrow() (whole map) and exportToParquetStreaming() (one batch per call).
std::shared_ptr<arrow::Table> buildFeatureTableRange(
  const ConsensusMap& cmap,
  const FeatureIdLookups& id_lut,
  size_t range_begin,
  size_t range_end)
{
  // -- Simple column builders --
  arrow::StringBuilder sequence_builder, peptidoform_builder, anchor_protein_builder;
  arrow::StringBuilder run_file_name_builder, id_run_file_name_builder;
  arrow::Int16Builder charge_builder, missed_cleavages_builder;
  arrow::BooleanBuilder is_decoy_builder, unique_builder;
  arrow::FloatBuilder calculated_mz_builder, observed_mz_builder, rt_builder;
  arrow::FloatBuilder mass_error_ppm_builder;
  arrow::FloatBuilder rt_start_builder, rt_stop_builder;
  arrow::FloatBuilder ion_mobility_builder, predicted_rt_builder;
  arrow::FloatBuilder im_start_builder, im_stop_builder;
  arrow::DoubleBuilder pep_builder, pg_qvalue_builder;

  // -- List<utf8> builders for gg_accessions, gg_names --
  auto gg_acc_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_accessions_builder(arrow::default_memory_pool(), gg_acc_vb);

  auto gg_names_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_names_builder(arrow::default_memory_pool(), gg_names_vb);

  // -- scan: list<int32> --
  auto scan_val_b = std::make_shared<arrow::Int32Builder>();
  arrow::ListBuilder scan_builder(arrow::default_memory_pool(), scan_val_b);

  // -- Modifications: list<struct{name, accession, positions: list<struct{position: int32, amino_acid, scores: list<struct>}>}> --
  // scores within position
  auto mod_score_name_b = std::make_shared<arrow::StringBuilder>();
  auto mod_score_value_b = std::make_shared<arrow::DoubleBuilder>();
  auto mod_score_hb_b = std::make_shared<arrow::BooleanBuilder>();
  auto mod_score_struct_type = arrow::struct_({
    arrow::field("score_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("score_value", arrow::float64(), /*nullable=*/false),
    arrow::field("higher_better", arrow::boolean())
  });
  auto mod_score_struct_b = std::make_shared<arrow::StructBuilder>(
    mod_score_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{mod_score_name_b, mod_score_value_b, mod_score_hb_b});
  auto mod_scores_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), mod_score_struct_b);

  auto pos_position_b = std::make_shared<arrow::Int32Builder>();
  auto pos_amino_acid_b = std::make_shared<arrow::StringBuilder>();
  auto position_struct_type = arrow::struct_({
    arrow::field("position", arrow::int32(), /*nullable=*/false),
    arrow::field("amino_acid", arrow::utf8()),
    arrow::field("scores", arrow::list(mod_score_struct_type))
  });
  auto position_struct_b = std::make_shared<arrow::StructBuilder>(
    position_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pos_position_b, pos_amino_acid_b, mod_scores_list_b});
  auto positions_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), position_struct_b);

  auto mod_name_b = std::make_shared<arrow::StringBuilder>();
  auto mod_acc_b = std::make_shared<arrow::StringBuilder>();
  auto modification_struct_type = arrow::struct_({
    arrow::field("name", arrow::utf8(), /*nullable=*/false),
    arrow::field("accession", arrow::utf8()),
    arrow::field("positions", arrow::list(position_struct_type), /*nullable=*/false)
  });
  auto modification_struct_b = std::make_shared<arrow::StructBuilder>(
    modification_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{mod_name_b, mod_acc_b, positions_list_b});
  arrow::ListBuilder modifications_builder(arrow::default_memory_pool(), modification_struct_b);

  // -- additional_scores: list<struct{score_name, score_value, higher_better}> --
  auto as_name_b = std::make_shared<arrow::StringBuilder>();
  auto as_value_b = std::make_shared<arrow::DoubleBuilder>();
  auto as_hb_b = std::make_shared<arrow::BooleanBuilder>();
  auto as_struct_type = arrow::struct_({
    arrow::field("score_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("score_value", arrow::float64(), /*nullable=*/false),
    arrow::field("higher_better", arrow::boolean())
  });
  auto as_struct_b = std::make_shared<arrow::StructBuilder>(
    as_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{as_name_b, as_value_b, as_hb_b});
  arrow::ListBuilder additional_scores_builder(arrow::default_memory_pool(), as_struct_b);

  // -- cv_params: list<struct{cv_name, cv_value}> --
  auto cv_name_b = std::make_shared<arrow::StringBuilder>();
  auto cv_value_b = std::make_shared<arrow::StringBuilder>();
  auto cv_struct_type = arrow::struct_({
    arrow::field("cv_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("cv_value", arrow::utf8(), /*nullable=*/false)
  });
  auto cv_struct_b = std::make_shared<arrow::StructBuilder>(
    cv_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{cv_name_b, cv_value_b});
  arrow::ListBuilder cv_params_builder(arrow::default_memory_pool(), cv_struct_b);

  // -- intensities: list<struct{label, intensity}> --
  auto intensity_label_b = std::make_shared<arrow::StringBuilder>();
  auto intensity_val_b = std::make_shared<arrow::FloatBuilder>();
  auto intensity_struct_type = arrow::struct_({
    arrow::field("label", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensity", arrow::float32(), /*nullable=*/false)
  });
  auto intensity_struct_b = std::make_shared<arrow::StructBuilder>(
    intensity_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{intensity_label_b, intensity_val_b});
  arrow::ListBuilder intensities_builder(arrow::default_memory_pool(), intensity_struct_b);

  // -- additional_intensities: list<struct{label, intensities: list<struct{intensity_name, intensity_value}>}> --
  auto ai_iname_b = std::make_shared<arrow::StringBuilder>();
  auto ai_ival_b = std::make_shared<arrow::FloatBuilder>();
  auto ai_entry_type = arrow::struct_({
    arrow::field("intensity_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensity_value", arrow::float32(), /*nullable=*/false)
  });
  auto ai_entry_b = std::make_shared<arrow::StructBuilder>(
    ai_entry_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ai_iname_b, ai_ival_b});
  auto ai_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), ai_entry_b);
  auto ai_label_b = std::make_shared<arrow::StringBuilder>();
  auto ai_outer_type = arrow::struct_({
    arrow::field("label", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensities", arrow::list(ai_entry_type), /*nullable=*/false)
  });
  auto ai_outer_b = std::make_shared<arrow::StructBuilder>(
    ai_outer_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ai_label_b, ai_list_b});
  arrow::ListBuilder additional_intensities_builder(arrow::default_memory_pool(), ai_outer_b);

  // -- pg_accessions: list<struct{accession, start, end, pre, post}> --
  auto pga_acc_b = std::make_shared<arrow::StringBuilder>();
  auto pga_start_b = std::make_shared<arrow::Int32Builder>();
  auto pga_end_b = std::make_shared<arrow::Int32Builder>();
  auto pga_pre_b = std::make_shared<arrow::StringBuilder>();
  auto pga_post_b = std::make_shared<arrow::StringBuilder>();
  auto pga_struct_type = arrow::struct_({
    arrow::field("accession", arrow::utf8(), /*nullable=*/false),
    arrow::field("start", arrow::int32()),
    arrow::field("end", arrow::int32()),
    arrow::field("pre", arrow::utf8()),
    arrow::field("post", arrow::utf8())
  });
  auto pga_struct_b = std::make_shared<arrow::StructBuilder>(
    pga_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pga_acc_b, pga_start_b, pga_end_b, pga_pre_b, pga_post_b});
  arrow::ListBuilder pg_accessions_builder(arrow::default_memory_pool(), pga_struct_b);

  // -- pg_positions: list<struct{protein_accession, start, end}> --
  auto pgp_acc_b = std::make_shared<arrow::StringBuilder>();
  auto pgp_start_b = std::make_shared<arrow::Int32Builder>();
  auto pgp_end_b = std::make_shared<arrow::Int32Builder>();
  auto pgp_struct_type = arrow::struct_({
    arrow::field("protein_accession", arrow::utf8(), /*nullable=*/false),
    arrow::field("start", arrow::int32(), /*nullable=*/false),
    arrow::field("end", arrow::int32(), /*nullable=*/false)
  });
  auto pgp_struct_b = std::make_shared<arrow::StructBuilder>(
    pgp_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pgp_acc_b, pgp_start_b, pgp_end_b});
  arrow::ListBuilder pg_positions_builder(arrow::default_memory_pool(), pgp_struct_b);

  arrow::Status status;

  // -- Melt: one output row per (ConsensusFeature, run) --
  //
  // QPX's feature view is long -- "a quantified peptidoform in a specific run file"
  // (docs/spec/feature.md) -- so a ConsensusFeature spanning N runs expands to N rows, each
  // carrying only that run's intensities. Emitting one wide row per feature instead made
  // run_file_name an arbitrary pick among the contributing runs and forced every intensity
  // into a single list distinguishable only by a non-conformant label.
  //
  // The coordinate rule differs by acquisition mode:
  //  - label-free: each FeatureHandle holds that run's own RT / m/z / charge / width, so the
  //    row takes its quantitative coordinates from the handle.
  //  - isobaric: the handles are reporter ions -- their m/z is the reporter mass and their RT
  //    the quantification scan -- while precursor m/z and RT live on the ConsensusFeature.
  //    Keep the consensus coordinates and use the handles only for channel intensities.
  const bool is_isobaric = id_lut.is_isobaric;

  struct RowSpec
  {
    size_t feat_idx;
    std::string run;                            ///< stemmed run_file_name
    std::vector<const FeatureHandle*> handles;  ///< this run's handles, in map-index order
  };

  std::vector<RowSpec> rows;
  rows.reserve(range_end - range_begin);
  Size dropped_handleless = 0;
  for (size_t feat_idx = range_begin; feat_idx < range_end; ++feat_idx)
  {
    // Group by the run's SOURCE path, not by its stemmed name. Two runs whose paths differ
    // only in directory (/a/run1.mzML, /b/run1.mzML) share a stem; grouping on the stem
    // would silently merge them into one row carrying both runs' intensities -- strictly
    // worse than the ambiguous-but-separate rows #9818 warns about. Channels of one isobaric
    // run share ColumnHeader::filename, so they still group together, which is the point.
    // std::map keeps the order deterministic, so output is stable for any thread count.
    std::map<std::string, std::vector<const FeatureHandle*>> by_source;
    for (const auto& fh : cmap[feat_idx].getFeatures())
    {
      auto it = id_lut.run_columns.find(fh.getMapIndex());
      if (it == id_lut.run_columns.end()) { continue; }
      by_source[it->second.source].push_back(&fh);
    }
    if (by_source.empty()) { ++dropped_handleless; continue; }
    for (auto& [source, handles] : by_source)
    {
      // The emitted key is the stem; the grouping key was the full path.
      rows.push_back(RowSpec{feat_idx,
                             id_lut.run_columns.at(handles.front()->getMapIndex()).run_name,
                             std::move(handles)});
    }
  }
  if (dropped_handleless > 0)
  {
    OPENMS_LOG_INFO << "ConsensusMapArrowExport: skipped " << dropped_handleless
                    << " consensus feature(s) with no feature handles (no quantification and "
                       "no source run to key a QPX feature row on)." << std::endl;
  }

  const Size num_features = rows.size();

  // Reserve capacity for simple column builders
  #define RESERVE_OR_FAIL(builder) \
    status = builder.Reserve(num_features); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: " #builder " Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  RESERVE_OR_FAIL(sequence_builder)
  RESERVE_OR_FAIL(peptidoform_builder)
  RESERVE_OR_FAIL(anchor_protein_builder)
  RESERVE_OR_FAIL(run_file_name_builder)
  RESERVE_OR_FAIL(id_run_file_name_builder)
  RESERVE_OR_FAIL(charge_builder)
  RESERVE_OR_FAIL(missed_cleavages_builder)
  RESERVE_OR_FAIL(is_decoy_builder)
  RESERVE_OR_FAIL(unique_builder)
  RESERVE_OR_FAIL(calculated_mz_builder)
  RESERVE_OR_FAIL(observed_mz_builder)
  RESERVE_OR_FAIL(mass_error_ppm_builder)
  RESERVE_OR_FAIL(rt_builder)
  RESERVE_OR_FAIL(rt_start_builder)
  RESERVE_OR_FAIL(rt_stop_builder)
  RESERVE_OR_FAIL(ion_mobility_builder)
  RESERVE_OR_FAIL(predicted_rt_builder)
  RESERVE_OR_FAIL(im_start_builder)
  RESERVE_OR_FAIL(im_stop_builder)
  RESERVE_OR_FAIL(pep_builder)
  RESERVE_OR_FAIL(pg_qvalue_builder)

  #undef RESERVE_OR_FAIL


  // Per-map protein-group q-value and gene-info lookups (prebuilt once, shared across batches)
  const auto& pg_qvalue_lookup = id_lut.pg_qvalue;
  const auto& gg_accessions_lookup = id_lut.gg_accessions;
  const auto& gg_names_lookup = id_lut.gg_names;

  // Metavalue keys excluded from additional_scores on PeptideHit
  static const std::unordered_set<std::string> excluded_hit_mvs = {
    "target_decoy", "predicted_RT", "predicted_rt"
  };

  IDScoreSwitcherAlgorithm idsa;

  // Lazily-filled MetaInfo index->name cache shared across all features in this range
  // (see cachedMetaName_). Per-call (hence per OpenMP partition) → thread-safe.
  std::vector<std::string> meta_name_cache;

  // Pre-resolve the fixed dedicated-column metavalue names to registry indices ONCE for this range.
  // String-keyed exists()/getValue() take the MetaInfoRegistry omp-critical lock via getIndex() on
  // every call; resolving once and using the index overloads (lock-free flat_map lookups) keeps the
  // parallel build from serializing on that lock. UInt(-1) = name never registered on any object =>
  // treated as absent (matches the sentinel check in the string-keyed exists()).
  auto& mreg = MetaInfoInterface::metaRegistry();
  const UInt idx_target_decoy       = mreg.getIndex("target_decoy");
  const UInt idx_predicted_RT       = mreg.getIndex("predicted_RT");
  const UInt idx_predicted_rt       = mreg.getIndex("predicted_rt");
  const UInt idx_ion_mobility       = mreg.getIndex("ion_mobility");
  const UInt idx_IM                 = mreg.getIndex("IM");
  const UInt idx_spectrum_reference = mreg.getIndex("spectrum_reference");
  const UInt idx_missed_cleavages   = mreg.getIndex("missed_cleavages");
  const UInt idx_start_ion_mobility = mreg.getIndex("start_ion_mobility");
  const UInt idx_stop_ion_mobility  = mreg.getIndex("stop_ion_mobility");
  const UInt idx_rt_start           = mreg.getIndex("rt_start");
  const UInt idx_rt_stop            = mreg.getIndex("rt_stop");
  // Presence check against a pre-resolved index; UInt(-1) is guarded because the index overload of
  // metaValueExists() does not special-case the sentinel the way the string overload does.
  auto metaHas = [](const MetaInfoInterface& m, UInt idx) -> bool
  { return idx != static_cast<UInt>(-1) && m.metaValueExists(idx); };

  bool warned_mixed_scores = false;
  for (const auto& row : rows)
  {
    const auto& cf = cmap[row.feat_idx];
    const auto& pep_ids = cf.getPeptideIdentifications();

    // The feature's identity: the best-scoring hit across ALL of its identifications, not
    // pep_ids[0].getHits()[0]. The positional pick discarded a feature's identity whenever its
    // *first* identification happened to carry no hits, even if a later one did.
    Size best_id_index = 0;
    bool inconsistent_scores = false;
    const PeptideHit* best_hit = selectFeatureIdentity(pep_ids, best_id_index, inconsistent_scores);
    const bool has_id = (best_hit != nullptr);
    if (inconsistent_scores && !warned_mixed_scores)
    {
      warned_mixed_scores = true;
      OPENMS_LOG_WARN << "ConsensusMapArrowExport: consensus feature carries identifications that "
                         "disagree on score type or orientation, so their scores are not "
                         "comparable. Harmonize them upstream (e.g. IDScoreSwitcher); the best "
                         "hit was selected using the first identification's orientation."
                      << std::endl;
    }

    // The row's quantitative coordinates, fixed here so that charge, observed_mz,
    // calculated_mz, mass_error_ppm, rt and the rt bounds are all mutually consistent.
    // Computing some against the consensus centroid and others against the run's handle
    // would report a ppm error measured against an m/z the row does not contain.
    //  - label-free: the run's own FeatureHandle carries that run's RT / m/z / charge / width.
    //  - isobaric: handles are reporter ions (m/z = reporter mass, RT = quantification scan),
    //    so precursor coordinates must come from the ConsensusFeature.
    const FeatureHandle* coord = (!is_isobaric && !row.handles.empty()) ? row.handles.front() : nullptr;
    const double row_mz     = coord ? coord->getMZ()     : cf.getMZ();
    const double row_rt     = coord ? coord->getRT()     : cf.getRT();
    const double row_width  = coord ? coord->getWidth()  : cf.getWidth();
    const int    row_charge = coord ? coord->getCharge() : cf.getCharge();

    // EVERY identification-level column must come from this one identification. Reading some
    // columns from pep_ids[0] and others from the winner would emit a row describing one
    // peptide's sequence next to another's scan number and ion mobility.
    const PeptideIdentification* best_pid =
      (has_id && best_id_index < pep_ids.size()) ? &pep_ids[best_id_index] : nullptr;

    // === sequence / peptidoform ===
    if (best_hit)
    {
      const auto& seq = best_hit->getSequence();
      (void)sequence_builder.Append(seq.toUnmodifiedString());
      auto pf = ProForma::fromAASequence(seq);
      (void)peptidoform_builder.Append(ProForma::toString(pf, ProForma::WriteMode::CANONICAL));

      // calculated_mz
      int charge = row_charge;
      if (charge > 0)
      {
        float calc_mz = static_cast<float>(seq.getMZ(charge));
        (void)calculated_mz_builder.Append(calc_mz);
        float obs_mz = static_cast<float>(row_mz);
        // mass_error_ppm = (observed - calculated) / calculated * 1e6
        if (calc_mz > 0)
        {
          (void)mass_error_ppm_builder.Append((obs_mz - calc_mz) / calc_mz * 1e6f);
        }
        else
        {
          (void)mass_error_ppm_builder.AppendNull();
        }
      }
      else
      {
        (void)calculated_mz_builder.Append(static_cast<float>(row_mz));
        (void)mass_error_ppm_builder.AppendNull();
      }

      // === modifications ===
      (void)modifications_builder.Append();
      if (seq.isModified())
      {
        // Collect modifications keyed by name to group positions
        // name -> (accession, [(1-based-position, amino_acid_char)])
        std::map<std::string, std::pair<std::string, std::vector<std::pair<int, std::string>>>> mod_map;

        // N-terminal
        if (seq.hasNTerminalModification())
        {
          const auto* mod = seq.getNTerminalModification();
          std::string name = seq.getNTerminalModificationName();
          std::string acc_str;
          if (mod && mod->getUniModRecordId() > 0)
          {
            acc_str = "UNIMOD:" + std::to_string(mod->getUniModRecordId());
          }
          mod_map[name] = {acc_str, {}};
          mod_map[name].second.push_back({0, ""});
        }

        // Residue modifications
        for (Size pos = 0; pos < seq.size(); ++pos)
        {
          const auto& residue = seq[pos];
          if (residue.isModified())
          {
            const auto* mod = residue.getModification();
            std::string name = residue.getModificationName();
            std::string acc_str;
            if (mod && mod->getUniModRecordId() > 0)
            {
              acc_str = "UNIMOD:" + std::to_string(mod->getUniModRecordId());
            }
            if (!mod_map.contains(name))
            {
              mod_map[name] = {acc_str, {}};
            }
            mod_map[name].second.push_back({static_cast<int>(pos + 1), std::string(residue.getOneLetterCode())});
          }
        }

        // C-terminal
        if (seq.hasCTerminalModification())
        {
          const auto* mod = seq.getCTerminalModification();
          std::string name = seq.getCTerminalModificationName();
          std::string acc_str;
          if (mod && mod->getUniModRecordId() > 0)
          {
            acc_str = "UNIMOD:" + std::to_string(mod->getUniModRecordId());
          }
          if (!mod_map.contains(name))
          {
            mod_map[name] = {acc_str, {}};
          }
          mod_map[name].second.push_back({static_cast<int>(seq.size() + 1), ""});
        }

        // Write modification structs
        for (const auto& [name, acc_positions] : mod_map)
        {
          (void)modification_struct_b->Append();
          (void)mod_name_b->Append(name);
          if (acc_positions.first.empty())
          {
            (void)mod_acc_b->AppendNull();
          }
          else
          {
            (void)mod_acc_b->Append(acc_positions.first);
          }
          (void)positions_list_b->Append();
          for (const auto& [pos_idx, aa] : acc_positions.second)
          {
            (void)position_struct_b->Append();
            (void)pos_position_b->Append(pos_idx);
            if (aa.empty())
            {
              (void)pos_amino_acid_b->AppendNull();
            }
            else
            {
              (void)pos_amino_acid_b->Append(aa);
            }
            (void)mod_scores_list_b->Append(); // empty scores list
          }
        }
      }
      // else: empty modifications list (Append() already called)
    }
    else
    {
      (void)sequence_builder.Append("");
      (void)peptidoform_builder.Append("");
      (void)calculated_mz_builder.Append(static_cast<float>(row_mz));
      (void)mass_error_ppm_builder.AppendNull();
      // empty modifications list
      (void)modifications_builder.Append();
    }

    // Quantitative coordinates: this run's handle for label-free, the consensus centroid for
    // isobaric (whose handles carry reporter-ion m/z and the quantification scan's RT).
    (void)charge_builder.Append(static_cast<int16_t>(row_charge));
    (void)observed_mz_builder.Append(static_cast<float>(row_mz));
    (void)rt_builder.Append(static_cast<float>(row_rt));

    // === PEP ===
    if (has_id)
    {
      auto result = idsa.findScoreType(*best_pid, IDScoreSwitcherAlgorithm::ScoreType::PEP);
      // Resolve the (dynamic) PEP score name to an index once per feature (lock-free use below).
      const UInt idx_pep = result.score_name.empty() ? static_cast<UInt>(-1) : mreg.getIndex(result.score_name);
      if (!result.score_name.empty() && best_hit)
      {
        if (result.is_main_score_type)
        {
          (void)pep_builder.Append(best_hit->getScore());
        }
        else if (metaHas(*best_hit, idx_pep))
        {
          (void)pep_builder.Append(static_cast<double>(best_hit->getMetaValue(idx_pep)));
        }
        else
        {
          (void)pep_builder.AppendNull();
        }
      }
      else
      {
        (void)pep_builder.AppendNull();
      }
    }
    else
    {
      (void)pep_builder.AppendNull();
    }

    // === is_decoy ===
    if (best_hit && metaHas(*best_hit, idx_target_decoy))
    {
      std::string td = best_hit->getMetaValue(idx_target_decoy).toString();
      (void)is_decoy_builder.Append(StringUtils::substr(td, 0, 5) == "decoy");
    }
    else
    {
      (void)is_decoy_builder.Append(false);
    }

    // === additional_scores from best hit ===
    (void)additional_scores_builder.Append();
    if (best_hit)
    {
      // Single index-based pass over the hit's metavalues (avoids getKeys()/getMetaValue(string),
      // which take the MetaInfoRegistry lock per key); names resolved via the shared cache. Same
      // iteration order (by index) and same filtering as the previous getKeys() loop.
      for (auto mv_it = best_hit->metaBegin(); mv_it != best_hit->metaEnd(); ++mv_it)
      {
        const std::string& key = cachedMetaName_(mv_it->first, meta_name_cache);
        if (excluded_hit_mvs.contains(key)) continue;
        const DataValue& val = mv_it->second;
        const auto vt = val.valueType();
        if ((vt == DataValue::INT_VALUE || vt == DataValue::DOUBLE_VALUE)
            && Scores::isKnownScoreType(key))
        {
          (void)as_struct_b->Append();
          (void)as_name_b->Append(key);
          (void)as_value_b->Append(static_cast<double>(val));
          try
          {
            auto st = IDScoreSwitcherAlgorithm::toScoreTypeEnum(key);
            (void)as_hb_b->Append(idsa.isScoreTypeHigherBetter(st));
          }
          catch (...)
          {
            (void)as_hb_b->AppendNull();
          }
        }
      }
    }

    // === predicted_rt ===
    if (best_hit && metaHas(*best_hit, idx_predicted_RT))
    {
      (void)predicted_rt_builder.Append(static_cast<float>(static_cast<double>(best_hit->getMetaValue(idx_predicted_RT))));
    }
    else if (best_hit && metaHas(*best_hit, idx_predicted_rt))
    {
      (void)predicted_rt_builder.Append(static_cast<float>(static_cast<double>(best_hit->getMetaValue(idx_predicted_rt))));
    }
    else
    {
      (void)predicted_rt_builder.AppendNull();
    }

    // === run_file_name ===
    // The run this melted row belongs to, already stemmed. QPX defines the column as the
    // spectrum file name without path or extension, and joins the psm, feature and pg views
    // on it, so all three must agree on the same spelling.
    // For isobaric input all channels of one file share ColumnHeader::filename (the channel
    // identity lives in ColumnHeader::label), so a run groups all of its channels.
    (void)run_file_name_builder.Append(row.run);

    // === scan (list<int32>) ===
    {
      std::string spec_ref;
      if (best_pid)   // NOT !pep_ids.empty(): all of them may be hitless
      {
        // Prefer the dedicated member, fall back to metavalue
        spec_ref = best_pid->getSpectrumReference();
        if (spec_ref.empty() && metaHas(*best_pid, idx_spectrum_reference))
        {
          spec_ref = best_pid->getMetaValue(idx_spectrum_reference).toString();
        }
      }

      (void)scan_builder.Append();
      if (!spec_ref.empty())
      {
        Int scan_num = extractScan(spec_ref);
        if (scan_num >= 0)
        {
          (void)scan_val_b->Append(scan_num);
        }
      }
    }

    // === cv_params (empty list for now) ===
    (void)cv_params_builder.AppendNull();

    // === ion_mobility, start/stop ===
    if (best_pid && metaHas(*best_pid, idx_ion_mobility))
    {
      (void)ion_mobility_builder.Append(static_cast<float>(static_cast<double>(best_pid->getMetaValue(idx_ion_mobility))));
    }
    else if (best_pid && metaHas(*best_pid, idx_IM))
    {
      (void)ion_mobility_builder.Append(static_cast<float>(static_cast<double>(best_pid->getMetaValue(idx_IM))));
    }
    else
    {
      (void)ion_mobility_builder.AppendNull();
    }

    if (metaHas(cf, idx_start_ion_mobility))
    {
      (void)im_start_builder.Append(static_cast<float>(static_cast<double>(cf.getMetaValue(idx_start_ion_mobility))));
    }
    else
    {
      (void)im_start_builder.AppendNull();
    }

    if (metaHas(cf, idx_stop_ion_mobility))
    {
      (void)im_stop_builder.Append(static_cast<float>(static_cast<double>(cf.getMetaValue(idx_stop_ion_mobility))));
    }
    else
    {
      (void)im_stop_builder.AppendNull();
    }

    // === missed_cleavages ===
    if (best_hit && metaHas(*best_hit, idx_missed_cleavages))
    {
      (void)missed_cleavages_builder.Append(static_cast<int16_t>(static_cast<int>(best_hit->getMetaValue(idx_missed_cleavages))));
    }
    else
    {
      (void)missed_cleavages_builder.AppendNull();
    }

    // === additional_intensities (empty list for now) ===
    (void)additional_intensities_builder.AppendNull();

    // === id_run_file_name ===
    // Spec: "run file containing best PSM". This may legitimately differ from run_file_name
    // — a match-between-runs feature is quantified in one run but identified in another, and
    // that divergence is exactly what this column records.
    // Resolved from the identification that actually holds the best hit, not pep_ids[0] --
    // a feature's identifications can come from different runs, which is the whole point.
    if (best_pid)
    {
      const std::string id_run =
        ArrowIOHelpers::qpxRunFileName(id_lut.run_mapper.getPrimaryMSRunPath(*best_pid));
      if (id_run.empty()) { (void)id_run_file_name_builder.AppendNull(); }
      else                { (void)id_run_file_name_builder.Append(id_run); }
    }
    else
    {
      (void)id_run_file_name_builder.AppendNull();
    }

    // === Protein accessions (now list<struct>) ===
    // Collect ALL peptide evidences (including repeated accessions with different positions)
    // and separately track unique accessions for anchor_protein/unique
    std::vector<std::string> protein_accs;
    std::unordered_set<std::string> seen_accs;
    struct EvidenceInfo { std::string acc; int start; int end; char pre; char post; };
    std::vector<EvidenceInfo> evidences;
    if (best_hit)
    {
      for (const auto& ev : best_hit->getPeptideEvidences())
      {
        const std::string& acc = ev.getProteinAccession();
        // Track unique accessions for anchor_protein / unique
        if (seen_accs.insert(acc).second)
        {
          protein_accs.push_back(acc);
        }
        // Emit every evidence (including repeated accessions at different positions)
        EvidenceInfo ei;
        ei.acc = acc;
        ei.start = ev.getStart();
        ei.end = ev.getEnd();
        ei.pre = ev.getAABefore();
        ei.post = ev.getAAAfter();
        evidences.push_back(ei);
      }
    }

    (void)pg_accessions_builder.Append();
    for (const auto& ei : evidences)
    {
      (void)pga_struct_b->Append();
      (void)pga_acc_b->Append(ei.acc);
      if (ei.start != PeptideEvidence::UNKNOWN_POSITION)
      {
        (void)pga_start_b->Append(ei.start);
      }
      else
      {
        (void)pga_start_b->AppendNull();
      }
      if (ei.end != PeptideEvidence::UNKNOWN_POSITION)
      {
        (void)pga_end_b->Append(ei.end);
      }
      else
      {
        (void)pga_end_b->AppendNull();
      }
      if (ei.pre != PeptideEvidence::UNKNOWN_AA && ei.pre != 0)
      {
        (void)pga_pre_b->Append(std::string(1, ei.pre));
      }
      else
      {
        (void)pga_pre_b->AppendNull();
      }
      if (ei.post != PeptideEvidence::UNKNOWN_AA && ei.post != 0)
      {
        (void)pga_post_b->Append(std::string(1, ei.post));
      }
      else
      {
        (void)pga_post_b->AppendNull();
      }
    }

    // anchor_protein — null, not "", when the feature has no protein mapping. bigbio/qpx#212
    // made this nullable for de novo workflows: "null when no protein mapping was performed".
    // An empty string would satisfy the column but is neither an accession nor an absence.
    if (protein_accs.empty())
    {
      (void)anchor_protein_builder.AppendNull();
    }
    else
    {
      (void)anchor_protein_builder.Append(protein_accs[0]);
    }

    // unique: true if exactly one unique protein accession
    (void)unique_builder.Append(protein_accs.size() == 1);

    // === pg_global_qvalue ===
    double min_qval = std::numeric_limits<double>::max();
    bool found_qval = false;
    for (const auto& acc : protein_accs)
    {
      auto it = pg_qvalue_lookup.find(acc);
      if (it != pg_qvalue_lookup.end())
      {
        if (it->second < min_qval)
        {
          min_qval = it->second;
        }
        found_qval = true;
      }
    }
    if (found_qval)
    {
      (void)pg_qvalue_builder.Append(min_qval);
    }
    else
    {
      (void)pg_qvalue_builder.AppendNull();
    }

    // === pg_positions (empty for now — no positional data available at consensus level) ===
    (void)pg_positions_builder.AppendNull();

    // === gg_accessions / gg_names ===
    std::unordered_set<std::string> seen_gg_acc, seen_gg_nam;
    (void)gg_accessions_builder.Append();
    (void)gg_names_builder.Append();
    for (const auto& acc : protein_accs)
    {
      auto it_a = gg_accessions_lookup.find(acc);
      if (it_a != gg_accessions_lookup.end())
      {
        for (const auto& ga : it_a->second)
        {
          if (seen_gg_acc.insert(ga).second)
          {
            (void)gg_acc_vb->Append(ga);
          }
        }
      }
      auto it_n = gg_names_lookup.find(acc);
      if (it_n != gg_names_lookup.end())
      {
        for (const auto& gn : it_n->second)
        {
          if (seen_gg_nam.insert(gn).second)
          {
            (void)gg_names_vb->Append(gn);
          }
        }
      }
    }

    // === rt_start / rt_stop ===
    // NOTE: the explicit rt_start/rt_stop metavalues are authored on the ConsensusFeature,
    // so a melted per-run row inherits the consensus window even when it reports the run
    // handle's RT. Where the two RTs diverge that window does not bracket the row's own RT.
    // Left as-is deliberately: overriding an authored value is a semantic change beyond this
    // PR, and existing coverage pins the current behaviour. The width fallback below does use
    // the row's own coordinates.
    if (metaHas(cf, idx_rt_start))
    {
      (void)rt_start_builder.Append(static_cast<float>(static_cast<double>(cf.getMetaValue(idx_rt_start))));
    }
    else if (row_width > 0)
    {
      (void)rt_start_builder.Append(static_cast<float>(row_rt - row_width / 2.0));
    }
    else
    {
      (void)rt_start_builder.AppendNull();
    }

    if (metaHas(cf, idx_rt_stop))
    {
      (void)rt_stop_builder.Append(static_cast<float>(static_cast<double>(cf.getMetaValue(idx_rt_stop))));
    }
    else if (row_width > 0)
    {
      (void)rt_stop_builder.Append(static_cast<float>(row_rt + row_width / 2.0));
    }
    else
    {
      (void)rt_stop_builder.AppendNull();
    }

    // === intensities: list<struct{label, intensity}> ===
    // Only this run's handles: one entry for label-free ("LFQ"), one per reporter channel for
    // isobaric ("TMT126", ...). The label is a join key against run.samples[].label, so it
    // must be a canonical channel token rather than the source file name.
    (void)intensities_builder.Append();
    for (const auto* fh : row.handles)
    {
      auto it = id_lut.run_columns.find(fh->getMapIndex());
      if (it == id_lut.run_columns.end()) { continue; }
      (void)intensity_struct_b->Append();
      (void)intensity_label_b->Append(it->second.label);
      (void)intensity_val_b->Append(static_cast<float>(fh->getIntensity()));
    }
  } // end melted row loop

  // Finalize all arrays
  #define FINISH_OR_FAIL(builder, arr) \
    status = builder.Finish(&arr); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: " #builder " Finish failed: " << status.ToString() << std::endl; return nullptr; }

  std::shared_ptr<arrow::Array> arr_sequence, arr_peptidoform, arr_modifications;
  std::shared_ptr<arrow::Array> arr_charge, arr_pep, arr_is_decoy;
  std::shared_ptr<arrow::Array> arr_calc_mz, arr_obs_mz, arr_mass_error_ppm;
  std::shared_ptr<arrow::Array> arr_additional_scores, arr_predicted_rt;
  std::shared_ptr<arrow::Array> arr_run_file, arr_cv_params, arr_scan;
  std::shared_ptr<arrow::Array> arr_rt, arr_ion_mobility, arr_missed_cleavages;
  std::shared_ptr<arrow::Array> arr_intensities, arr_additional_intensities;
  std::shared_ptr<arrow::Array> arr_pg_acc, arr_anchor, arr_unique, arr_pg_qval;
  std::shared_ptr<arrow::Array> arr_pg_positions;
  std::shared_ptr<arrow::Array> arr_im_start, arr_im_stop;
  std::shared_ptr<arrow::Array> arr_gg_acc, arr_gg_names;
  std::shared_ptr<arrow::Array> arr_id_run_file;
  std::shared_ptr<arrow::Array> arr_rt_start, arr_rt_stop;

  FINISH_OR_FAIL(sequence_builder, arr_sequence)
  FINISH_OR_FAIL(peptidoform_builder, arr_peptidoform)
  FINISH_OR_FAIL(modifications_builder, arr_modifications)
  FINISH_OR_FAIL(charge_builder, arr_charge)
  FINISH_OR_FAIL(pep_builder, arr_pep)
  FINISH_OR_FAIL(is_decoy_builder, arr_is_decoy)
  FINISH_OR_FAIL(calculated_mz_builder, arr_calc_mz)
  FINISH_OR_FAIL(observed_mz_builder, arr_obs_mz)
  FINISH_OR_FAIL(mass_error_ppm_builder, arr_mass_error_ppm)
  FINISH_OR_FAIL(additional_scores_builder, arr_additional_scores)
  FINISH_OR_FAIL(predicted_rt_builder, arr_predicted_rt)
  FINISH_OR_FAIL(run_file_name_builder, arr_run_file)
  FINISH_OR_FAIL(cv_params_builder, arr_cv_params)
  FINISH_OR_FAIL(scan_builder, arr_scan)
  FINISH_OR_FAIL(rt_builder, arr_rt)
  FINISH_OR_FAIL(ion_mobility_builder, arr_ion_mobility)
  FINISH_OR_FAIL(missed_cleavages_builder, arr_missed_cleavages)
  FINISH_OR_FAIL(intensities_builder, arr_intensities)
  FINISH_OR_FAIL(additional_intensities_builder, arr_additional_intensities)
  FINISH_OR_FAIL(pg_accessions_builder, arr_pg_acc)
  FINISH_OR_FAIL(anchor_protein_builder, arr_anchor)
  FINISH_OR_FAIL(unique_builder, arr_unique)
  FINISH_OR_FAIL(pg_qvalue_builder, arr_pg_qval)
  FINISH_OR_FAIL(pg_positions_builder, arr_pg_positions)
  FINISH_OR_FAIL(im_start_builder, arr_im_start)
  FINISH_OR_FAIL(im_stop_builder, arr_im_stop)
  FINISH_OR_FAIL(gg_accessions_builder, arr_gg_acc)
  FINISH_OR_FAIL(gg_names_builder, arr_gg_names)
  FINISH_OR_FAIL(id_run_file_name_builder, arr_id_run_file)
  FINISH_OR_FAIL(rt_start_builder, arr_rt_start)
  FINISH_OR_FAIL(rt_stop_builder, arr_rt_stop)

  #undef FINISH_OR_FAIL

  // Build schema from registry (QPX format)
  auto schema = QPXFeatureSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_sequence, arr_peptidoform, arr_modifications,
    arr_charge, arr_pep, arr_is_decoy,
    arr_calc_mz, arr_obs_mz, arr_mass_error_ppm,
    arr_additional_scores, arr_predicted_rt,
    arr_run_file, arr_cv_params, arr_scan,
    arr_rt, arr_ion_mobility, arr_missed_cleavages,
    arr_intensities, arr_additional_intensities,
    arr_pg_acc, arr_anchor, arr_unique, arr_pg_qval,
    arr_pg_positions, arr_im_start, arr_im_stop,
    arr_gg_acc, arr_gg_names, arr_id_run_file,
    arr_rt_start, arr_rt_stop,
  });

  // Validate table against registry schema (strict — write path must match exactly)
  auto validation = ArrowSchemaValidation::validate(table, QPXFeatureSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}

/// Build the canonical QPX "feature" file metadata for @p cmap. Each call mints a fresh
/// uuid/creation_date, so build it once per output file and reuse it for the writer schema
/// and every batch. Returns nullptr when QPX defines no token for the configured compression.
std::shared_ptr<const arrow::KeyValueMetadata> qpxFeatureMetadata(
  const ConsensusMap& cmap,
  const ParquetWriteConfig& config)
{
  // scan_format describes the `scan` column, which is extracted from the spectrum native
  // IDs of the features' peptide identifications (and the unassigned ones, which share the
  // same source runs).
  std::vector<std::string> refs;
  for (const auto& cf : cmap)
  {
    for (const auto& pid : cf.getPeptideIdentifications())
    {
      if (!pid.getSpectrumReference().empty()) { refs.push_back(pid.getSpectrumReference()); }
    }
  }
  for (const auto& pid : cmap.getUnassignedPeptideIdentifications())
  {
    if (!pid.getSpectrumReference().empty()) { refs.push_back(pid.getSpectrumReference()); }
  }

  std::map<std::string, std::string> extra;
  const std::string scan_format = ArrowIOHelpers::qpxScanFormat(refs);
  if (!scan_format.empty()) { extra["scan_format"] = scan_format; }
  return ArrowIOHelpers::qpxFileMetadata("feature_file", config, extra);
}

/// Refuse a merged run whose PSMs cannot be attributed to one of its files.
///
/// The feature view resolves id_run_file_name through IdentifierMSRunMapper, which falls back to
/// the run's FIRST file when 'id_merge_index' is missing -- so the column would name a run the
/// best PSM did not come from, silently. The psm and pg views already refuse this input; matching
/// them matters beyond consistency, because ProteomicsLFQ and ProteinQuantifier write the feature
/// view FIRST and pg third. Accepting here and refusing there leaves a wrong feature.parquet on
/// disk next to a half-written collection; refusing here fails before this view is written.
/// Note the pg ConsensusMap overload does NOT check merge indices -- only the psm view does --
/// so this is not redundant with it. QPXCollectionExport::requireExportable() is what makes the
/// whole collection fail before the FIRST file.
///
/// Must run as a preflight -- single-threaded, before any OpenMP region and before the output
/// file is opened -- for the same reason as requireUnambiguousIdentities().
void requireResolvableIdRunsWith(const ConsensusMap& cmap, const IdentifierMSRunMapper& mapper)
{
  for (Size i = 0; i < cmap.size(); ++i)
  {
    // Only the identification the row actually uses. Validating every attached one refused a
    // feature whose winning hit resolves perfectly, because a sibling identification -- hitless,
    // or simply not selected -- lacked an index that is never asked for.
    const auto& pep_ids = cmap[i].getPeptideIdentifications();
    Size best_id_index = 0;
    bool inconsistent_scores = false;
    if (selectFeatureIdentity(pep_ids, best_id_index, inconsistent_scores) == nullptr) { continue; }
    mapper.validateMergeIndex(pep_ids[best_id_index], i);
  }
}

} // anonymous namespace (feature range builder + shared per-map lookups)


void ConsensusMapArrowExport::requireUnambiguousIdentities(const ConsensusMap& cmap)
{
  Size divergent = 0;
  Size first_index = 0;
  for (Size i = 0; i < cmap.size(); ++i)
  {
    if (cmap[i].getAnnotationState() != BaseFeature::AnnotationState::FEATURE_ID_MULTIPLE_DIVERGENT) { continue; }
    if (divergent++ == 0) { first_index = i; }
  }
  if (divergent == 0) { return; }

  throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
    std::to_string(divergent) + " consensus feature(s), the first at index "
    + std::to_string(first_index) + ", carry peptide identifications that disagree on the top "
    "peptide. The QPX feature and protein-group views record one peptide per feature and cannot "
    "represent this. Resolve the conflicts first (IDConflictResolver, ConsensusID with "
    "algorithm 'best', or FileFilter with id:keep_best_score_id), or export the unresolved data "
    "as consensusparquet, which preserves every identification.");
}

void ConsensusMapArrowExport::requireResolvableIdRuns(const ConsensusMap& cmap)
{
  // Builds its own mapper so the check is usable without the rest of buildFeatureIdLookups().
  // The in-exporter paths pass the lookup's mapper instead, to avoid building it twice.
  std::set<std::string> seen_identifiers;
  for (const auto& prot_id : cmap.getProteinIdentifications())
  {
    if (!seen_identifiers.insert(prot_id.getIdentifier()).second)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "ConsensusMapArrowExport: several identification runs share the identifier '"
        + prot_id.getIdentifier() + "'. Their MS run paths cannot be told apart, so a feature "
        "would be attributed to the wrong run file.", prot_id.getIdentifier());
    }
  }
  IdentifierMSRunMapper mapper;
  try { mapper.create(cmap.getProteinIdentifications()); }
  catch (const Exception::InvalidValue&) { /* duplicate path lists; forward map is complete */ }
  requireResolvableIdRunsWith(cmap, mapper);
}

std::shared_ptr<arrow::Table> ConsensusMapArrowExport::exportToArrow(const ConsensusMap& cmap)
{
  requireUnambiguousIdentities(cmap);   // preflight: single-threaded, before any OpenMP region
  const auto id_lut = buildFeatureIdLookups(cmap);
  requireResolvableIdRunsWith(cmap, id_lut.run_mapper);
  // A channel QPX cannot name would be written into intensities[].label, a documented join
  // key. Refuse rather than emit a guessed or empty token (the pg exporter does the same).
  if (id_lut.has_unlabelable_channel) { return nullptr; }
  return buildFeatureTableRange(cmap, id_lut, 0, cmap.size());
}


bool ConsensusMapArrowExport::exportToParquet(
  const ConsensusMap& cmap,
  const std::string& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportToArrow(cmap);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to create Arrow table" << std::endl;
    return false;
  }

  auto meta = qpxFeatureMetadata(cmap, config);
  if (!meta) { return false; } // unsupported compression for QPX; already logged
  table = table->ReplaceSchemaMetadata(meta);

  QPXValueValidation value_validator(QPXValueValidation::View::FEATURE);
  const auto value_validation = value_validator.validate(table);
  if (!value_validation.valid)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: refusing invalid QPX feature table ("
                     << filename << "): " << value_validation.toString() << std::endl;
    return false;
  }

  // Open output file
  auto result = arrow::io::FileOutputStream::Open(filename);
  if (!result.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to open file: " << filename << std::endl;
    return false;
  }
  const auto& outfile = *result;

  auto writer_props = makeFeatureWriterProperties(config);
  auto arrow_props = parquet::ArrowWriterProperties::Builder().store_schema()->build();

  // Write
  auto status = parquet::arrow::WriteTable(
    *table, arrow::default_memory_pool(), outfile,
    config.row_group_size, writer_props, arrow_props);

  if (!status.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to write Parquet: "
                     << status.ToString() << std::endl;
    return false;
  }

  return true;
}


bool ConsensusMapArrowExport::exportToParquetStreaming(
  const ConsensusMap& cmap,
  const std::string& filename,
  size_t batch_size,
  const ParquetWriteConfig& config,
  int n_threads)
{
  if (batch_size == 0) { batch_size = 1000000; } // guard: a zero step would never advance

  // Resolve the number of OpenMP threads used to build each batch's partitions.
  // 1 = serial (default), 0 = all cores, N = fixed. Without OpenMP, always serial.
  int threads = n_threads;
#ifdef _OPENMP
  if (threads <= 0) { threads = omp_get_max_threads(); }
#endif
  if (threads < 1) { threads = 1; }

  // Prebuild the per-map protein-group/gene lookups once; reused for every batch.
  // Preflight before the writer is opened and before the OpenMP batch build: a throw from
  // inside that region would be captured by the worker's catch-all and flattened to a bool.
  requireUnambiguousIdentities(cmap);

  const auto id_lut = buildFeatureIdLookups(cmap);
  requireResolvableIdRunsWith(cmap, id_lut.run_mapper); // preflight, before the OpenMP region below
  if (id_lut.has_unlabelable_channel) { return false; } // already logged

  // Pre-warm lazily-initialized singletons/statics touched by the per-feature build, so first-touch
  // cannot race across the OpenMP workers below. Runs once, serially. (ProForma and the
  // SpectrumNativeIDParser scan path resolve through ModificationsDB/ResidueDB/AASequence, which are
  // warmed here — the same set the PSM streaming path pre-warms.)
  ModificationsDB::getInstance();
  ResidueDB::getInstance();
  (void)Scores::isKnownScoreType("q-value");
  (void)AASequence::fromString("PEPTIDEK").getMZ(2); // inits AASequence static residue path

  // One metadata identity (uuid/creation_date) for the whole file, reused for the writer
  // schema and every batch so their schemas are identical.
  auto meta = qpxFeatureMetadata(cmap, config);
  if (!meta) { return false; } // unsupported compression for QPX; already logged

  // The writer's schema must match each batch table's schema exactly. buildFeatureTableRange()
  // stamps QPXFeatureSchema::schema() on every batch, so open the writer with the same schema
  // plus the shared metadata (re-stamped per batch below).
  const auto schema = QPXFeatureSchema::schema()->WithMetadata(meta);

  auto file_result = arrow::io::FileOutputStream::Open(filename);
  if (!file_result.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to open " << filename << " for writing: "
                     << file_result.status().ToString() << std::endl;
    return false;
  }
  auto outfile = file_result.ValueOrDie();

  auto writer_props = makeFeatureWriterProperties(config);
  auto arrow_props = parquet::ArrowWriterProperties::Builder().store_schema()->build();

  auto writer_result = parquet::arrow::FileWriter::Open(
    *schema, arrow::default_memory_pool(), outfile, writer_props, arrow_props);
  if (!writer_result.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to open Parquet writer for " << filename << ": "
                     << writer_result.status().ToString() << std::endl;
    return false;
  }
  auto writer = std::move(writer_result).ValueOrDie();

  const int64_t rg = static_cast<int64_t>(config.row_group_size);
  QPXValueValidation value_validator(QPXValueValidation::View::FEATURE);

  // Build one feature batch [b, e): partition it into W contiguous sub-ranges, build each sub-table
  // in parallel (OpenMP), then write them in index order (serial — the FileWriter is not thread-safe).
  // Peak memory stays batch-bounded (the W sub-tables together hold one batch's rows). Row content and
  // order are identical to the serial path (contiguous, in-order), so output is deterministic for any
  // thread count.
  auto write_batch = [&](size_t b, size_t e) -> bool
  {
    const size_t rows = e - b;
    const int W = static_cast<int>(std::min<size_t>(static_cast<size_t>(threads), std::max<size_t>(rows, 1)));
    std::vector<std::shared_ptr<arrow::Table>> parts(W);
    std::vector<std::exception_ptr> errs(W);

    #pragma omp parallel for num_threads(W) schedule(static)
    for (int t = 0; t < W; ++t)
    {
      try
      {
        const size_t base = rows / W, rem = rows % W;
        const size_t sb = b + t * base + std::min<size_t>(static_cast<size_t>(t), rem);
        const size_t se = sb + base + (static_cast<size_t>(t) < rem ? 1 : 0);
        parts[t] = buildFeatureTableRange(cmap, id_lut, sb, se);
      }
      catch (...)
      {
        errs[t] = std::current_exception(); // never let an exception escape the OpenMP region
      }
    }

    // Surface any worker exception as a logged failure (convert to the bool error contract).
    for (int t = 0; t < W; ++t)
    {
      if (errs[t])
      {
        try { std::rethrow_exception(errs[t]); }
        catch (const std::exception& ex) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: feature batch build threw: " << ex.what() << std::endl; }
        catch (...)                      { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: feature batch build threw (unknown exception)" << std::endl; }
        return false;
      }
    }

    // Write partitions in index order (serial — the FileWriter is not thread-safe).
    for (int t = 0; t < W; ++t)
    {
      if (!parts[t])
      {
        OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to build feature batch partition " << t << std::endl;
        return false;
      }
      const auto value_validation = value_validator.validate(parts[t]);
      if (!value_validation.valid)
      {
        OPENMS_LOG_ERROR << "ConsensusMapArrowExport: refusing invalid QPX feature batch for "
                         << filename << ": " << value_validation.toString() << std::endl;
        return false;
      }
      auto st = writer->WriteTable(*parts[t]->ReplaceSchemaMetadata(meta), rg);
      if (!st.ok())
      {
        OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to write feature batch to " << filename << ": "
                         << st.ToString() << std::endl;
        return false;
      }
    }
    return true;
  };

  bool ok = true;
  const size_t N = cmap.size();
  if (N == 0)
  {
    ok = write_batch(0, 0); // valid 0-row file (schema only), matching the one-shot path
  }
  else
  {
    for (size_t b = 0; b < N && ok; b += batch_size)
    {
      ok = write_batch(b, std::min(b + batch_size, N));
    }
  }

  // Close the writer (flushes the footer) then the sink, regardless of write outcome.
  auto writer_close = writer->Close();
  if (!writer_close.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to close Parquet writer for " << filename << ": "
                     << writer_close.ToString() << std::endl;
    ok = false;
  }
  auto outfile_close = outfile->Close();
  if (!outfile_close.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to close output stream for " << filename << ": "
                     << outfile_close.ToString() << std::endl;
    ok = false;
  }
  return ok;
}

} // namespace OpenMS
