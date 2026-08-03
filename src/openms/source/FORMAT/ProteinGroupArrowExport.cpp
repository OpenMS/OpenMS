// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport_impl.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/IdentifierMSRunMapper.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <arrow/api.h>
#include <arrow/builder.h>

#include <algorithm>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{

using namespace OpenMS::Internal;

namespace // anonymous
{
  /// accession -> runs in which that accession has identification evidence.
  using AccessionRuns = std::unordered_map<std::string, std::set<std::string>>;

  /// Build the identifier -> spectra_data mapping used to resolve each PSM's origin run.
  ///
  /// Refuses duplicate identifiers: IdentifierMSRunMapper keys the forward map on the identifier
  /// and overwrites silently (IdentifierMSRunMapper.cpp), so two runs sharing one would leave
  /// only the last one's paths -- and every PSM of the first would then resolve to a file it
  /// never came from, stamping a foreign run name on its protein groups. Nothing downstream
  /// could detect that, so it has to be rejected here.
  ///
  /// A duplicate *path list* is different: create() throws for it, but populates the forward map
  /// -- the only direction used here -- completely first, so that one is downgraded to a warning.
  /// Deliberately narrow: any other exception means the forward map itself is incomplete, and
  /// silently continuing on it would turn malformed metadata into a plausible-looking export.
  IdentifierMSRunMapper buildRunMapper(const std::vector<ProteinIdentification>& prot_ids,
                                       const std::string& context)
  {
    std::set<std::string> seen;
    for (const auto& prot_id : prot_ids)
    {
      if (!seen.insert(prot_id.getIdentifier()).second)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          context + ": several identification runs share the identifier '" + prot_id.getIdentifier()
          + "'. Their MS run paths cannot be told apart, so protein groups would be attributed to "
            "the wrong run file.", prot_id.getIdentifier());
      }
    }

    IdentifierMSRunMapper mapper;
    try
    {
      mapper.create(prot_ids);
    }
    catch (const Exception::InvalidValue& e)
    {
      OPENMS_LOG_WARN << context << ": ambiguous MS run paths (" << e.getMessage()
                      << "); run attribution continues on the identifier -> path mapping."
                      << std::endl;
    }
    return mapper;
  }

  /// Record, for every accession a PSM supports, the run the PSM came from.
  ///
  /// These are identifications that contribute no quantified feature: in a ConsensusMap the
  /// unassigned ones, and in identification-only input every one of them. They establish that a
  /// group was seen in a run even when nothing was quantified there -- exactly the rows the melt
  /// exists to emit -- so they feed the run set without affecting any count.
  void indexEvidenceRuns(const IdentifierMSRunMapper& mapper,
                         const PeptideIdentificationList& pep_ids,
                         AccessionRuns& acc_runs)
  {
    for (const auto& pid : pep_ids)
    {
      if (pid.getHits().empty()) { continue; }
      const std::string run = ArrowIOHelpers::qpxRunFileName(mapper.getPrimaryMSRunPath(pid));
      if (run.empty()) { continue; }
      for (const auto& ev : pid.getHits()[0].getPeptideEvidences())
      {
        if (!ev.getProteinAccession().empty()) { acc_runs[ev.getProteinAccession()].insert(run); }
      }
    }
  }

  /// Runs where any member of @p group has identification evidence.
  std::set<std::string> evidenceRuns(const AccessionRuns& acc_runs,
                                     const ProteinIdentification::ProteinGroup& group)
  {
    std::set<std::string> runs;
    for (const auto& acc : group.accessions)
    {
      auto it = acc_runs.find(acc);
      if (it != acc_runs.end()) { runs.insert(it->second.begin(), it->second.end()); }
    }
    return runs;
  }

  /// Track emitted (anchor_protein, grouped_runs, label) tuples -- the QPX pg primary key -- and report
  /// duplicates once per export.
  ///
  /// Two ProteinIdentifications whose paths share a stem, a repeated accession group, or two
  /// identifiers mapping to one path list (which IdentifierMSRunMapper::create() overwrites
  /// silently, throwing only on duplicate path lists) all produce byte-different rows under one
  /// key. Emitting both is deliberate -- dropping a row loses data and no merge rule exists for
  /// two group probabilities -- but it must not be silent, because nothing downstream enforces
  /// the key either: qpxc validate reports duplicate PKs at severity "warning" and its is_valid
  /// ignores warnings.
  class PrimaryKeyWatch
  {
  public:
    void add(const std::string& anchor, const std::vector<std::string>& runs,
             const std::optional<std::string>& label)
    {
      // QPX stores fraction order but compares grouped_runs set-wise for primary-key identity.
      // Canonicalize only the watch key; never reorder the emitted list.
      std::vector<std::string> canonical_runs = runs;
      std::sort(canonical_runs.begin(), canonical_runs.end());
      const std::string key =
        ListUtils::concatenate(StringList(canonical_runs.begin(), canonical_runs.end()), "+");
      const auto pk = std::make_tuple(anchor, key, label);
      if (!seen_.insert(pk).second && examples_.size() < 3)
      {
        examples_.push_back(anchor + " / " + key + " / " + label.value_or("<null>"));
      }
      ++emitted_;
    }

    void report(const std::string& context) const
    {
      const size_t duplicates = emitted_ - seen_.size();
      if (duplicates == 0) { return; }
      OPENMS_LOG_WARN << context << ": " << duplicates << " of " << emitted_
                      << " rows repeat an (anchor_protein, grouped_runs, label) primary key, e.g. "
                      << ListUtils::concatenate(StringList(examples_.begin(), examples_.end()), "; ")
                      << ". The rows are kept, but consumers keying on the tuple will see collisions."
                      << std::endl;
    }

  private:
    std::set<std::tuple<std::string, std::string, std::optional<std::string>>> seen_;
    std::vector<std::string> examples_;
    size_t emitted_ = 0;
  };

} // anonymous namespace

std::shared_ptr<arrow::Table> ProteinGroupArrowExport::exportToArrow(const ConsensusMap& cmap)
{
  return exportToArrow(cmap, ExperimentalDesign::fromConsensusMap(cmap));
}

std::shared_ptr<arrow::Table> ProteinGroupArrowExport::exportToArrow(const ConsensusMap& cmap,
                                                                    const ExperimentalDesign& design)
{
  if (cmap.getProteinIdentifications().empty())
  {
    OPENMS_LOG_WARN << "ProteinGroupArrowExport: No protein identifications found" << std::endl;
    return nullptr;
  }

  // Same contract as the feature view: the pg view records one peptide per feature, so a
  // feature whose identifications disagree cannot be counted. Preflight, single-threaded.
  ConsensusMapArrowExport::requireUnambiguousIdentities(cmap);

  const auto& prot_id = cmap.getProteinIdentifications()[0];
  const auto& groups = prot_id.getIndistinguishableProteins();

  if (groups.empty())
  {
    OPENMS_LOG_WARN << "ProteinGroupArrowExport: No indistinguishable protein groups found" << std::endl;
    return nullptr;
  }

  // Build accession -> ProteinHit lookup
  std::unordered_map<std::string, const ProteinHit*> hit_lookup;
  for (const auto& hit : prot_id.getHits())
  {
    hit_lookup[hit.getAccession()] = &hit;
  }

  // -- Per-accession, per-run evidence index --
  //
  // Serves two purposes:
  //  * the run set for a group that has no quantification, so it can be melted onto the runs
  //    where it was actually identified instead of onto an empty grouped_runs (a QPX
  //    primary-key component that MUST NOT be null) or onto every run in the experiment;
  //  * the peptide/feature counts, which QPX defines per (group, run).
  //
  // The pg table is written after FDR filtering, so whatever PSMs survive in the map already
  // are the filtered set -- no extra score threshold is applied here.
  struct FeatureInfo
  {
    std::string sequence;
    std::string peptidoform;   ///< unmodified sequence + charge is not enough; see below
    int charge = 0;
  };
  std::vector<FeatureInfo> feature_info(cmap.size());
  // accession -> run -> contributing consensus feature indices
  std::unordered_map<std::string, std::map<std::string, std::set<Size>>> acc_run_features;
  // accession -> runs with identification evidence (features OR unassigned IDs)
  std::unordered_map<std::string, std::set<std::string>> acc_runs;

  {
    // Same warn-only policy as the PSM and feature exporters.
    std::vector<std::string> header_paths;
    for (const auto& [map_index, header] : cmap.getColumnHeaders()) { header_paths.push_back(header.filename); }
    (void)ArrowIOHelpers::qpxWarnOnRunNameCollisions("ProteinGroupArrowExport", header_paths);

    const IdentifierMSRunMapper run_mapper =
      buildRunMapper(cmap.getProteinIdentifications(), "ProteinGroupArrowExport");

    const auto& headers = cmap.getColumnHeaders();
    for (Size fi = 0; fi < cmap.size(); ++fi)
    {
      const auto& cf = cmap[fi];

      // The runs this feature was quantified in.
      std::set<std::string> runs;
      for (const auto& fh : cf.getFeatures())
      {
        auto it = headers.find(fh.getMapIndex());
        if (it != headers.end()) { runs.insert(ArrowIOHelpers::qpxRunFileName(it->second.filename)); }
      }
      if (runs.empty()) { continue; }

      for (const auto& pid : cf.getPeptideIdentifications())
      {
        if (pid.getHits().empty()) { continue; }
        const auto& hit = pid.getHits()[0];
        if (feature_info[fi].sequence.empty())
        {
          feature_info[fi].sequence = hit.getSequence().toUnmodifiedString();
          feature_info[fi].peptidoform = hit.getSequence().toString();
          feature_info[fi].charge = hit.getCharge();
        }
        for (const auto& ev : hit.getPeptideEvidences())
        {
          const std::string& acc = ev.getProteinAccession();
          if (acc.empty()) { continue; }
          for (const auto& run : runs)
          {
            acc_run_features[acc][run].insert(fi);
            acc_runs[acc].insert(run);
          }
        }
      }
    }

    indexEvidenceRuns(run_mapper, cmap.getUnassignedPeptideIdentifications(), acc_runs);
  }

  // -- Quantification units --
  //
  // Active QPX 1.1 emits one row per (group, fraction_group, label). grouped_runs is exactly the
  // set of run files in that experimental-design fraction group; sample identity does not join
  // fraction groups.
  // (run, label) cells the map actually has a column header for. Cell-level, not run-level: a
  // design channel with no header would otherwise enter the unit and be given an invented
  // intensity label by labelFor() below -- a documented join key that resolves to nothing.
  std::set<std::pair<std::string, unsigned>> header_cells;
  for (const auto& [map_index, header] : cmap.getColumnHeaders())
  {
    (void)map_index;
    header_cells.emplace(ArrowIOHelpers::qpxRunFileName(header.filename),
                         header.getLabelAsUInt(cmap.getExperimentType()));
  }

  QuantificationUnits units(design, header_cells);
  if (units.inconsistent())
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: the experimental design assigns several samples "
                        "to one (run file, label) combination, so a run does not belong to a "
                        "single QPX quantification unit and grouped_runs cannot be keyed."
                     << std::endl;
    return nullptr;
  }
  if (!units.invalidSample().empty())
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: " << units.invalidSample()
                     << ". The experimental design cannot index the quantified abundances."
                     << std::endl;
    return nullptr;
  }
  if (!units.runInSeveralGroups().empty())
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: " << units.runInSeveralGroups()
                     << ". QPX requires grouped_runs to be disjoint within one label."
                     << std::endl;
    return nullptr;
  }
  if (!units.missingCell().empty())
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: " << units.missingCell()
                     << ". The experimental design and the data disagree about the channel "
                        "layout -- an isobaric method with fewer channels than the design lists "
                        "is the usual cause. Writing that cell would put an invented label in "
                        "the row, which joins to nothing." << std::endl;
    return nullptr;
  }
  if (!units.ambiguousLabel().empty())
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: " << units.ambiguousLabel()
                     << ". A pg row carries one scalar label and intensity, so a label that means two "
                        "different samples across the fractions of one fraction group cannot be "
                        "written. Correct the design so each label has one sample within that "
                        "fraction group." << std::endl;
    return nullptr;
  }
  if (!units.raggedUnit().empty())
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: " << units.raggedUnit()
                     << ". grouped_runs names the files aggregated into the row's quantities, and "
                        "QPX resolves a sample from ANY of them, so a run that has no row for a "
                        "label must not appear beside that label's sample. A design with a missing "
                        "row is the usual cause: it bridges runs whose remaining samples were seen "
                        "in fewer files." << std::endl;
    return nullptr;
  }
  if (!units.droppedRuns().empty())
  {
    OPENMS_LOG_INFO << "ProteinGroupArrowExport: " << units.droppedRuns().size()
                    << " file(s) of the experimental design have no column header in this map and "
                       "are left out of grouped_runs, e.g. '" << *units.droppedRuns().begin()
                    << "'. They contributed no quantification." << std::endl;
  }

  // Quantified groups emit one row per (fraction_group, label). Identification-only groups emit
  // one row per evidence unit, which is bounded well enough by the same minimum reservation.
  Size quantities_per_group = 0;
  for (const auto& unit : units.units()) { quantities_per_group += unit.channels.size(); }
  const Size estimated_rows = groups.size() * std::max<Size>(quantities_per_group, 1);

  // -- Simple column builders --
  arrow::StringBuilder anchor_protein_builder;
  auto grouped_runs_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder grouped_runs_builder(arrow::default_memory_pool(), grouped_runs_vb);
  arrow::DoubleBuilder global_qvalue_builder, pg_qvalue_builder, gg_qvalue_builder;
  arrow::FloatBuilder sequence_coverage_builder, molecular_weight_builder;
  arrow::BooleanBuilder is_decoy_builder, contaminant_builder;

  // -- list<utf8> builders --
  auto pg_acc_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder pg_accessions_builder(arrow::default_memory_pool(), pg_acc_vb);

  auto pg_names_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder pg_names_builder(arrow::default_memory_pool(), pg_names_vb);

  auto gg_acc_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_accessions_builder(arrow::default_memory_pool(), gg_acc_vb);

  auto gg_names_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_names_builder(arrow::default_memory_pool(), gg_names_vb);

  // -- scalar quantity columns; nullable for identification-only evidence rows --
  arrow::StringBuilder label_builder;
  arrow::FloatBuilder intensity_builder;

  // -- additional_intensities: list<struct{label: utf8 not null, intensities: list<struct{intensity_name: utf8 not null, intensity_value: float32 not null}> not null}> --
  auto int_pair_type = arrow::struct_({
    arrow::field("intensity_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensity_value", arrow::float32(), /*nullable=*/false)
  });
  auto ip_name_b = std::make_shared<arrow::StringBuilder>();
  auto ip_value_b = std::make_shared<arrow::FloatBuilder>();
  auto ip_struct_b = std::make_shared<arrow::StructBuilder>(
    int_pair_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ip_name_b, ip_value_b});
  auto ip_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), ip_struct_b);

  auto add_int_struct_type = arrow::struct_({
    arrow::field("label", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensities", arrow::list(int_pair_type), /*nullable=*/false)
  });
  auto ai_label_b = std::make_shared<arrow::StringBuilder>();
  auto ai_struct_b = std::make_shared<arrow::StructBuilder>(
    add_int_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ai_label_b, ip_list_b});
  arrow::ListBuilder additional_intensities_builder(arrow::default_memory_pool(), ai_struct_b);

  // -- peptides: list<struct{protein_name: utf8 not null, peptide_count: int32 not null}> --
  auto pep_struct_type = arrow::struct_({
    arrow::field("protein_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("peptide_count", arrow::int32(), /*nullable=*/false)
  });
  auto pep_name_b = std::make_shared<arrow::StringBuilder>();
  auto pep_count_b = std::make_shared<arrow::Int32Builder>();
  auto pep_struct_b = std::make_shared<arrow::StructBuilder>(
    pep_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pep_name_b, pep_count_b});
  arrow::ListBuilder peptides_builder(arrow::default_memory_pool(), pep_struct_b);

  // -- peptide_counts: struct{unique_sequences: int32 not null, total_sequences: int32 not null} --
  auto pep_counts_type = arrow::struct_({
    arrow::field("unique_sequences", arrow::int32(), /*nullable=*/false),
    arrow::field("total_sequences", arrow::int32(), /*nullable=*/false)
  });
  auto pc_unique_b = std::make_shared<arrow::Int32Builder>();
  auto pc_total_b = std::make_shared<arrow::Int32Builder>();
  arrow::StructBuilder peptide_counts_builder(
    pep_counts_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pc_unique_b, pc_total_b});

  // -- feature_counts: struct{unique_features: int32 not null, total_features: int32 not null} --
  auto feat_counts_type = arrow::struct_({
    arrow::field("unique_features", arrow::int32(), /*nullable=*/false),
    arrow::field("total_features", arrow::int32(), /*nullable=*/false)
  });
  auto fc_unique_b = std::make_shared<arrow::Int32Builder>();
  auto fc_total_b = std::make_shared<arrow::Int32Builder>();
  arrow::StructBuilder feature_counts_builder(
    feat_counts_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{fc_unique_b, fc_total_b});

  // -- additional_scores: list<struct{score_name: utf8 not null, score_value: float64 not null, higher_better: bool}> --
  auto as_struct_type = arrow::struct_({
    arrow::field("score_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("score_value", arrow::float64(), /*nullable=*/false),
    arrow::field("higher_better", arrow::boolean())
  });
  auto as_name_b = std::make_shared<arrow::StringBuilder>();
  auto as_value_b = std::make_shared<arrow::DoubleBuilder>();
  auto as_hb_b = std::make_shared<arrow::BooleanBuilder>();
  auto as_struct_b = std::make_shared<arrow::StructBuilder>(
    as_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{as_name_b, as_value_b, as_hb_b});
  arrow::ListBuilder additional_scores_builder(arrow::default_memory_pool(), as_struct_b);

  // -- cv_params: list<struct{cv_name: utf8 not null, cv_value: utf8 not null}> --
  auto cv_struct_type = arrow::struct_({
    arrow::field("cv_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("cv_value", arrow::utf8(), /*nullable=*/false)
  });
  auto cv_name_b = std::make_shared<arrow::StringBuilder>();
  auto cv_value_b = std::make_shared<arrow::StringBuilder>();
  auto cv_struct_b = std::make_shared<arrow::StructBuilder>(
    cv_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{cv_name_b, cv_value_b});
  arrow::ListBuilder cv_params_builder(arrow::default_memory_pool(), cv_struct_b);

  // Reserve capacity
  arrow::Status status;
  status = anchor_protein_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: anchor_protein Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = grouped_runs_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: grouped_runs Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = global_qvalue_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: global_qvalue Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_qvalue_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: pg_qvalue Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = label_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: label Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = intensity_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: intensity Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_qvalue_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: gg_qvalue Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = sequence_coverage_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: sequence_coverage Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = molecular_weight_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: molecular_weight Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: is_decoy Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = contaminant_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: contaminant Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  // Helper lambda to emit one row for a protein group + fraction-group/label quantity.
  PrimaryKeyWatch pk_watch;
  auto emitRow = [&](const ProteinIdentification::ProteinGroup& group,
                     const std::vector<std::string>& runs,
                     const std::optional<std::string>& label,
                     const std::optional<float>& intensity)
  {
    const std::string& anchor = group.accessions.empty() ? "" : group.accessions[0];
    pk_watch.add(anchor, runs, label);

    // anchor_protein
    (void)anchor_protein_builder.Append(anchor);

    // grouped_runs -- the raw files aggregated into this quantity.
    // Deliberately no File::stemName() here: the values are already QPX run names, and
    // FileHandler::stripExtension is not idempotent - for an unrecognized extension it strips at
    // the last dot, so stemming twice would turn a legitimate stem like 'HeLa_1.2ug' into 'HeLa_1'.
    (void)grouped_runs_builder.Append();
    for (const auto& run : runs) { (void)grouped_runs_vb->Append(run); }

    // pg_accessions
    (void)pg_accessions_builder.Append();
    for (const auto& acc : group.accessions)
    {
      (void)pg_acc_vb->Append(acc);
    }

    // pg_names
    (void)pg_names_builder.Append();
    for (const auto& acc : group.accessions)
    {
      auto it = hit_lookup.find(acc);
      if (it != hit_lookup.end())
      {
        (void)pg_names_vb->Append(it->second->getDescription());
      }
      else
      {
        (void)pg_names_vb->Append("");
      }
    }

    // gg_accessions, gg_names (gene groups -- one entry per group member, keep parallel to
    // pg_accessions so consumers can zip positionally; empty string when the meta value is
    // missing). Appending only for annotated members would shorten these lists relative to
    // pg_accessions and shift every later gene onto the wrong protein.
    (void)gg_accessions_builder.Append();
    (void)gg_names_builder.Append();
    for (const auto& acc : group.accessions)
    {
      auto it = hit_lookup.find(acc);
      bool has_ga = (it != hit_lookup.end()) && it->second->metaValueExists("gene_accession");
      bool has_gn = (it != hit_lookup.end()) && it->second->metaValueExists("gene_name");
      (void)gg_acc_vb->Append(has_ga ? it->second->getMetaValue("gene_accession").toString() : "");
      (void)gg_names_vb->Append(has_gn ? it->second->getMetaValue("gene_name").toString() : "");
    }

    // gg_qvalue - not available
    (void)gg_qvalue_builder.AppendNull();

    // global_qvalue
    (void)global_qvalue_builder.Append(group.probability);

    // pg_qvalue - not per-run in OpenMS
    (void)pg_qvalue_builder.AppendNull();

    // label/intensity are both null for a row backed only by identification evidence.
    if (!label.has_value())
    {
      (void)label_builder.AppendNull();
      (void)intensity_builder.AppendNull();
    }
    else
    {
      (void)label_builder.Append(*label);
      (void)intensity_builder.Append(*intensity);
    }

    // additional_intensities - empty for quantified rows, null for identification-only rows.
    if (label.has_value()) { (void)additional_intensities_builder.Append(); }
    else                   { (void)additional_intensities_builder.AppendNull(); }

    // is_decoy
    auto anchor_it = hit_lookup.find(anchor);
    if (anchor_it != hit_lookup.end())
    {
      (void)is_decoy_builder.Append(anchor_it->second->isDecoy());
    }
    else
    {
      (void)is_decoy_builder.Append(false);
    }

    // contaminant - not tracked
    (void)contaminant_builder.AppendNull();

    // peptides / peptide_counts / feature_counts, all per (group, quantification unit) -- the
    // same unit the row's intensity describes, so the counts and the quantity cover one set of
    // runs rather than two.
    //
    // Semantics follow quantms (our interop target): "unique" means deduplicated, not
    // proteotypic. DIA-NN's converter reads the same fields as proteotypic counts, so these
    // numbers are only comparable within quantms-derived collections. See
    // qpx/converters/quantms/pg_adapter.py.
    {
      /// The (run, consensus feature) pairs @p acc contributes across the unit's runs.
      ///
      /// Keyed on the pair, not the feature alone, because that is the grain of the feature
      /// view: ConsensusMapArrowExport melts one ConsensusFeature into one row per run it was
      /// quantified in, so a feature spanning two fractions of this unit IS two feature rows.
      /// Counting bare feature indices here would make feature_counts.total_features disagree
      /// with the number of feature.parquet rows the same (group, unit) has.
      auto featuresOf = [&](const std::string& acc)
      {
        std::set<std::pair<std::string, Size>> features;
        auto ait = acc_run_features.find(acc);
        if (ait == acc_run_features.end()) { return features; }
        for (const auto& run : runs)
        {
          auto rit = ait->second.find(run);
          if (rit == ait->second.end()) { continue; }
          for (Size fi : rit->second) { features.emplace(run, fi); }
        }
        return features;
      };

      // The group's contributing feature rows in this unit, unioned over its members so a
      // feature shared by two members is counted once.
      std::set<std::pair<std::string, Size>> group_features;

      // peptides: one {protein_name, peptide_count} per group member, counting that member's
      // distinct peptide sequences in this unit.
      (void)peptides_builder.Append();
      int total_sequences = 0;   // sum of the per-member counts, i.e. what peptides[] adds up to
      for (const auto& acc : group.accessions)
      {
        const auto member_features = featuresOf(acc);
        group_features.insert(member_features.begin(), member_features.end());

        // Sequences deduplicate across the unit's runs: one peptide seen in two fractions is one
        // peptide. That is the difference between peptide_counts and feature_counts.
        std::set<std::string> seqs;
        for (const auto& [run, fi] : member_features)
        {
          (void)run;
          if (!feature_info[fi].sequence.empty()) { seqs.insert(feature_info[fi].sequence); }
        }
        (void)pep_struct_b->Append();
        (void)pep_name_b->Append(acc);
        (void)pep_count_b->Append(static_cast<int>(seqs.size()));
        total_sequences += static_cast<int>(seqs.size());
      }

      if (group_features.empty())
      {
        // Identification-only row: nothing was quantified in this unit, so there is nothing to
        // count. Null is the honest answer, distinct from a genuine zero.
        (void)peptide_counts_builder.AppendNull();
        (void)feature_counts_builder.AppendNull();
      }
      else
      {
        std::set<std::string> unique_sequences;
        std::set<std::pair<std::string, int>> unique_features;
        for (const auto& [run, fi] : group_features)
        {
          (void)run;
          if (!feature_info[fi].sequence.empty()) { unique_sequences.insert(feature_info[fi].sequence); }
          unique_features.emplace(feature_info[fi].peptidoform, feature_info[fi].charge);
        }
        // peptide_counts counts SEQUENCES: unique across the group, total as the sum of the
        // per-member counts emitted in peptides[] above (so the two agree, and shared peptides
        // are counted once per member). Distinguishing peptidoform and charge is what
        // feature_counts is for -- using group_features.size() here made total_sequences a
        // feature count, identical to total_features and contradicting peptides[] in the
        // same row. unique_features stays per peptidoform+charge, so it does NOT grow when a
        // feature spans two runs of the unit; only total_features does.
        (void)peptide_counts_builder.Append();
        (void)pc_unique_b->Append(static_cast<int>(unique_sequences.size()));
        (void)pc_total_b->Append(total_sequences);

        (void)feature_counts_builder.Append();
        (void)fc_unique_b->Append(static_cast<int>(unique_features.size()));
        (void)fc_total_b->Append(static_cast<int>(group_features.size()));
      }
    }

    // sequence_coverage
    if (anchor_it != hit_lookup.end() && anchor_it->second->getCoverage() != ProteinHit::COVERAGE_UNKNOWN)
    {
      (void)sequence_coverage_builder.Append(static_cast<float>(anchor_it->second->getCoverage()));
    }
    else
    {
      (void)sequence_coverage_builder.AppendNull();
    }

    // molecular_weight - not available
    (void)molecular_weight_builder.AppendNull();

    // additional_scores - anchor protein score
    (void)additional_scores_builder.Append();
    if (anchor_it != hit_lookup.end())
    {
      (void)as_struct_b->Append();
      (void)as_name_b->Append(prot_id.getScoreType());
      (void)as_value_b->Append(anchor_it->second->getScore());
      (void)as_hb_b->Append(prot_id.isHigherScoreBetter());
    }

    // cv_params - empty for now
    (void)cv_params_builder.Append();
  };

  // Build (run stem, 1-based channel) -> QPX intensity label from the column headers.
  //
  // A design label and a ColumnHeader channel are the same number: ColumnHeader::getLabelAsUInt()
  // returns channel_id + 1 and PeptideAndProteinQuant looks its abundances up under that same
  // 1-based number. Keying the lookup on the run as well as the channel keeps it correct for maps
  // whose runs use different methods, rather than assuming one global channel -> label mapping.
  std::map<std::pair<std::string, UInt>, std::string> channel_labels;
  for (const auto& [map_index, header] : cmap.getColumnHeaders())
  {
    const UInt channel = header.getLabelAsUInt(cmap.getExperimentType());
    const std::string channel_name = header.metaValueExists("channel_name")
                                   ? header.getMetaValue("channel_name").toString() : "";
    const std::string label = ArrowIOHelpers::qpxIntensityLabel(header.label, channel_name);
    if (label.empty()) { return nullptr; } // unidentifiable channel; already logged
    channel_labels[{ArrowIOHelpers::qpxRunFileName(header.filename), channel}] = label;
  }

  // Resolve a (run, channel) to its QPX label. This cannot miss: channel_labels above and the
  // header_cells that gate a design row into a unit are built from the same column headers under
  // the same key, {qpxRunFileName(filename), getLabelAsUInt(experimentType)}, so every channel a
  // unit carries has an entry here.
  //
  // It used to invent one -- "LFQ" for channel 1, silently, or the bare channel number, warned
  // once per export so every further channel was silent. The scalar label is a documented join
  // key against run.samples[].label, so an invented token resolves to nothing: a TMT10-
  // shaped design against a 6-plex map wrote a pg view carrying "7".."10" and exited 0. The
  // feature view refuses an unnameable channel (ConsensusMapArrowExport), and so does this one
  // now. Keeping the branch as a refusal rather than deleting it makes the coupling between the
  // two key sets a checked invariant instead of a comment: they live in separate places and
  // nothing but this check would notice if one drifted.
  auto labelFor = [&](const std::string& run, Int channel) -> std::string
  {
    auto it = channel_labels.find({run, static_cast<UInt>(channel)});
    if (it != channel_labels.end()) { return it->second; }
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: no column header describes channel " << channel
                     << " of run '" << run << "', so its label cannot be named and "
                        "would not join against run.parquet." << std::endl;
    return "";
  };
  if (!units.resolveLabels(labelFor, "ProteinGroupArrowExport")) { return nullptr; }

  std::set<std::pair<UInt, UInt>> expected_quantity_keys;
  for (const auto& unit : units.units())
  {
    for (const auto& channel : unit.channels)
    {
      expected_quantity_keys.emplace(unit.fraction_group, channel.label);
    }
  }

  // Iterate protein groups and emit rows
  Size unattributable_groups = 0;
  for (const auto& group : groups)
  {
    // The legacy sample abundance array identifies a quantified group and remains part of the
    // mzTab-facing annotation. QPX must read the new fraction-group/label arrays: a sample may
    // occur in multiple fraction groups, so its experiment-wide abundance cannot be split back
    // into the independent quantities after protein aggregation.
    const auto* sample_abundances = sampleAbundances(group);
    const FractionGroupAbundanceData quantities = fractionGroupAbundances(group);
    if (!quantities.valid)
    {
      OPENMS_LOG_ERROR << "ProteinGroupArrowExport: protein group '"
                       << (group.accessions.empty() ? "" : group.accessions[0]) << "' has invalid "
                       << "fraction-group abundance annotations: " << quantities.error << "."
                       << std::endl;
      return nullptr;
    }
    if (sample_abundances != nullptr && sample_abundances->size() != design.getNumberOfSamples())
    {
      OPENMS_LOG_ERROR << "ProteinGroupArrowExport: protein group '"
                       << (group.accessions.empty() ? "" : group.accessions[0]) << "' carries "
                       << sample_abundances->size() << " sample abundances but the experimental design "
                       << "describes " << design.getNumberOfSamples() << " samples. The design "
                       << "does not belong to this quantification, so the quantities cannot be "
                       << "attributed to their runs." << std::endl;
      return nullptr;
    }

    const bool quantified = sample_abundances != nullptr || quantities.present;
    if (quantified)
    {
      if (!quantities.present)
      {
        OPENMS_LOG_ERROR << "ProteinGroupArrowExport: protein group '"
                         << (group.accessions.empty() ? "" : group.accessions[0])
                         << "' carries sample abundances but no fraction-group/label abundances. "
                            "It was quantified before the QPX 1.1 unit-grain annotation was "
                            "created and cannot be exported without changing its quantities."
                         << std::endl;
        return nullptr;
      }
      if (units.empty())
      {
        OPENMS_LOG_ERROR << "ProteinGroupArrowExport: protein group '"
                         << (group.accessions.empty() ? "" : group.accessions[0])
                         << "' carries fraction-group abundances, but the experimental design "
                            "yields no matching quantification unit -- it lists no MS files this "
                            "map has a column header for." << std::endl;
        return nullptr;
      }
      if (quantities.values.size() != expected_quantity_keys.size())
      {
        OPENMS_LOG_ERROR << "ProteinGroupArrowExport: protein group '"
                         << (group.accessions.empty() ? "" : group.accessions[0]) << "' carries "
                         << quantities.values.size() << " fraction-group/label quantities but the "
                            "experimental design requires " << expected_quantity_keys.size()
                         << ". The annotations and design do not describe the same quantification."
                         << std::endl;
        return nullptr;
      }

      // One scalar row per label in each fraction group, matching the active QPX 1.1 primary key.
      for (const auto& unit : units.units())
      {
        for (const auto& channel : unit.channels)
        {
          const auto value = quantities.values.find({unit.fraction_group, channel.label});
          if (value == quantities.values.end())
          {
            OPENMS_LOG_ERROR << "ProteinGroupArrowExport: protein group '"
                             << (group.accessions.empty() ? "" : group.accessions[0])
                             << "' has no quantity for fraction group " << unit.fraction_group
                             << ", label " << channel.label << "." << std::endl;
            return nullptr;
          }
          emitRow(group, unit.runs, channel.qpx_label, value->second);
        }
      }
    }
    else
    {
      // No quantification for this group. grouped_runs is a QPX primary-key component that MUST
      // NOT be empty, and quantms silently drops rows that violate that, so a run-less row is
      // not an option. Melt onto the units where the group actually has identification evidence,
      // with null label/intensity. Fanning out across every unit in the experiment instead would
      // assert observations that were never made.

      // Evidence is per run; the row is per unit, so the runs of one unit collapse into a single
      // row rather than repeating the group once per fraction. A run the design does not list is
      // its own unit, which is also what an unfractionated design would have produced for it.
      std::set<std::vector<std::string>> evidence_units;   // ordered, so the rows are stable
      for (const auto& run : evidenceRuns(acc_runs, group))
      {
        const size_t index = units.indexOf(run);
        evidence_units.insert(index < units.units().size() ? units.units()[index].runs
                                                           : std::vector<std::string>{run});
      }
      if (evidence_units.empty())
      {
        ++unattributable_groups;
        continue;
      }
      for (const auto& runs : evidence_units) { emitRow(group, runs, std::nullopt, std::nullopt); }
    }
  }
  pk_watch.report("ProteinGroupArrowExport");
  if (unattributable_groups > 0)
  {
    OPENMS_LOG_INFO << "ProteinGroupArrowExport: skipped " << unattributable_groups
                    << " protein group(s) with neither quantification nor identification "
                       "evidence in any run (no grouped_runs to key a QPX row on)." << std::endl;
  }

  // Finalize all arrays
  std::shared_ptr<arrow::Array> arr_pg_acc, arr_pg_names, arr_gg_acc, arr_gg_names;
  std::shared_ptr<arrow::Array> arr_gg_qval, arr_anchor, arr_grouped_runs;
  std::shared_ptr<arrow::Array> arr_global_qval, arr_pg_qval;
  std::shared_ptr<arrow::Array> arr_label, arr_intensity, arr_add_intensities;
  std::shared_ptr<arrow::Array> arr_is_decoy, arr_contaminant;
  std::shared_ptr<arrow::Array> arr_peptides, arr_pep_counts, arr_feat_counts;
  std::shared_ptr<arrow::Array> arr_seq_cov, arr_mol_weight;
  std::shared_ptr<arrow::Array> arr_add_scores, arr_cv_params;

  status = pg_accessions_builder.Finish(&arr_pg_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: pg_accessions Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_names_builder.Finish(&arr_pg_names);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: pg_names Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_accessions_builder.Finish(&arr_gg_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: gg_accessions Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_names_builder.Finish(&arr_gg_names);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: gg_names Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_qvalue_builder.Finish(&arr_gg_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: gg_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = anchor_protein_builder.Finish(&arr_anchor);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: anchor_protein Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = grouped_runs_builder.Finish(&arr_grouped_runs);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: grouped_runs Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = global_qvalue_builder.Finish(&arr_global_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: global_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_qvalue_builder.Finish(&arr_pg_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: pg_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = label_builder.Finish(&arr_label);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: label Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = intensity_builder.Finish(&arr_intensity);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: intensity Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_intensities_builder.Finish(&arr_add_intensities);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: additional_intensities Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Finish(&arr_is_decoy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: is_decoy Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = contaminant_builder.Finish(&arr_contaminant);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: contaminant Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptides_builder.Finish(&arr_peptides);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: peptides Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptide_counts_builder.Finish(&arr_pep_counts);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: peptide_counts Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = feature_counts_builder.Finish(&arr_feat_counts);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: feature_counts Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = sequence_coverage_builder.Finish(&arr_seq_cov);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: sequence_coverage Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = molecular_weight_builder.Finish(&arr_mol_weight);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: molecular_weight Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_scores_builder.Finish(&arr_add_scores);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: additional_scores Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = cv_params_builder.Finish(&arr_cv_params);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: cv_params Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // One schema definition for both overloads: this one used to restate QPXPgSchema field by
  // field, which is what let the two drift apart. exportToParquet() validates against the
  // registry anyway, and Table::Make + Validate() below catches any type the builders emit that
  // the registry does not declare.
  auto table = arrow::Table::Make(QPXPgSchema::schema(), {
    arr_pg_acc, arr_pg_names, arr_gg_acc, arr_gg_names,
    arr_gg_qval, arr_anchor, arr_grouped_runs,
    arr_global_qval, arr_pg_qval,
    arr_label, arr_intensity, arr_add_intensities,
    arr_is_decoy, arr_contaminant,
    arr_peptides, arr_pep_counts, arr_feat_counts,
    arr_seq_cov, arr_mol_weight,
    arr_add_scores, arr_cv_params,
  });

  // Appends discard their arrow::Status throughout this file, so a failed one (e.g. under memory
  // pressure) leaves that column shorter than the rest while every Finish() still succeeds
  // independently. Validate() catches the resulting length mismatch, so the contract of returning
  // nullptr on failure holds instead of handing back a malformed table.
  if (auto valid = table->Validate(); !valid.ok())
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: built an inconsistent table: "
                     << valid.ToString() << std::endl;
    return nullptr;
  }
  return table;
}


namespace
{
  /// Attach the canonical QPX "pg" file metadata to a table's schema.
  /// Returns nullptr if QPX does not define a token for the configured compression.
  std::shared_ptr<arrow::Table> attachQPXPgMetadata(
    const std::shared_ptr<arrow::Table>& table,
    const ParquetWriteConfig& config)
  {
    // The pg view has no scan column, so no scan_format key.
    auto metadata = ArrowIOHelpers::qpxFileMetadata("pg_file", config);
    if (!metadata) { return nullptr; }
    return table->ReplaceSchemaMetadata(metadata);
  }
}

bool ProteinGroupArrowExport::exportToParquet(
  const ConsensusMap& cmap,
  const std::string& filename,
  const ParquetWriteConfig& config)
{
  return exportToParquet(cmap, ExperimentalDesign::fromConsensusMap(cmap), filename, config);
}

bool ProteinGroupArrowExport::exportToParquet(
  const ConsensusMap& cmap,
  const ExperimentalDesign& design,
  const std::string& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportToArrow(cmap, design);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: Failed to create Arrow table" << std::endl;
    return false;
  }
  return exportToParquet(table, filename, config);
}

bool ProteinGroupArrowExport::exportToParquet(
  const std::shared_ptr<arrow::Table>& table,
  const std::string& filename,
  const ParquetWriteConfig& config)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: null table passed to exportToParquet (" << filename << ")" << std::endl;
    return false;
  }

  // Guard: this overload stamps file_type="pg_file", so the caller must pass a
  // QPXPgSchema table.
  auto validation = ArrowSchemaValidation::validate(table, QPXPgSchema::schema(), ArrowSchemaValidation::Mode::Strict);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: table schema does not match QPXPgSchema ("
                     << filename << "): " << validation.toString() << std::endl;
    return false;
  }

  auto stamped = attachQPXPgMetadata(table, config);
  if (!stamped) { return false; } // unsupported compression for QPX; already logged
  return ArrowIOHelpers::writeTableToParquet(stamped, filename, config);
}


std::shared_ptr<arrow::Table> ProteinGroupArrowExport::exportToArrow(
  const std::vector<ProteinIdentification>& protein_identifications,
  const PeptideIdentificationList& peptide_identifications)
{
  auto target_schema = QPXPgSchema::schema();

  if (protein_identifications.empty())
  {
    OPENMS_LOG_WARN << "ProteinGroupArrowExport (id-only): No protein identifications found" << std::endl;
    return arrow::Table::MakeEmpty(target_schema).ValueOrDie();
  }

  // Compute total estimated rows across all protein identifications
  Size estimated_rows = 0;
  for (const auto& pid : protein_identifications)
  {
    estimated_rows += pid.getIndistinguishableProteins().size();
  }

  // -- Simple column builders --
  //
  // Identification input carries no experimental design, so nothing aggregates several runs into
  // one quantity: every grouped_runs list holds exactly the one run that keyed the row -- the
  // single-element form the QPX 1.1 schema documents for unfractionated input.
  arrow::StringBuilder anchor_protein_builder;
  auto grouped_runs_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder grouped_runs_builder(arrow::default_memory_pool(), grouped_runs_vb,
                                          QPXPgSchema::groupedRunsType());
  arrow::DoubleBuilder global_qvalue_builder, pg_qvalue_builder, gg_qvalue_builder;
  arrow::FloatBuilder sequence_coverage_builder, molecular_weight_builder;
  arrow::BooleanBuilder is_decoy_builder, contaminant_builder;

  // -- list<utf8> builders --
  auto pg_acc_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder pg_accessions_builder(arrow::default_memory_pool(), pg_acc_vb);

  auto pg_names_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder pg_names_builder(arrow::default_memory_pool(), pg_names_vb);

  auto gg_acc_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_accessions_builder(arrow::default_memory_pool(), gg_acc_vb);

  auto gg_names_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_names_builder(arrow::default_memory_pool(), gg_names_vb);

  // -- scalar label/intensity: null (not available for id-only) --
  arrow::StringBuilder label_builder;
  arrow::FloatBuilder intensity_builder;

  // -- additional_intensities: null (not available for id-only) --
  auto add_int_value_type = std::dynamic_pointer_cast<arrow::ListType>(QPXPgSchema::additionalIntensitiesType())->value_type();
  // inner list type for the "intensities" field of the struct
  auto int_pair_type = arrow::struct_({
    arrow::field("intensity_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensity_value", arrow::float32(), /*nullable=*/false)
  });
  auto ip_name_b = std::make_shared<arrow::StringBuilder>();
  auto ip_value_b = std::make_shared<arrow::FloatBuilder>();
  auto ip_struct_b = std::make_shared<arrow::StructBuilder>(
    int_pair_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ip_name_b, ip_value_b});
  auto ip_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), ip_struct_b, arrow::list(int_pair_type));
  auto ai_label_b = std::make_shared<arrow::StringBuilder>();
  auto ai_struct_b = std::make_shared<arrow::StructBuilder>(
    add_int_value_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ai_label_b, ip_list_b});
  arrow::ListBuilder additional_intensities_builder(arrow::default_memory_pool(), ai_struct_b, QPXPgSchema::additionalIntensitiesType());

  // -- peptides: list<struct{protein_name: utf8, peptide_count: int32}> --
  auto pep_value_type = std::dynamic_pointer_cast<arrow::ListType>(QPXPgSchema::peptidesType())->value_type();
  auto pep_name_b = std::make_shared<arrow::StringBuilder>();
  auto pep_count_b = std::make_shared<arrow::Int32Builder>();
  auto pep_struct_b = std::make_shared<arrow::StructBuilder>(
    pep_value_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pep_name_b, pep_count_b});
  arrow::ListBuilder peptides_builder(arrow::default_memory_pool(), pep_struct_b, QPXPgSchema::peptidesType());

  // -- peptide_counts: struct{unique_sequences, total_sequences} --
  auto pep_counts_type = QPXPgSchema::peptideCountsType();
  auto pc_unique_b = std::make_shared<arrow::Int32Builder>();
  auto pc_total_b = std::make_shared<arrow::Int32Builder>();
  arrow::StructBuilder peptide_counts_builder(
    pep_counts_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pc_unique_b, pc_total_b});

  // -- feature_counts: null (not available for id-only) --
  auto feat_counts_type = QPXPgSchema::featureCountsType();
  auto fc_unique_b = std::make_shared<arrow::Int32Builder>();
  auto fc_total_b = std::make_shared<arrow::Int32Builder>();
  arrow::StructBuilder feature_counts_builder(
    feat_counts_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{fc_unique_b, fc_total_b});

  // -- additional_scores: list<struct{score_name, score_value, higher_better}> --
  auto as_type = QPXPgSchema::additionalScoresType();
  auto as_value_type = std::dynamic_pointer_cast<arrow::ListType>(as_type)->value_type();
  auto as_name_b = std::make_shared<arrow::StringBuilder>();
  auto as_value_b = std::make_shared<arrow::DoubleBuilder>();
  auto as_hb_b = std::make_shared<arrow::BooleanBuilder>();
  auto as_struct_b = std::make_shared<arrow::StructBuilder>(
    as_value_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{as_name_b, as_value_b, as_hb_b});
  arrow::ListBuilder additional_scores_builder(arrow::default_memory_pool(), as_struct_b, as_type);

  // -- cv_params: list<struct{cv_name, cv_value}> --
  auto cv_type = QPXPgSchema::cvParamsType();
  auto cv_value_type = std::dynamic_pointer_cast<arrow::ListType>(cv_type)->value_type();
  auto cv_name_b = std::make_shared<arrow::StringBuilder>();
  auto cv_value_b = std::make_shared<arrow::StringBuilder>();
  auto cv_struct_b = std::make_shared<arrow::StructBuilder>(
    cv_value_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{cv_name_b, cv_value_b});
  arrow::ListBuilder cv_params_builder(arrow::default_memory_pool(), cv_struct_b, cv_type);

  // Reserve capacity
  arrow::Status status;
  status = anchor_protein_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): anchor_protein Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = grouped_runs_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): grouped_runs Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = global_qvalue_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): global_qvalue Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = label_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): label Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = intensity_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): intensity Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): is_decoy Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  // Resolve each PSM's origin run, so a group can be melted onto the runs where it was actually
  // identified.
  const IdentifierMSRunMapper mapper =
    buildRunMapper(protein_identifications, "ProteinGroupArrowExport (id-only)");
  {
    // Warn-only, as in the psm and feature views: two paths sharing a stem are a legitimate
    // layout that the spec's representation simply cannot tell apart.
    std::vector<std::string> all_paths;
    for (const auto& identifier : mapper.getIdentifiers())
    {
      const StringList& paths = mapper.getMSRunPaths(identifier);
      all_paths.insert(all_paths.end(), paths.begin(), paths.end());
    }
    (void)ArrowIOHelpers::qpxWarnOnRunNameCollisions("ProteinGroupArrowExport (id-only)", all_paths);
  }

  // Bucket the PSMs by identification run once, instead of rescanning the whole list per
  // ProteinIdentification. Keeps the per-prot_id scoping the counts below rely on. The original
  // index travels along so diagnostics can still name the offending PSM.
  std::unordered_map<std::string, std::vector<std::pair<size_t, const PeptideIdentification*>>> pep_by_identifier;
  for (size_t i = 0; i < peptide_identifications.size(); ++i)
  {
    pep_by_identifier[peptide_identifications[i].getIdentifier()].emplace_back(i, &peptide_identifications[i]);
  }

  // Iterate over all protein identifications; for each, build scoped lookup
  // structures and append one row per (indistinguishable protein group, run).
  Size total_groups_found = 0;
  Size groups_without_run_name = 0;   // run carries no spectra_data at all
  Size groups_without_evidence = 0;   // merged run, but nothing identifies which file
  Size unresolved_psms = 0;           // PSM in a run with files, but resolving to none of them
  PrimaryKeyWatch pk_watch;
  for (const auto& prot_id : protein_identifications)
  {
    const auto& groups = prot_id.getIndistinguishableProteins();
    if (groups.empty()) continue;
    total_groups_found += groups.size();

    // A single-file run has exactly one possible origin, so its groups can be keyed on that file
    // even without peptide evidence. Only merged runs need evidence to disambiguate.
    StringList run_paths;
    prot_id.getPrimaryMSRunPath(run_paths);
    const std::string sole_run =
      (run_paths.size() == 1) ? ArrowIOHelpers::qpxRunFileName(run_paths[0]) : std::string();

    // Build accession -> ProteinHit lookup scoped to this prot_id
    std::unordered_map<std::string, const ProteinHit*> hit_lookup;
    for (const auto& hit : prot_id.getHits())
    {
      hit_lookup[hit.getAccession()] = &hit;
    }

    // Peptide evidence for this prot_id, split by origin run: QPX defines the counts per
    // (group, run), so the run dimension has to be carried alongside the accession.
    // Keyed on the UNMODIFIED sequence, matching the ConsensusMap overload -- peptide_counts
    // counts sequences, while distinguishing peptidoforms is what feature_counts is for.
    std::unordered_map<std::string, std::map<std::string, std::unordered_set<std::string>>> acc_run_peptides;
    AccessionRuns acc_runs;
    auto bucket = pep_by_identifier.find(prot_id.getIdentifier());
    if (bucket != pep_by_identifier.end())
    {
      // A merged run whose PSMs carry no usable 'id_merge_index' resolves every one of them to
      // the run's FIRST file, silently mislabelling a primary-key column -- refuse instead, the
      // same guard the PSM view applies. Scoped to the PSMs this export actually reads: input
      // that contributes no protein group must not be able to abort an export it has no say in.
      for (const auto& [psm_index, pep_ptr] : bucket->second)
      {
        mapper.validateMergeIndex(*pep_ptr, psm_index);
      }

      for (const auto& [psm_index, pep_ptr] : bucket->second)
      {
        const auto& pep_id = *pep_ptr;
        if (pep_id.getHits().empty()) continue;
        const std::string run = ArrowIOHelpers::qpxRunFileName(mapper.getPrimaryMSRunPath(pep_id));
        if (run.empty())
        {
          // The run has files but this PSM resolves to none of them -- e.g. an 'id_merge_index'
          // left over from a merge that was later subset. Its evidence is dropped, which would
          // otherwise surface only as a silently zero peptide_count.
          if (!run_paths.empty()) { ++unresolved_psms; }
          continue;
        }
        // Use the best hit (first after assuming sorted by score)
        const auto& best_hit = pep_id.getHits()[0];
        std::string seq = best_hit.getSequence().toUnmodifiedString();
        for (const auto& ev : best_hit.getPeptideEvidences())
        {
          const std::string& acc = ev.getProteinAccession();
          if (!acc.empty())
          {
            acc_run_peptides[acc][run].insert(seq);
            acc_runs[acc].insert(run);
          }
        }
      }
    }

    /// Unmodified peptide sequences of @p acc seen in @p run, or nullptr when there are none.
    auto peptidesIn = [&acc_run_peptides](const std::string& acc, const std::string& run)
                        -> const std::unordered_set<std::string>*
    {
      auto ait = acc_run_peptides.find(acc);
      if (ait == acc_run_peptides.end()) { return nullptr; }
      auto rit = ait->second.find(run);
      return (rit == ait->second.end()) ? nullptr : &rit->second;
    };

    // Iterate protein groups for this prot_id and emit one row per run with evidence
    for (const auto& group : groups)
    {
      // grouped_runs is a QPX primary-key component that MUST NOT be null, and quantms
      // silently drops rows that violate that -- so a group that cannot be attributed to any
      // run is skipped with a diagnostic rather than emitted with an empty key.
      std::set<std::string> runs = evidenceRuns(acc_runs, group);
      if (runs.empty() && !sole_run.empty()) { runs.insert(sole_run); }
      if (runs.empty())
      {
        // Two distinct causes, reported separately below because the remedies differ.
        if (run_paths.empty()) { ++groups_without_run_name; }
        else                   { ++groups_without_evidence; }
        continue;
      }

      for (const auto& run_file : runs)
      {
        const std::string& anchor = group.accessions.empty() ? "" : group.accessions[0];
        pk_watch.add(anchor, {run_file}, std::nullopt);

        // anchor_protein
        (void)anchor_protein_builder.Append(anchor);

        // grouped_runs
        (void)grouped_runs_builder.Append();
        (void)grouped_runs_vb->Append(run_file);

        // pg_accessions
        (void)pg_accessions_builder.Append();
        for (const auto& acc : group.accessions)
        {
          (void)pg_acc_vb->Append(acc);
        }

        // pg_names
        (void)pg_names_builder.Append();
        for (const auto& acc : group.accessions)
        {
          auto it = hit_lookup.find(acc);
          if (it != hit_lookup.end())
          {
            (void)pg_names_vb->Append(it->second->getDescription());
          }
          else
          {
            (void)pg_names_vb->Append("");
          }
        }

        // gg_accessions, gg_names (gene groups — one entry per group member, keep
        // parallel to pg_accessions so consumers can zip positionally; empty string
        // when the meta value is missing).
        (void)gg_accessions_builder.Append();
        (void)gg_names_builder.Append();
        for (const auto& acc : group.accessions)
        {
          auto it = hit_lookup.find(acc);
          bool has_ga = (it != hit_lookup.end()) && it->second->metaValueExists("gene_accession");
          bool has_gn = (it != hit_lookup.end()) && it->second->metaValueExists("gene_name");
          (void)gg_acc_vb->Append(has_ga ? it->second->getMetaValue("gene_accession").toString() : "");
          (void)gg_names_vb->Append(has_gn ? it->second->getMetaValue("gene_name").toString() : "");
        }

        // gg_qvalue - not available
        (void)gg_qvalue_builder.AppendNull();

        // global_qvalue — use group probability
        (void)global_qvalue_builder.Append(group.probability);

        // pg_qvalue - not available for id-only
        (void)pg_qvalue_builder.AppendNull();

        // label/intensity — null for id-only
        (void)label_builder.AppendNull();
        (void)intensity_builder.AppendNull();

        // additional_intensities — null for id-only
        (void)additional_intensities_builder.AppendNull();

        // is_decoy
        auto anchor_it = hit_lookup.find(anchor);
        if (anchor_it != hit_lookup.end())
        {
          (void)is_decoy_builder.Append(anchor_it->second->isDecoy());
        }
        else
        {
          (void)is_decoy_builder.Append(false);
        }

        // contaminant - not tracked
        (void)contaminant_builder.AppendNull();

        // peptides — one {protein_name, peptide_count} per group member
        (void)peptides_builder.Append();
        int total_sequences = 0;
        std::unordered_set<std::string> group_unique_peptides;
        for (const auto& acc : group.accessions)
        {
          const auto* peps = peptidesIn(acc, run_file);
          int count = peps ? static_cast<int>(peps->size()) : 0;
          (void)pep_struct_b->Append();
          (void)pep_name_b->Append(acc);
          (void)pep_count_b->Append(count);
          total_sequences += count;
          if (peps) { group_unique_peptides.insert(peps->begin(), peps->end()); }
        }

        // peptide_counts — unique sequences across the union of group members,
        // total_sequences is the sum of per-protein counts (double-counts shared peptides).
        {
          (void)peptide_counts_builder.Append();
          (void)pc_unique_b->Append(static_cast<int>(group_unique_peptides.size()));
          (void)pc_total_b->Append(total_sequences);
        }

        // feature_counts - not available for id-only
        (void)feature_counts_builder.AppendNull();

        // sequence_coverage
        if (anchor_it != hit_lookup.end() && anchor_it->second->getCoverage() != ProteinHit::COVERAGE_UNKNOWN)
        {
          (void)sequence_coverage_builder.Append(static_cast<float>(anchor_it->second->getCoverage()));
        }
        else
        {
          (void)sequence_coverage_builder.AppendNull();
        }

        // molecular_weight - not available
        (void)molecular_weight_builder.AppendNull();

        // additional_scores — anchor protein score
        (void)additional_scores_builder.Append();
        if (anchor_it != hit_lookup.end())
        {
          (void)as_struct_b->Append();
          (void)as_name_b->Append(prot_id.getScoreType());
          (void)as_value_b->Append(anchor_it->second->getScore());
          (void)as_hb_b->Append(prot_id.isHigherScoreBetter());
        }

        // cv_params - empty for now
        (void)cv_params_builder.Append();
      } // per-run melt
    }
  }

  pk_watch.report("ProteinGroupArrowExport (id-only)");
  // WARN, not INFO: for a protein-inference result this is typically *every* group, and an empty
  // output file that was announced at INFO reads as "nothing to export" rather than "your input
  // could not be keyed".
  if (groups_without_run_name > 0)
  {
    OPENMS_LOG_WARN << "ProteinGroupArrowExport (id-only): skipped " << groups_without_run_name
                    << " protein group(s) whose identification run carries no MS run path, so no "
                       "grouped_runs can be derived -- and QPX makes it a primary-key component "
                       "that must not be empty. Protein inference commonly drops the path: set it "
                       "with ProteinIdentification::setPrimaryMSRunPath() before exporting."
                    << std::endl;
  }
  if (unresolved_psms > 0)
  {
    OPENMS_LOG_WARN << "ProteinGroupArrowExport (id-only): ignored " << unresolved_psms
                    << " PSM(s) that could not be resolved to any file of their identification "
                       "run, so they contribute to no peptide count. A stale 'id_merge_index' "
                       "left over from a merge is the usual cause." << std::endl;
  }
  if (groups_without_evidence > 0)
  {
    OPENMS_LOG_WARN << "ProteinGroupArrowExport (id-only): skipped " << groups_without_evidence
                    << " protein group(s) in a merged run with no peptide evidence identifying "
                       "which of its files they came from, so no row can be keyed. Annotate the "
                       "PSMs with 'id_merge_index', or export the runs separately." << std::endl;
  }

  if (total_groups_found == 0)
  {
    OPENMS_LOG_WARN << "ProteinGroupArrowExport (id-only): No indistinguishable protein groups found" << std::endl;
    return arrow::Table::MakeEmpty(target_schema).ValueOrDie();
  }

  // Finalize all arrays
  std::shared_ptr<arrow::Array> arr_pg_acc, arr_pg_names, arr_gg_acc, arr_gg_names;
  std::shared_ptr<arrow::Array> arr_gg_qval, arr_anchor, arr_grouped_runs;
  std::shared_ptr<arrow::Array> arr_global_qval, arr_pg_qval;
  std::shared_ptr<arrow::Array> arr_label, arr_intensity, arr_add_intensities;
  std::shared_ptr<arrow::Array> arr_is_decoy, arr_contaminant;
  std::shared_ptr<arrow::Array> arr_peptides, arr_pep_counts, arr_feat_counts;
  std::shared_ptr<arrow::Array> arr_seq_cov, arr_mol_weight;
  std::shared_ptr<arrow::Array> arr_add_scores, arr_cv_params;

  status = pg_accessions_builder.Finish(&arr_pg_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): pg_accessions Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_names_builder.Finish(&arr_pg_names);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): pg_names Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_accessions_builder.Finish(&arr_gg_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): gg_accessions Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_names_builder.Finish(&arr_gg_names);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): gg_names Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_qvalue_builder.Finish(&arr_gg_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): gg_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = anchor_protein_builder.Finish(&arr_anchor);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): anchor_protein Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = grouped_runs_builder.Finish(&arr_grouped_runs);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): grouped_runs Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = global_qvalue_builder.Finish(&arr_global_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): global_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_qvalue_builder.Finish(&arr_pg_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): pg_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = label_builder.Finish(&arr_label);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): label Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = intensity_builder.Finish(&arr_intensity);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): intensity Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_intensities_builder.Finish(&arr_add_intensities);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): additional_intensities Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Finish(&arr_is_decoy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): is_decoy Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = contaminant_builder.Finish(&arr_contaminant);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): contaminant Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptides_builder.Finish(&arr_peptides);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): peptides Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptide_counts_builder.Finish(&arr_pep_counts);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): peptide_counts Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = feature_counts_builder.Finish(&arr_feat_counts);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): feature_counts Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = sequence_coverage_builder.Finish(&arr_seq_cov);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): sequence_coverage Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = molecular_weight_builder.Finish(&arr_mol_weight);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): molecular_weight Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_scores_builder.Finish(&arr_add_scores);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): additional_scores Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = cv_params_builder.Finish(&arr_cv_params);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): cv_params Finish failed: " << status.ToString() << std::endl; return nullptr; }

  auto table = arrow::Table::Make(target_schema, {
    arr_pg_acc, arr_pg_names, arr_gg_acc, arr_gg_names,
    arr_gg_qval, arr_anchor, arr_grouped_runs,
    arr_global_qval, arr_pg_qval,
    arr_label, arr_intensity, arr_add_intensities,
    arr_is_decoy, arr_contaminant,
    arr_peptides, arr_pep_counts, arr_feat_counts,
    arr_seq_cov, arr_mol_weight,
    arr_add_scores, arr_cv_params,
  });

  // Appends discard their arrow::Status throughout this file, so a failed one (e.g. under memory
  // pressure) leaves that column shorter than the rest while every Finish() still succeeds
  // independently. Validate() catches the resulting length mismatch, so the contract of returning
  // nullptr on failure holds instead of handing back a malformed table.
  if (auto valid = table->Validate(); !valid.ok())
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: built an inconsistent table: "
                     << valid.ToString() << std::endl;
    return nullptr;
  }
  return table;
}


bool ProteinGroupArrowExport::exportToParquet(
  const std::vector<ProteinIdentification>& protein_identifications,
  const PeptideIdentificationList& peptide_identifications,
  const std::string& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportToArrow(protein_identifications, peptide_identifications);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): Failed to create Arrow table" << std::endl;
    return false;
  }
  return exportToParquet(table, filename, config);
}

} // namespace OpenMS
