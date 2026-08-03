// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QPXCollectionExport.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport_impl.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <map>
#include <set>
#include <string>

namespace OpenMS
{

bool QPXCollectionExport::requireExportable(const ConsensusMap& cmap, const ExperimentalDesign& design)
{
  const auto refuse = [](const std::string& what)
  {
    OPENMS_LOG_ERROR << "QPXCollectionExport: " << what
                     << ". The protein-group view cannot be written, so no QPX file was."
                     << std::endl;
    return false;
  };

  // The pg exporter requires the first identification run and its inferred groups. Catch both
  // before the feature and PSM writers can leave a partial collection behind.
  if (cmap.getProteinIdentifications().empty())
  {
    return refuse("the ConsensusMap contains no protein identification run");
  }
  if (cmap.getProteinIdentifications()[0].getIndistinguishableProteins().empty())
  {
    return refuse("the first protein identification run contains no indistinguishable protein groups");
  }

  // -- feature view --
  // Both throw. Order matches the exporter's, so the first complaint a user sees is the same
  // one they would have seen without the preflight.
  ConsensusMapArrowExport::requireUnambiguousIdentities(cmap);
  ConsensusMapArrowExport::requireResolvableIdRuns(cmap);

  // -- psm view --
  // The psm exporter validates every identification the caller hands it: the assigned ones plus
  // the unassigned ones. Its set strictly contains the feature view's, which only checks each
  // feature's winning identification, so this is not redundant.
  //
  // Keep pointers into the map rather than deep-copying every PSM. The map outlives this
  // synchronous preflight, and QPXFile's pointer overload is the same check used by its streaming
  // writer.
  {
    std::vector<const PeptideIdentification*> all_pep_ids;
    size_t n = cmap.getUnassignedPeptideIdentifications().size();
    for (const auto& feature : cmap) { n += feature.getPeptideIdentifications().size(); }
    all_pep_ids.reserve(n);
    for (const auto& feature : cmap)
    {
      for (const auto& pep_id : feature.getPeptideIdentifications()) { all_pep_ids.push_back(&pep_id); }
    }
    for (const auto& pep_id : cmap.getUnassignedPeptideIdentifications()) { all_pep_ids.push_back(&pep_id); }
    QPXFile::requireResolvableMergeIndices(cmap.getProteinIdentifications(), all_pep_ids);
  }

  // -- channel labels, shared by the feature and pg views --
  // Feature intensities[].label and the pg scalar label both join against run.samples[].label.
  std::set<std::pair<std::string, unsigned>> header_cells;
  std::map<std::pair<std::string, UInt>, std::string> channel_labels;
  for (const auto& [map_index, header] : cmap.getColumnHeaders())
  {
    (void)map_index;
    const std::string run = ArrowIOHelpers::qpxRunFileName(header.filename);
    const UInt channel = header.getLabelAsUInt(cmap.getExperimentType());
    header_cells.emplace(run, channel);
    const std::string channel_name = header.metaValueExists("channel_name")
                                   ? header.getMetaValue("channel_name").toString() : "";
    const std::string label = ArrowIOHelpers::qpxIntensityLabel(header.label, channel_name);
    if (label.empty())
    {
      // qpxIntensityLabel already logged which header and why.
      OPENMS_LOG_ERROR << "QPXCollectionExport: a column header names a channel QPX cannot "
                          "label, so its quantity would not join against run.parquet. "
                          "No QPX file was written." << std::endl;
      return false;
    }
    channel_labels[{run, channel}] = label;
  }

  // -- pg view --
  // Same QuantificationUnits the exporter builds, from the same inputs, so the two cannot
  // disagree about which designs are representable.
  Internal::QuantificationUnits units(design, header_cells);
  if (units.inconsistent()) { return refuse("the experimental design assigns several samples to one (run file, label) combination"); }
  if (!units.invalidSample().empty()) { return refuse(units.invalidSample()); }
  if (!units.runInSeveralGroups().empty()) { return refuse(units.runInSeveralGroups()); }
  if (!units.missingCell().empty())
  {
    return refuse(units.missingCell()
                  + ", so the design and the data disagree about the channel layout");
  }
  if (!units.ambiguousLabel().empty()) { return refuse(units.ambiguousLabel()); }
  if (!units.raggedUnit().empty())     { return refuse(units.raggedUnit()); }
  if (!units.resolveLabels(
        [&channel_labels](const std::string& run, Int channel) -> std::string
        {
          auto it = channel_labels.find({run, static_cast<UInt>(channel)});
          return it == channel_labels.end() ? std::string() : it->second;
        },
        "QPXCollectionExport"))
  {
    return false;
  }

  std::set<std::pair<UInt, UInt>> expected_quantity_keys;
  for (const auto& unit : units.units())
  {
    for (const auto& channel : unit.channels)
    {
      expected_quantity_keys.emplace(unit.fraction_group, channel.label);
    }
  }

  // Validate every group's complete annotation contract. No early exit on a well-formed group:
  // a malformed later group must be found before the first collection file is opened.
  for (const auto& group : cmap.getProteinIdentifications()[0].getIndistinguishableProteins())
  {
    const auto* sample_abundances = Internal::sampleAbundances(group);
    const auto quantities = Internal::fractionGroupAbundances(group);
    if (!quantities.valid) { return refuse(quantities.error); }
    if (sample_abundances != nullptr && sample_abundances->size() != design.getNumberOfSamples())
    {
      return refuse("a protein group carries " + std::to_string(sample_abundances->size())
                    + " sample abundances but the experimental design describes "
                    + std::to_string(design.getNumberOfSamples())
                    + " samples, so the design does not belong to this quantification");
    }

    const bool quantified = sample_abundances != nullptr || quantities.present;
    if (!quantified) { continue; }
    if (!quantities.present)
    {
      return refuse("a protein group carries sample abundances but no fraction-group/label "
                    "abundances from the QPX 1.1 quantification path");
    }
    if (units.empty())
    {
      return refuse("a protein group carries fraction-group abundances, but the experimental "
                    "design yields no matching quantification unit");
    }
    if (quantities.values.size() != expected_quantity_keys.size())
    {
      return refuse("a protein group's fraction-group/label abundance count does not match the "
                    "experimental design");
    }
    for (const auto& key : expected_quantity_keys)
    {
      if (!quantities.values.count(key))
      {
        return refuse("a protein group has no quantity for fraction group "
                      + std::to_string(key.first) + ", label " + std::to_string(key.second));
      }
    }
  }

  return true;
}

} // namespace OpenMS
