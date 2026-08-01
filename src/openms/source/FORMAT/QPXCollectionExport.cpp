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

#include <set>
#include <string>

namespace OpenMS
{

bool QPXCollectionExport::requireExportable(const ConsensusMap& cmap, const ExperimentalDesign& design)
{
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
  // The list is materialised rather than iterated in place because QPXFile's check takes one,
  // and the copies are of PeptideIdentification values the caller already holds. For the
  // streaming caller (IsobaricWorkflow) that is a real cost at millions of PSMs -- which is
  // exactly the caller that benefits most from failing before writing two huge files.
  {
    PeptideIdentificationList all_pep_ids;
    size_t n = cmap.getUnassignedPeptideIdentifications().size();
    for (const auto& feature : cmap) { n += feature.getPeptideIdentifications().size(); }
    all_pep_ids.reserve(n);
    for (const auto& feature : cmap)
    {
      for (const auto& pep_id : feature.getPeptideIdentifications()) { all_pep_ids.push_back(pep_id); }
    }
    for (const auto& pep_id : cmap.getUnassignedPeptideIdentifications()) { all_pep_ids.push_back(pep_id); }
    QPXFile::requireResolvableMergeIndices(cmap.getProteinIdentifications(), all_pep_ids);
  }

  // -- channel labels, shared by the feature and pg views --
  // intensities[].label is a documented join key against run.samples[].label. The feature view
  // refuses a channel it cannot name; the pg view invents a token. Refusing here makes the
  // collection consistent either way.
  std::set<std::pair<std::string, unsigned>> header_cells;
  for (const auto& [map_index, header] : cmap.getColumnHeaders())
  {
    (void)map_index;
    header_cells.emplace(ArrowIOHelpers::qpxRunFileName(header.filename),
                         header.getLabelAsUInt(cmap.getExperimentType()));
    const std::string channel_name = header.metaValueExists("channel_name")
                                   ? header.getMetaValue("channel_name").toString() : "";
    if (ArrowIOHelpers::qpxIntensityLabel(header.label, channel_name).empty())
    {
      // qpxIntensityLabel already logged which header and why.
      OPENMS_LOG_ERROR << "QPXCollectionExport: a column header names a channel QPX cannot "
                          "label, so intensities[].label would not join against run.parquet. "
                          "No QPX file was written." << std::endl;
      return false;
    }
  }

  // -- pg view --
  // Same QuantificationUnits the exporter builds, from the same inputs, so the two cannot
  // disagree about which designs are representable.
  const Internal::QuantificationUnits units(design, header_cells);
  const auto refuse = [](const std::string& what)
  {
    OPENMS_LOG_ERROR << "QPXCollectionExport: " << what
                     << ". The protein-group view cannot be written, so no QPX file was."
                     << std::endl;
    return false;
  };
  if (!units.missingCell().empty()) { return refuse(units.missingCell() + ", so the design and the data disagree about the channel layout"); }
  if (units.inconsistent()) { return refuse("the experimental design assigns several samples to one (run file, label) combination"); }
  if (!units.ambiguousLabel().empty()) { return refuse(units.ambiguousLabel()); }
  if (!units.raggedUnit().empty())     { return refuse(units.raggedUnit()); }
  if (!units.duplicatedSample().empty()) { return refuse(units.duplicatedSample()); }

  // The abundance array is indexed by the design's sample space; a length mismatch means the
  // caller passed a design other than the one that produced the quantification.
  if (!cmap.getProteinIdentifications().empty())
  {
    for (const auto& group : cmap.getProteinIdentifications()[0].getIndistinguishableProteins())
    {
      const auto* abundances = Internal::sampleAbundances(group);
      if (abundances == nullptr) { continue; }
      if (units.empty())
      {
        // Quantified, but the design yields no unit to attribute it to. The exporter refuses
        // rather than writing a null-intensity row that would deny the quantification.
        return refuse("a protein group carries sample abundances, but the experimental design "
                      "yields no quantification unit for them -- it lists no MS files this map "
                      "has a column header for");
      }
      if (abundances->size() != design.getNumberOfSamples())
      {
        return refuse("a protein group carries " + std::to_string(abundances->size())
                      + " sample abundances but the experimental design describes "
                      + std::to_string(design.getNumberOfSamples())
                      + " samples, so the design does not belong to this quantification");
      }
      // No early exit: the exporter checks EVERY group, so stopping at the first well-sized one
      // would let a later mismatch refuse the export after two files were already written --
      // exactly the failure this preflight exists to prevent.
    }
  }

  return true;
}

} // namespace OpenMS
