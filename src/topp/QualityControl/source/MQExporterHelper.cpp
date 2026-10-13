// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Virginia Rossow, Lenny Kovac, Hendrik Beschorner$
// --------------------------------------------------------------------------

#include "MQExporterHelper.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>

using namespace OpenMS;

Size MQExporterHelper::proteinGroupID_(std::map<std::string, OpenMS::Size>& database,
                                       const std::string& protein_accession)
{
  if (auto it = database.find(protein_accession); it == database.end())
  {
    database.emplace(protein_accession, database.size() + 1);
    return database.size();
  }
  else
  {
    return it->second;
  }
}

std::map<Size, Size> MQExporterHelper::makeFeatureUIDtoConsensusMapIndex_(const ConsensusMap& cmap)
{
  std::map<Size, Size> f_to_ci;
  for (Size i = 0; i < cmap.size(); ++i)
  {
    for (const auto& fh : cmap[i].getFeatures())
    {
      auto[it, was_created_newly] = f_to_ci.emplace(fh.getUniqueId(), i);
      if (!was_created_newly)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Adding [" + StringUtils::toStr(it->first) + "," + StringUtils::toStr(it->second) +  "] failed. FeatureHandle exists twice in ConsensusMap!");
      }
    }
  }
  return f_to_ci;
}


std::vector<MQExporterHelper::LinkedIdentification> MQExporterHelper::linkedIdentifications(const BaseFeature& feature, const IdentificationData& data)
{
  std::vector<LinkedIdentification> result;
  std::vector<double> scores; // of the best matches, NaN without one
  std::optional<bool> higher_better;
  for (const auto& entry : feature.getLinkedIdentifications(data))
  {
    LinkedIdentification item;
    item.observation = entry.query->getObservation();
    for (const auto& source : entry.run->getSources())
    {
      const auto& queries = source.identifications;
      if (!queries.empty() && entry.query >= &queries.front() && entry.query <= &queries.back())
      {
        item.file = source.file.path;
        break;
      }
    }
    double score = std::numeric_limits<double>::quiet_NaN();
    if (const auto* best = entry.getBestMatch())
    {
      const auto primary = *entry.run->getPrimaryScore();
      item.best = IdentificationDataAdapter::materializePeptide(*entry.run, *best, primary);
      score = *entry.run->getScore(best->getId(), primary);
      // the direction of the first scored identification, as for a list of peptide identifications
      if (!higher_better) higher_better = entry.run->getScoreDefinition(primary).higher_better;
    }
    result.push_back(std::move(item));
    scores.push_back(score);
  }
  std::vector<Size> order(result.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](Size a, Size b) {
    if (std::isnan(scores[a]) || std::isnan(scores[b])) return !std::isnan(scores[a]) && std::isnan(scores[b]);
    return *higher_better ? scores[a] > scores[b] : scores[a] < scores[b];
  });
  std::vector<LinkedIdentification> sorted;
  sorted.reserve(result.size());
  for (Size index : order) sorted.push_back(std::move(result[index]));
  return sorted;
}

const MQExporterHelper::LinkedIdentification* MQExporterHelper::bestIdentification(const std::vector<LinkedIdentification>& f_ids, const std::vector<LinkedIdentification>& cf_ids, std::string& type)
{
  if (!f_ids.empty() && f_ids[0].best)
  {
    const auto& best = f_ids[0];
    const bool in_consensus_feature = std::any_of(cf_ids.begin(), cf_ids.end(), [&](const LinkedIdentification& id) {
      return id.file == best.file && id.observation.data_id == best.observation.data_id;
    });
    if (in_consensus_feature)
    {
      type = "MULTI-MSMS";
      return &best;
    }
  }
  if (!cf_ids.empty() && cf_ids[0].best)
  {
    type = "MULTI-MATCH";
    return &cf_ids[0];
  }
  return nullptr;
}

bool MQExporterHelper::isValid(const std::string& filename)
{
  return File::writable(filename);
}

std::string MQExporterHelper::extractGeneName(const std::string& prot_description)
{
  std::string gene_name;
  auto pos_gene = prot_description.find("GN=");
  // description might not contain a gene name ...
  if (pos_gene == std::string::npos) return gene_name;

  gene_name = StringUtils::substr(prot_description, pos_gene +3, prot_description.find(" ", pos_gene +3));
  return gene_name;
}

MQExporterHelper::MQCommonOutputs::MQCommonOutputs(
  const OpenMS::Feature& f,
  const OpenMS::PeptideHit& best_hit,
  const OpenMS::MSExperiment& exp,
  const std::map<std::string,std::string>& prot_mapper)
{
  const OpenMS::PeptideHit* ptr_best_hit = &best_hit; // the best hit referring to score
  const OpenMS::AASequence& pep_seq = ptr_best_hit->getSequence();
  if (pep_seq.empty())
  {
    return; // empty AASequence; nothing to export
  }

  std::map<std::string, OpenMS::Size> modifications_temp;
  if (pep_seq.hasNTerminalModification())
  {
    const std::string& n_terminal_modification = pep_seq.getNTerminalModificationName();
    modifications_temp.emplace(std::make_pair(n_terminal_modification, 1));
  }
  if (pep_seq.hasCTerminalModification())
  {
    modifications_temp.emplace(std::make_pair(pep_seq.getCTerminalModificationName(), 1));
  }
  for (OpenMS::Size i = 0; i < pep_seq.size(); ++i)
  {
    if (pep_seq.getResidue(i).isModified())
    {
      ++modifications_temp[pep_seq.getResidue(i).getModification()->getFullId()];
    }
  }

  if (modifications_temp.empty())
  {
    modifications.str("Unmodified");
  }
  else
  {
    auto it = modifications_temp.begin();
    modifications.clear();
    modifications << it->first;
    ++it;
    for (; it != modifications_temp.end(); ++it)
    {
      modifications << ";" << it->first;
    }
  }

  acetyl = '0';
  if (pep_seq.hasNTerminalModification() && StringUtils::hasSubstring(pep_seq.getNTerminalModificationName(), "Acetyl"))
  {
    acetyl = '1'; // Acetyl (Protein N-term)
  }

  oxidation << modifications_temp["Oxidation (M)"];
  const std::set<std::string>& accessions = ptr_best_hit->extractProteinAccessionsSet();
  std::vector<std::string> gene_names_temp;
  std::vector<std::string> protein_names_temp;
  for(const auto& prot_access : accessions)
  {
    if (const auto& prot_mapper_it = prot_mapper.find(prot_access); prot_mapper_it == prot_mapper.end())
    {
      continue;
    }
    else
    {
      auto protein_description = prot_mapper_it->second;
      auto gn = extractGeneName(protein_description);
      if(!gn.empty())
      {
        gene_names_temp.push_back(std::move(gn));
      }
      protein_names_temp.push_back(std::move(protein_description));
    }
  }
  gene_names.str(ListUtils::concatenate(gene_names_temp, ";"));     //Gene Names
  protein_names.str(ListUtils::concatenate(protein_names_temp, ";"));  //Protein Names

  if (f.metaValueExists("spectrum_index") && !exp.empty() && exp.getNrSpectra() >= (OpenMS::Size)f.getMetaValue("spectrum_index") && !exp[f.getMetaValue("spectrum_index")].empty())
  {
    const OpenMS::MSSpectrum& ms2_spec = exp[f.getMetaValue("spectrum_index")];
    if (!ms2_spec.getPrecursors().empty())
    {
      msms_mz << ms2_spec.getPrecursors()[0].getMZ(); // MS/MS m/z
    }
  }
  const double& uncalibrated_mz_error_ppm = ptr_best_hit->getMetaValue("uncalibrated_mz_error_ppm", NAN);
  const double& calibrated_mz_error_ppm = ptr_best_hit->getMetaValue("calibrated_mz_error_ppm", NAN);

  if (std::isnan(uncalibrated_mz_error_ppm) && std::isnan(calibrated_mz_error_ppm))
  {
    uncalibrated_calibrated_mz_ppm.str("NA"); // Uncalibrated - Calibrated m/z [ppm]
    uncalibrated_calibrated_mz_mda.str("NA"); // Uncalibrated - Calibrated m/z [mDa]
    mass_error_ppm.str("NA");                 // Mass error [ppm]
    mass_error_da.str("NA");                  // Mass error [Da]
    uncalibrated_mass_error_ppm.str("NA");    // Uncalibrated Mass error [ppm]
    uncalibrated_mass_error_da.str("NA");     // Uncalibrated Mass error [Da]
  }
  else if (std::isnan(calibrated_mz_error_ppm))
  {
    uncalibrated_calibrated_mz_ppm.str("NA"); // Uncalibrated - Calibrated m/z [ppm]
    uncalibrated_calibrated_mz_mda.str("NA"); // Uncalibrated - Calibrated m/z [mDa]
    mass_error_ppm.str("NA");                 // Mass error [ppm]
    mass_error_da.str("NA");                  // Mass error [Da]
    uncalibrated_mass_error_ppm.clear();
    uncalibrated_mass_error_da.clear();
    uncalibrated_mass_error_ppm << uncalibrated_mz_error_ppm;                                    // Uncalibrated Mass error [ppm]
    uncalibrated_mass_error_da << OpenMS::Math::ppmToMass(uncalibrated_mz_error_ppm, f.getMZ()); // Uncalibrated Mass error [Da]
  }
  else if (std::isnan(uncalibrated_mz_error_ppm))
  {
    uncalibrated_calibrated_mz_ppm.str("NA"); // Uncalibrated - Calibrated m/z [ppm]
    uncalibrated_calibrated_mz_mda.str("NA"); // Uncalibrated - Calibrated m/z [mDa]
    mass_error_ppm.clear();
    mass_error_da.clear();
    mass_error_ppm << calibrated_mz_error_ppm;                                    // Mass error [ppm]
    mass_error_da << OpenMS::Math::ppmToMass(calibrated_mz_error_ppm, f.getMZ()); // Mass error [Da]
    uncalibrated_mass_error_ppm.str("NA");                                        // Uncalibrated Mass error [ppm]
    uncalibrated_mass_error_da.str("NA");                                         // Uncalibrated Mass error [Da]
  }
  else
  {
    uncalibrated_calibrated_mz_ppm.clear(); // Uncalibrated - Calibrated m/z [ppm]
    uncalibrated_calibrated_mz_mda.clear(); // Uncalibrated - Calibrated m/z [mDa]
    mass_error_ppm.clear();                 // Mass error [ppm]
    mass_error_da.clear();                  // Mass error [Da]
    uncalibrated_mass_error_ppm.clear();    // Uncalibrated Mass error [ppm]
    uncalibrated_mass_error_da.clear();     // Uncalibrated Mass error [Da]

    uncalibrated_calibrated_mz_ppm << uncalibrated_mz_error_ppm - calibrated_mz_error_ppm;                                       // Uncalibrated - Calibrated m/z [ppm]
    uncalibrated_calibrated_mz_mda << OpenMS::Math::ppmToMass((uncalibrated_mz_error_ppm - calibrated_mz_error_ppm), f.getMZ()); // Uncalibrated - Calibrated m/z [Da]
    mass_error_ppm << calibrated_mz_error_ppm;                                                                                   // Mass error [ppm]
    mass_error_da << OpenMS::Math::ppmToMass(calibrated_mz_error_ppm, f.getMZ());                                                // Mass error [Da]
    uncalibrated_mass_error_ppm << uncalibrated_mz_error_ppm;                                                                    // Uncalibrated Mass error [ppm]
    uncalibrated_mass_error_da << OpenMS::Math::ppmToMass(uncalibrated_mz_error_ppm, f.getMZ());                                 // Uncalibrated Mass error [Da]
  }

  base_peak_fraction.clear();
  if(f.metaValueExists("spectrum_index") && f.metaValueExists("base_peak_intensity") && !exp.empty())
  {
    const MSSpectrum& ms2_spec = exp[f.getMetaValue("spectrum_index")];
    if (!ms2_spec.getPrecursors().empty())
    {
      base_peak_fraction << (ms2_spec.getPrecursors()[0].getIntensity() / (double)f.getMetaValue("base_peak_intensity")); // Base peak fraction
    }
  }
}
