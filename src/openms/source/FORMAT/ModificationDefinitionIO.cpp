// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ModificationDefinitionIO.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <algorithm>
#include <exception>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
  bool ModificationDefinitionIO::isDefinition(const ResidueModification* mod)
  {
    // the Id guard excludes anonymous mass shifts that merely inherited the DEFINED default
    return mod != nullptr && mod->getProvenance() == ResidueModification::DEFINED && !mod->getId().empty();
  }

  namespace
  {
    using DefsByRun = std::map<std::string, std::set<const ResidueModification*>>;

    // definitions named by the search space, whether or not any hit carries them
    void addSearchSpace_(DefsByRun& out, const std::vector<ProteinIdentification>& protein_ids)
    {
      const ModificationsDB* db = ModificationsDB::getInstance();
      for (const ProteinIdentification& prot : protein_ids)
      {
        const ProteinIdentification::SearchParameters& sp = prot.getSearchParameters();
        for (const std::vector<std::string>* names : {&sp.fixed_modifications, &sp.variable_modifications})
        {
          for (const std::string& name : *names)
          {
            if (!db->has(name)) continue; // has() is silent; searchModifications logs on a miss
            std::set<const ResidueModification*> matches; // a bare Id names every site it is defined on
            db->searchModifications(matches, name);
            for (const ResidueModification* mod : matches)
            {
              if (ModificationDefinitionIO::isDefinition(mod)) out[prot.getIdentifier()].insert(mod);
            }
          }
        }
      }
    }

    // definitions sitting on the peptidoforms of one identification. The set of its run is looked up once per
    // identification, and only if a hit carries a definition (the lookup by run name for every modified residue took
    // most of the 0.13 s that the 250,000 identifications of a timsTOF run took); a modification inserted just before is
    // not inserted again.
    void addHits_(DefsByRun& out, const PeptideIdentification& pep)
    {
      std::set<const ResidueModification*>* run_defs = nullptr;
      const ResidueModification* last = nullptr;
      const auto add = [&](const ResidueModification* mod)
      {
        if (mod == last || !ModificationDefinitionIO::isDefinition(mod)) return;
        if (run_defs == nullptr) run_defs = &out[pep.getIdentifier()];
        run_defs->insert(mod);
        last = mod;
      };
      for (const PeptideHit& hit : pep.getHits())
      {
        const AASequence& seq = hit.getSequence();
        if (seq.hasNTerminalModification()) add(seq.getNTerminalModification());
        for (Size i = 0; i < seq.size(); ++i)
        {
          if (seq[i].isModified()) add(seq[i].getModification());
        }
        if (seq.hasCTerminalModification()) add(seq.getCTerminalModification());
      }
    }

    void addFeature_(DefsByRun& out, const Feature& feature)
    {
      for (const PeptideIdentification& pep : feature.getPeptideIdentifications()) addHits_(out, pep);
      for (const Feature& sub : feature.getSubordinates()) addFeature_(out, sub);
    }
  } // namespace

  std::map<std::string, std::set<const ResidueModification*>> ModificationDefinitionIO::collect(
    const std::vector<ProteinIdentification>& protein_ids,
    const PeptideIdentificationList& peptide_ids)
  {
    DefsByRun out;
    addSearchSpace_(out, protein_ids);
    // The identifications in parallel, each thread a contiguous range of them, and the definitions of each run joined:
    // a union of sets, the same whatever the ranges. Walking the hits is bound by the memory latency (0.13 s on one
    // thread for the 250,000 identifications of a timsTOF run).
    const SignedSize n = static_cast<SignedSize>(peptide_ids.size());
    int num_parts = 1;
#ifdef _OPENMP
    num_parts = static_cast<int>(std::clamp<SignedSize>(n / 4096, 1, omp_get_max_threads()));
#endif
    if (num_parts == 1)
    {
      for (const PeptideIdentification& pep : peptide_ids) addHits_(out, pep);
      return out;
    }
    std::vector<DefsByRun> part_out(num_parts);
    std::vector<std::exception_ptr> errors(num_parts); // an exception must not leave the parallel region
#pragma omp parallel for schedule(static, 1) num_threads(num_parts)
    for (int part = 0; part < num_parts; ++part)
    {
      try
      {
        for (SignedSize i = n * part / num_parts; i < n * (part + 1) / num_parts; ++i) addHits_(part_out[part], peptide_ids[i]);
      }
      catch (...)
      {
        errors[part] = std::current_exception();
      }
    }
    for (const std::exception_ptr& error : errors)
    {
      if (error) std::rethrow_exception(error);
    }
    for (const DefsByRun& part : part_out)
    {
      for (const auto& [run, defs] : part) out[run].insert(defs.begin(), defs.end());
    }
    return out;
  }

  std::map<std::string, std::set<const ResidueModification*>> ModificationDefinitionIO::collect(const FeatureMap& map)
  {
    DefsByRun out;
    addSearchSpace_(out, map.getProteinIdentifications());
    for (const Feature& feature : map) addFeature_(out, feature);
    for (const PeptideIdentification& pep : map.getUnassignedPeptideIdentifications()) addHits_(out, pep);
    return out;
  }

  std::map<std::string, std::set<const ResidueModification*>> ModificationDefinitionIO::collect(const ConsensusMap& map)
  {
    DefsByRun out;
    addSearchSpace_(out, map.getProteinIdentifications());
    for (const ConsensusFeature& feature : map)
    {
      for (const PeptideIdentification& pep : feature.getPeptideIdentifications()) addHits_(out, pep);
    }
    for (const PeptideIdentification& pep : map.getUnassignedPeptideIdentifications()) addHits_(out, pep);
    return out;
  }

  std::map<std::string, std::string> ModificationDefinitionIO::encodeByRun(
    const std::vector<ProteinIdentification>& protein_ids,
    const std::map<std::string, std::set<const ResidueModification*>>& defs)
  {
    std::map<std::string, std::string> out;
    for (const ProteinIdentification& prot : protein_ids)
    {
      ProteinIdentification::SearchParameters sp = prot.getSearchParameters();
      const auto it = defs.find(prot.getIdentifier());
      attach(sp, it == defs.end() ? std::set<const ResidueModification*>() : it->second);
      if (sp.metaValueExists(Constants::UserParam::MODIFICATION_DEFINITIONS))
      {
        out[prot.getIdentifier()] = sp.getMetaValue(Constants::UserParam::MODIFICATION_DEFINITIONS).toString();
      }
    }
    return out;
  }

  std::string ModificationDefinitionIO::encode(const std::set<const ResidueModification*>& defs)
  {
    // sorted, so the same definitions always produce the same bytes
    std::vector<std::string> records;
    records.reserve(defs.size());
    for (const ResidueModification* mod : defs)
    {
      records.push_back(mod->toDefinitionString());
    }
    std::sort(records.begin(), records.end());
    std::string out;
    for (const std::string& r : records)
    {
      if (!out.empty()) out += ';';
      out += r;
    }
    return out;
  }

  void ModificationDefinitionIO::attach(ProteinIdentification::SearchParameters& sp,
                                        const std::set<const ResidueModification*>& defs)
  {
    const std::string& key = Constants::UserParam::MODIFICATION_DEFINITIONS;
    std::set<std::string> records;
    if (sp.metaValueExists(key))
    {
      for (const std::string& r : ResidueModification::splitDefinitionRecords(sp.getMetaValue(key).toString()))
      {
        records.insert(r);
      }
    }
    for (const ResidueModification* mod : defs)
    {
      records.insert(mod->toDefinitionString());
    }
    if (records.empty()) return;
    std::string out;
    for (const std::string& r : records) // std::set iterates sorted: deterministic
    {
      if (!out.empty()) out += ';';
      out += r;
    }
    sp.setMetaValue(key, out);
  }

  Size ModificationDefinitionIO::registerFrom(const std::string& records)
  {
    const ModificationsDB* db = ModificationsDB::getInstance();
    Size n = 0;
    for (const std::string& record : ResidueModification::splitDefinitionRecords(records))
    {
      try
      {
        db->registerDefinition(ResidueModification::fromDefinitionString(record));
        ++n;
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_WARN << "Skipping an unreadable modification definition record (" << e.getMessage()
                        << "): '" << record << "'" << std::endl;
      }
    }
    return n;
  }

  Size ModificationDefinitionIO::registerFrom(const ProteinIdentification::SearchParameters& sp)
  {
    const std::string& key = Constants::UserParam::MODIFICATION_DEFINITIONS;
    if (!sp.metaValueExists(key)) return 0;
    return registerFrom(sp.getMetaValue(key).toString());
  }
} // namespace OpenMS
