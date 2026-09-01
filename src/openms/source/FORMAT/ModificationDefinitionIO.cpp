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

#include <algorithm>

namespace OpenMS
{
  bool ModificationDefinitionIO::isDefinition(const ResidueModification* mod)
  {
    // the Id guard excludes anonymous mass shifts that merely inherited the DEFINED default
    return mod != nullptr && mod->getProvenance() == ResidueModification::DEFINED && !mod->getId().empty();
  }

  std::map<std::string, std::set<const ResidueModification*>> ModificationDefinitionIO::collect(
    const std::vector<ProteinIdentification>& protein_ids,
    const PeptideIdentificationList& peptide_ids)
  {
    std::map<std::string, std::set<const ResidueModification*>> out;
    const ModificationsDB* db = ModificationsDB::getInstance();

    // Definitions named by the search space, whether or not any hit carries them.
    for (const ProteinIdentification& prot : protein_ids)
    {
      const ProteinIdentification::SearchParameters& sp = prot.getSearchParameters();
      for (const std::vector<std::string>* names : {&sp.fixed_modifications, &sp.variable_modifications})
      {
        for (const std::string& name : *names)
        {
          if (!db->has(name)) continue; // has() is silent; searchModificationsFast logs on a miss
          bool multiple_matches = false;
          const ResidueModification* mod = db->searchModificationsFast(name, multiple_matches);
          if (isDefinition(mod)) out[prot.getIdentifier()].insert(mod);
        }
      }
    }

    // Definitions actually sitting on peptidoforms.
    for (const PeptideIdentification& pep : peptide_ids)
    {
      const std::string& run = pep.getIdentifier();
      for (const PeptideHit& hit : pep.getHits())
      {
        const AASequence& seq = hit.getSequence();
        if (seq.hasNTerminalModification() && isDefinition(seq.getNTerminalModification()))
        {
          out[run].insert(seq.getNTerminalModification());
        }
        for (Size i = 0; i < seq.size(); ++i)
        {
          if (seq[i].isModified() && isDefinition(seq[i].getModification()))
          {
            out[run].insert(seq[i].getModification());
          }
        }
        if (seq.hasCTerminalModification() && isDefinition(seq.getCTerminalModification()))
        {
          out[run].insert(seq.getCTerminalModification());
        }
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
