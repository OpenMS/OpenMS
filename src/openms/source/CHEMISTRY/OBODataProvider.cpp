// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/OBODataProvider.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

namespace OpenMS
{

  OBODataProvider::OBODataProvider(const String& filename, bool cross_links_only)
    : filename_(filename),
      cross_links_only_(cross_links_only)
  {
  }

  std::vector<std::unique_ptr<ResidueModification>> OBODataProvider::loadModifications()
  {
    ResidueModification mod;
    multimap<String, ResidueModification> all_mods;

    String resolved = File::find(filename_);
    ifstream is(resolved.c_str());
    if (!is)
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, resolved);
    }
    String line, line_wo_spaces, id;
    String origin = "";

    bool skip_current_term = false;

    // Phase 1: parse OBO file
    while (getline(is, line, '\n'))
    {
      line.trim();
      line_wo_spaces = line;
      line_wo_spaces.removeWhitespaces();

      if (line.empty() || line[0] == '!') //skip empty lines and comments
      {
        continue;
      }

      if (line_wo_spaces == "[Term]")       //new term
      {
        if (!id.empty() && !skip_current_term) //store last term
        {
          // split into single residues and make unique (for XL-MS, where equal specificities for both sides are possible)
          vector<String> origins;
          origin.split(",", origins);

          std::sort(origins.begin(), origins.end());
          vector<String>::iterator unique_end = unique(origins.begin(), origins.end());
          origins.resize(distance(origins.begin(), unique_end));

          for (vector<String>::iterator orig_it = origins.begin(); orig_it != origins.end(); ++orig_it)
          {
            // we don't allow modifications with ambiguity codes as origin (except "X"):
            if ((orig_it->size() == 1) && (*orig_it != "B") && (*orig_it != "J") && (*orig_it != "Z"))
            {
              mod.setOrigin((*orig_it)[0]);
              all_mods.insert(make_pair(id, mod));
            }
          }

          // for mono-links from XLMOD.obo / cross-linker terminal specificities:
          if (origin.hasSubstring("ProteinN-term"))
          {
            mod.setTermSpecificity(cross_links_only_ ? ResidueModification::N_TERM : ResidueModification::PROTEIN_N_TERM);
            mod.setOrigin('X');
            all_mods.insert(make_pair(id, mod));
          }
          if (origin.hasSubstring("ProteinC-term"))
          {
            mod.setTermSpecificity(cross_links_only_ ? ResidueModification::C_TERM : ResidueModification::PROTEIN_C_TERM);
            mod.setOrigin('X');
            all_mods.insert(make_pair(id, mod));
          }

          id = "";
          origin = "";
          mod = ResidueModification();
        }
        else if (skip_current_term) // re-initialize before reading next [Term]
        {
          id = "";
          origin = "";
          mod = ResidueModification();
          skip_current_term = false;
        }
      }

      //new id line
      else if (line_wo_spaces.hasPrefix("id:"))
      {
        id = line.substr(line.find(':') + 1).trim();
        mod.setId(id);
        mod.setPSIMODAccession(id);
      }
      else if (line_wo_spaces.hasPrefix("name:"))
      {
        String name = line.substr(line.find(':') + 1).trim();
        mod.setFullName(name);
        if (mod.getId().hasSubstring("XLMOD"))
        {
          mod.setName(name);
          mod.setId(name);
          mod.setFullName(name);
        }
      }
      else if (line_wo_spaces.hasPrefix("is_a:"))
      {
        // TODO
      }
      else if (line_wo_spaces.hasPrefix("def:"))
      {
        line.remove('[');
        line.remove(']');
        line.remove(',');
        vector<String> split;
        line.split(' ', split);
        for (Size i = 0; i != split.size(); ++i)
        {
          if (split[i].hasPrefix("UniMod:"))
          {
            // Parse UniMod identifier to int
            String identifier = split[i].substr(7, split[i].size());
            mod.setUniModRecordId(identifier.toInt());
          }
        }
      }
      else if (line_wo_spaces.hasPrefix("comment:"))
      {
        // TODO
      }
      else if (line_wo_spaces.hasPrefix("synonym:"))
      {
        vector<String> val_split;
        line.split('"', val_split);
        if (val_split.size() < 3)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "missing \" characters to enclose argument!");
        }
        mod.addSynonym(val_split[1]);

        if (line_wo_spaces.hasSubstring("PSI-MOD-label"))
        {
          mod.setName(val_split[1]);
        }
      }
      else if (line_wo_spaces.hasPrefix("property_value:"))
      {
        String val = line_wo_spaces.substr(15, line_wo_spaces.size() - 15);
        val.trim();

        if (val.hasSubstring("\"none\""))
        {
          continue;
        }

        vector<String> val_split;
        val.split('"', val_split);
        if (val_split.size() != 3)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "missing \" characters to enclose argument!");
        }
        if (val.hasPrefix("DiffAvg:"))
        {
          mod.setDiffAverageMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("DiffFormula:"))
        {
          vector<String> tmp_split;
          line.split('"', tmp_split);
          tmp_split[1].removeWhitespaces();
          mod.setDiffFormula(EmpiricalFormula(tmp_split[1]));
        }
        else if (val.hasPrefix("DiffMono:"))
        {
          mod.setDiffMonoMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("Formula:"))
        {
          mod.setFormula(val_split[1]);
        }
        else if (val.hasPrefix("MassAvg:"))
        {
          mod.setAverageMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("MassMono:"))
        {
          mod.setMonoMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("Origin:"))
        {
          //mod.setOrigin(val_split[1]);
          origin = val_split[1];
        }
        else if (val.hasPrefix("Source:"))
        {
          mod.setSourceClassification(val_split[1]);
        }
        else if (val.hasPrefix("TermSpec:"))
        {
          mod.setTermSpecificity(val_split[1]);
        }
        // XLMOD specific fields
        else if (val.hasPrefix("reactionSites:"))
        {
          if (val_split[1] == "2" && !cross_links_only_)
          {
            skip_current_term = true; // skip cross-linkers when loading for ModificationsDB
          }
          if (val_split[1] == "1" && cross_links_only_)
          {
            skip_current_term = true; // skip mono-links when loading for CrossLinksDB
          }
        }
        else if (val.hasPrefix("monoisotopicMass:"))
        {
          mod.setDiffMonoMass(val_split[1].toDouble());
        }
        else if (val.hasPrefix("specificities:"))
        {
          // TODO cross-linker specificities can be different for both chain sides, right now the union of both sides is used
          // Input parameters of the cross-link search tool make sure, that the chemistry is not violated
          origin = val_split[1];

          // remove brackets
          origin.remove('(');
          origin.remove(')');
          origin.substitute("&", ",");
        }
      }
    }

    if (!id.empty() && !skip_current_term) //store last term
    {
      // split into single residues and make unique (for XL-MS, where equal specificities for both sides are possible)
      vector<String> origins;
      origin.split(",", origins);

      std::sort(origins.begin(), origins.end());
      vector<String>::iterator unique_end = unique(origins.begin(), origins.end());
      origins.resize(distance(origins.begin(), unique_end));

      for (vector<String>::iterator orig_it = origins.begin(); orig_it != origins.end(); ++orig_it)
      {
        // we don't allow modifications with ambiguity codes as origin (except "X"):
        if ((orig_it->size() == 1) && (*orig_it != "B") && (*orig_it != "J") && (*orig_it != "Z"))
        {
          mod.setOrigin((*orig_it)[0]);
          all_mods.insert(make_pair(id, mod));
        }
      }

      // for mono-links from XLMOD.obo / cross-linker terminal specificities:
      if (origin.hasSubstring("ProteinN-term"))
      {
        mod.setTermSpecificity(cross_links_only_ ? ResidueModification::N_TERM : ResidueModification::PROTEIN_N_TERM);
        mod.setOrigin('X');
        all_mods.insert(make_pair(id, mod));
      }
      if (origin.hasSubstring("ProteinC-term"))
      {
        mod.setTermSpecificity(cross_links_only_ ? ResidueModification::C_TERM : ResidueModification::PROTEIN_C_TERM);
        mod.setOrigin('X');
        all_mods.insert(make_pair(id, mod));
      }
    }

    // Phase 2: build ResidueModification objects from parsed terms
    vector<unique_ptr<ResidueModification>> result;

    for (multimap<String, ResidueModification>::const_iterator it = all_mods.begin(); it != all_mods.end(); ++it)
    {
      // PSI-MOD entries that map to UniMod: return them as-is.
      // ModificationsDB's loadFromProviders_() will detect these and handle alias resolution.
      if (it->second.getUniModRecordId() > 0)
      {
        result.push_back(make_unique<ResidueModification>(it->second));
        continue;
      }

      // the mod has so far not been mapped to a unimod mod
      // first check whether the mod is specific
      if ((it->second.getOrigin() != 'X') ||
          ((it->second.getTermSpecificity() != ResidueModification::ANYWHERE) &&
           (it->second.getDiffMonoMass() != 0)))
      {
        auto mod_ptr = make_unique<ResidueModification>(it->second);

        // full ID is auto-generated based on (short) ID, but we want the name instead:
        mod_ptr->setId(it->second.getFullName());
        mod_ptr->setFullId();
        mod_ptr->setId(it->second.getId());

        result.push_back(std::move(mod_ptr));
      }
    }

    return result;
  }

} // namespace OpenMS
