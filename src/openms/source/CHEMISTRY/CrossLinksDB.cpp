// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/CrossLinksDB.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>

using namespace std;

namespace OpenMS
{
  CrossLinksDB::CrossLinksDB():
      ModificationsDB("","","") // if you clear the mods you don't need to read them in the first place! Also avoids easy memory leaks.
  {
    mods_.clear();
    modification_names_.clear();
    readFromOBOFile("CHEMISTRY/XLMOD.obo"); //TODO please comment why this is needed! Why not use the one from ModificationsDB
  }


  CrossLinksDB::~CrossLinksDB()
  {
    modification_names_.clear();
    for (vector<ResidueModification*>::iterator it = mods_.begin(); it != mods_.end(); ++it)
    {
      delete *it;
    }
  }

  void CrossLinksDB::readFromOBOFile(const String& filename)
  {
    ResidueModification mod;
    // add multiple mods for multiple specificities
    //Map<String, ResidueModification> all_mods;
    multimap<String, ResidueModification> all_mods;

    ifstream is(File::find(filename).c_str());
    String line, line_wo_spaces, id;
    String origin = "";

    bool reading_mono_link = false;

    //parse file
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
        // if the last [Term] was a moon-link, then it does not belong in CrossLinksDB
        if (!id.empty() && !reading_mono_link) //store last term
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

          if (origin.hasSubstring("ProteinN-term"))
          {
            mod.setTermSpecificity(ResidueModification::N_TERM);
            mod.setOrigin('X');
            all_mods.insert(make_pair(id, mod));
          }
          if (origin.hasSubstring("ProteinC-term"))
          {
            mod.setTermSpecificity(ResidueModification::C_TERM);
            mod.setOrigin('X');
            all_mods.insert(make_pair(id, mod));
          }

          id = "";
          origin = "";
          mod = ResidueModification();
        }
        else if (reading_mono_link) // re-initialize before reading next [Term]
        {
          id = "";
          origin = "";
          mod = ResidueModification();
          reading_mono_link = false;
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
          if (val_split[1] == "1")
          {
            reading_mono_link = true;
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

    if (!id.empty()) //store last term
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

      if (origin.hasSubstring("ProteinN-term"))
      {
        mod.setTermSpecificity(ResidueModification::N_TERM);
        mod.setOrigin('X');
        all_mods.insert(make_pair(id, mod));
      }
      if (origin.hasSubstring("ProteinC-term"))
      {
        mod.setTermSpecificity(ResidueModification::C_TERM);
        mod.setOrigin('X');
        all_mods.insert(make_pair(id, mod));
      }

      id = "";
      origin = "";
      mod = ResidueModification();
    }

    // now use the term and all synonyms to build the database
    for (multimap<String, ResidueModification>::const_iterator it = all_mods.begin(); it != all_mods.end(); ++it)
    {

      // check whether a unimod definition already exists, then simply add synonyms to it
      if (it->second.getUniModRecordId() > 0)
      {
        //cerr << "Found UniMod PSI-MOD mapping: " << it->second.getPSIMODAccession() << " " << it->second.getUniModAccession() << endl;
        set<const ResidueModification*> mods = modification_names_[it->second.getUniModAccession()];
        for (set<const ResidueModification*>::const_iterator mit = mods.begin(); mit != mods.end(); ++mit)
        {
          //cerr << "Adding PSIMOD accession: " << it->second.getPSIMODAccession() << " " << it->second.getUniModAccession() << endl;
          modification_names_[it->second.getPSIMODAccession()].insert(*mit);
        }
      }
      else
      {
        // the mod has so far not been mapped to a unimod mod
        // first check whether the mod is specific
        if ((it->second.getOrigin() != 'X') ||
            ((it->second.getTermSpecificity() != ResidueModification::ANYWHERE) &&
             (it->second.getDiffMonoMass() != 0)))
        {
          mods_.push_back(new ResidueModification(it->second));

          set<String> synonyms = it->second.getSynonyms();
          synonyms.insert(it->first);
          synonyms.insert(it->second.getFullName());
          //synonyms.insert(it->second.getUniModAccession());
          synonyms.insert(it->second.getPSIMODAccession());
          // full ID is auto-generated based on (short) ID, but we want the name instead:
          mods_.back()->setId(it->second.getFullName());
          mods_.back()->setFullId();
          mods_.back()->setId(it->second.getId());
          synonyms.insert(mods_.back()->getFullId());

          // now check each of the names and link it to the residue modification
          for (set<String>::const_iterator nit = synonyms.begin(); nit != synonyms.end(); ++nit)
          {
            modification_names_[*nit].insert(mods_.back());
          }
        }
      }
    }
  }

  void CrossLinksDB::getAllSearchModifications(vector<String>& modifications) const
  {
    modifications.clear();

    for (vector<ResidueModification*>::const_iterator it = mods_.begin(); it != mods_.end(); ++it)
    {
      if (!(*it)->getPSIMODAccession().empty())
      {
        modifications.push_back((*it)->getFullId());
      }
    }
    sort(modifications.begin(), modifications.end());
  }

}
