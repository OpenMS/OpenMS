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
#include <filesystem>
#include <future>
#include <thread>

using namespace std;

namespace OpenMS
{

  OBODataProvider::OBODataProvider(const std::string& filename, bool cross_links_only)
    : filename_(filename),
      cross_links_only_(cross_links_only)
  {
  }

  namespace
  {
    /// Phase 1 of OBODataProvider::loadModifications(): reads the lines of an OBO file, or of a part of one that starts
    /// with a [Term] line, and collects the terms in #all_mods
    class TermReader
    {
    public:
      explicit TermReader(bool cross_links_only) : cross_links_only_(cross_links_only) {}

      void line(std::string line)
      {
        StringUtils::trim(line);
        line_wo_spaces = line;
        StringUtils::removeWhitespaces(line_wo_spaces);

        if (line.empty() || line[0] == '!') //skip empty lines and comments
        {
          return;
        }

        if (line_wo_spaces == "[Term]")       //new term
        {
          if (!id.empty() && !skip_current_term) //store last term
          {
            storeTerm_();
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
        else if (StringUtils::hasPrefix(line_wo_spaces, "id:"))
        {
          id = StringUtils::trimmed(StringUtils::substr(line, line.find(':') + 1));
          mod.setId(id);
          mod.setPSIMODAccession(id);
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "name:"))
        {
          std::string name = StringUtils::trimmed(StringUtils::substr(line, line.find(':') + 1));
          mod.setFullName(name);
          if (StringUtils::hasSubstring(mod.getId(), "XLMOD"))
          {
            mod.setName(name);
            mod.setId(name);
            mod.setFullName(name);
          }
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "is_a:"))
        {
          // TODO
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "def:"))
        {
          StringUtils::remove(line, '[');
          StringUtils::remove(line, ']');
          StringUtils::remove(line, ',');
          vector<std::string> split;
          StringUtils::split(line, ' ', split);
          for (Size i = 0; i != split.size(); ++i)
          {
            if (StringUtils::hasPrefix(split[i], "UniMod:"))
            {
              // Parse UniMod identifier to int
              std::string identifier = split[i].substr(7, split[i].size());
              mod.setUniModRecordId(StringUtils::toInt32(identifier));
            }
          }
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "comment:"))
        {
          // TODO
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "synonym:"))
        {
          vector<std::string> val_split;
          StringUtils::split(line, '"', val_split);
          if (val_split.size() < 3)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "missing \" characters to enclose argument!");
          }
          mod.addSynonym(val_split[1]);

          if (StringUtils::hasSubstring(line_wo_spaces, "PSI-MOD-label"))
          {
            mod.setName(val_split[1]);
          }
        }
        else if (StringUtils::hasPrefix(line_wo_spaces, "property_value:"))
        {
          std::string val = StringUtils::substr(line_wo_spaces, 15, line_wo_spaces.size() - 15);
          StringUtils::trim(val);

          if (StringUtils::hasSubstring(val, "\"none\""))
          {
            return;
          }

          vector<std::string> val_split;
          StringUtils::split(val, '"', val_split);
          if (val_split.size() != 3)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "missing \" characters to enclose argument!");
          }
          if (StringUtils::hasPrefix(val, "DiffAvg:"))
          {
            mod.setDiffAverageMass(StringUtils::toDouble(val_split[1]));
          }
          else if (StringUtils::hasPrefix(val, "DiffFormula:"))
          {
            vector<std::string> tmp_split;
            StringUtils::split(line, '"', tmp_split);
            StringUtils::removeWhitespaces(tmp_split[1]);
            mod.setDiffFormula(EmpiricalFormula(tmp_split[1]));
          }
          else if (StringUtils::hasPrefix(val, "DiffMono:"))
          {
            mod.setDiffMonoMass(StringUtils::toDouble(val_split[1]));
          }
          else if (StringUtils::hasPrefix(val, "Formula:"))
          {
            mod.setFormula(val_split[1]);
          }
          else if (StringUtils::hasPrefix(val, "MassAvg:"))
          {
            mod.setAverageMass(StringUtils::toDouble(val_split[1]));
          }
          else if (StringUtils::hasPrefix(val, "MassMono:"))
          {
            mod.setMonoMass(StringUtils::toDouble(val_split[1]));
          }
          else if (StringUtils::hasPrefix(val, "Origin:"))
          {
            //mod.setOrigin(val_split[1]);
            origin = val_split[1];
          }
          else if (StringUtils::hasPrefix(val, "Source:"))
          {
            mod.setSourceClassification(val_split[1]);
          }
          else if (StringUtils::hasPrefix(val, "TermSpec:"))
          {
            mod.setTermSpecificity(val_split[1]);
          }
          // XLMOD specific fields
          else if (StringUtils::hasPrefix(val, "reactionSites:"))
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
          else if (StringUtils::hasPrefix(val, "monoisotopicMass:"))
          {
            mod.setDiffMonoMass(StringUtils::toDouble(val_split[1]));
          }
          else if (StringUtils::hasPrefix(val, "specificities:"))
          {
            // TODO cross-linker specificities can be different for both chain sides, right now the union of both sides is used
            // Input parameters of the cross-link search tool make sure, that the chemistry is not violated
            origin = val_split[1];

            // remove brackets
            StringUtils::remove(origin, '(');
            StringUtils::remove(origin, ')');
            StringUtils::substitute(origin, "&", ",");
          }
        }
    
      }

      /// After the last line: stores the last term (also at the end of a part: the next part begins with a [Term] line)
      void end()
      {
        if (!id.empty() && !skip_current_term) //store last term
        {
          storeTerm_();
        }
      }

      /// Is the state after the last line the one that a following [Term] line would start from? Not for lines without
      /// an id line since the last [Term] (e.g. the header of the file): those carry over into the next term.
      bool endsTerm() const { return !id.empty() || skip_current_term; }

      multimap<std::string, ResidueModification> all_mods;

    private:
      void storeTerm_()
      {
        // split into single residues and make unique (for XL-MS, where equal specificities for both sides are possible)
        vector<std::string> origins;
        StringUtils::split(origin, ",", origins);

        std::sort(origins.begin(), origins.end());
        vector<std::string>::iterator unique_end = unique(origins.begin(), origins.end());
        origins.resize(distance(origins.begin(), unique_end));

        for (vector<std::string>::iterator orig_it = origins.begin(); orig_it != origins.end(); ++orig_it)
        {
          // we don't allow modifications with ambiguity codes as origin (except "X"):
          if ((orig_it->size() == 1) && (*orig_it != "B") && (*orig_it != "J") && (*orig_it != "Z"))
          {
            mod.setOrigin((*orig_it)[0]);
            all_mods.insert(make_pair(id, mod));
          }
        }

        // for mono-links from XLMOD.obo / cross-linker terminal specificities:
        if (StringUtils::hasSubstring(origin, "ProteinN-term"))
        {
          mod.setTermSpecificity(cross_links_only_ ? ResidueModification::N_TERM : ResidueModification::PROTEIN_N_TERM);
          mod.setOrigin('X');
          all_mods.insert(make_pair(id, mod));
        }
        if (StringUtils::hasSubstring(origin, "ProteinC-term"))
        {
          mod.setTermSpecificity(cross_links_only_ ? ResidueModification::C_TERM : ResidueModification::PROTEIN_C_TERM);
          mod.setOrigin('X');
          all_mods.insert(make_pair(id, mod));
        }
      }

      bool cross_links_only_;
      ResidueModification mod;
      std::string line_wo_spaces, id;
      std::string origin;
      bool skip_current_term = false;
    };

    /**
      Phase 1 for a large regular file: reads parts of it that start with a [Term] line on several threads, each with
      its own TermReader, and collects their terms in file order (a multimap keeps equal keys in the order of
      insertion, also when merging). False if the file is small or not a regular file, a part does not end a term
      (lines without an id since the last [Term] line, which carry over into the next term), a part fails to parse or
      anything else fails (e.g. memory): then @p all_mods is unchanged and the caller reads the file as a whole, which
      also reports errors as before.
    */
    bool readInParts(const std::string& file, bool cross_links_only, multimap<std::string, ResidueModification>& all_mods)
    {
      const Size max_parts = 8;
      const Size min_part_bytes = 128 * 1024;
      const Size n_threads = std::thread::hardware_concurrency();
      if (n_threads < 2) return false;
      std::error_code ec;
      if (!std::filesystem::is_regular_file(file, ec)) return false;
      const auto size = std::filesystem::file_size(file, ec);
      if (ec) return false;
      const Size n_parts = std::min({max_parts, n_threads, Size(size) / min_part_bytes});
      if (n_parts < 2) return false;
      try
      {
        std::string text;
        {
          std::ifstream in(file, std::ios::binary);
          if (!in) return false;
          text.resize(size_t(size));
          if (!in.read(text.data(), std::streamsize(size)) || in.peek() != std::ifstream::traits_type::eof()) return false;
        }

        // parts start at lines that are exactly "[Term]" (PSI-MOD.obo has CRLF line ends)
        std::vector<size_t> starts{0};
        for (Size k = 1; k < n_parts; ++k)
        {
          size_t p = std::max(starts.back() + 1, text.size() * k / n_parts);
          for (p = text.find("\n[Term]", p - 1); p != std::string::npos; p = text.find("\n[Term]", p + 1))
          {
            const size_t after = p + 7;
            if (after == text.size() || text[after] == '\n' || (text[after] == '\r' && (after + 1 == text.size() || text[after + 1] == '\n'))) break;
          }
          if (p == std::string::npos) break;
          starts.push_back(p + 1);
        }
        starts.push_back(text.size());
        if (starts.size() < 3) return false;

        std::vector<TermReader> readers(starts.size() - 1, TermReader(cross_links_only));
        auto read_part = [&](Size i)
        {
          TermReader& reader = readers[i];
          for (size_t from = starts[i]; from < starts[i + 1];)
          {
            size_t to = text.find('\n', from);
            if (to == std::string::npos || to > starts[i + 1]) to = starts[i + 1];
            reader.line(text.substr(from, to - from));
            from = to + 1;
          }
          const bool ends_term = i + 1 == readers.size() || reader.endsTerm();
          reader.end();
          return ends_term;
        };
        std::vector<std::future<bool>> later;
        for (Size i = 1; i < readers.size(); ++i)
        {
          later.push_back(std::async(std::launch::async | std::launch::deferred, read_part, i));
        }
        bool ok = true;
        try
        {
          ok = read_part(0);
        }
        catch (...)
        {
          ok = false;
        }
        for (auto& f : later)
        {
          try
          {
            ok = f.get() && ok;
          }
          catch (...)
          {
            ok = false;
          }
        }
        if (!ok) return false;
        multimap<std::string, ResidueModification> merged;
        for (TermReader& reader : readers)
        {
          merged.merge(reader.all_mods); // moves the nodes, equal keys after those already there
        }
        all_mods = std::move(merged);
        return true;
      }
      catch (...)
      {
        return false;
      }
    }
  }

  std::vector<std::unique_ptr<ResidueModification>> OBODataProvider::loadModifications()
  {
    std::string resolved = File::find(filename_);
    ifstream is(resolved.c_str());
    if (!is)
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, resolved);
    }

    // Phase 1: parse OBO file
    multimap<std::string, ResidueModification> all_mods;
    if (!readInParts(resolved, cross_links_only_, all_mods))
    {
      TermReader reader(cross_links_only_);
      std::string line;
      while (getline(is, line, '\n'))
      {
        reader.line(line);
      }
      reader.end();
      all_mods = std::move(reader.all_mods);
    }

    // Phase 2: build ResidueModification objects from parsed terms
    vector<unique_ptr<ResidueModification>> result;

    for (multimap<std::string, ResidueModification>::const_iterator it = all_mods.begin(); it != all_mods.end(); ++it)
    {
      // PSI-MOD entries that map to UniMod: return them as-is.
      // ModificationsDB's loadFromProviders_() will detect these and handle alias resolution.
      if (it->second.getUniModRecordId() > 0)
      {
        auto alias_ptr = make_unique<ResidueModification>(it->second);
        alias_ptr->setProvenance(ResidueModification::CV);
        result.push_back(std::move(alias_ptr));
        continue;
      }

      // the mod has so far not been mapped to a unimod mod
      // first check whether the mod is specific
      if ((it->second.getOrigin() != 'X') ||
          ((it->second.getTermSpecificity() != ResidueModification::ANYWHERE) &&
           (it->second.getDiffMonoMass() != 0)))
      {
        auto mod_ptr = make_unique<ResidueModification>(it->second);
        mod_ptr->setProvenance(ResidueModification::CV); // PSI-MOD / XLMOD are shipped vocabularies

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
