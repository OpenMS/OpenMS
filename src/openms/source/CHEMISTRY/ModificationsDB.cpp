// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>

using namespace std;

namespace OpenMS
{
  bool ModificationsDB::is_instantiated_ = false;
  

  ModificationsDB::ModificationsDB(OpenMS::String unimod_file, OpenMS::String psimod_file, OpenMS::String xlmod_file)
  {
    if (!unimod_file.empty())
    {
      readFromUnimodXMLFile(unimod_file);
    }

    if (!psimod_file.empty())
    {
      readFromOBOFile(psimod_file);
    }

    if (!xlmod_file.empty())
    {
      readFromOBOFile(xlmod_file);
    }

    is_instantiated_ = true;
  }


  bool ModificationsDB::isInstantiated()
  {
    return is_instantiated_;
  }


  ModificationsDB::~ModificationsDB()
  {
    modification_names_.clear();
    for (vector<ResidueModification*>::iterator it = mods_.begin(); it != mods_.end(); ++it)
    {
      delete *it;
    }
  }


  Size ModificationsDB::getNumberOfModifications() const
  {
    return mods_.size();
  }


  const ResidueModification& ModificationsDB::getModification(Size index) const
  {
    if (index >= mods_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, mods_.size());
    }
    return *mods_[index];
  }


  void ModificationsDB::searchModifications(set<const ResidueModification*>& mods, const String& mod_name, const String& residue, ResidueModification::TermSpecificity term_spec) const
  {
    mods.clear();

    if (!modification_names_.has(mod_name))
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       mod_name);
    }

    const set<const ResidueModification*>& temp = modification_names_[mod_name];
    for (set<const ResidueModification*>::const_iterator it = temp.begin();
         it != temp.end(); ++it)
    {
      if (residuesMatch_(residue, (*it)->getOrigin()) &&
          (term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY ||
           (term_spec == (*it)->getTermSpecificity())))
      {
        mods.insert(*it);
      }
    }
  }


  const ResidueModification& ModificationsDB::getModification(const String& mod_name, const String& residue, ResidueModification::TermSpecificity term_spec) const
  {
    set<const ResidueModification*> mods;
    // if residue is specified, try residue-specific search first to avoid
    // ambiguities (e.g. "Carbamidomethyl (N-term)"/"Carbamidomethyl (C)"):
    if (!residue.empty() &&
        (term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY))
    {
      searchModifications(mods, mod_name, residue,
                          ResidueModification::ANYWHERE);
    }
    if (mods.empty()) searchModifications(mods, mod_name, residue, term_spec);

    if (mods.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Retrieving the modification failed. It is not available for the residue '" + String(residue) + "' and term specificity " + String(Int(term_spec)) + ".", mod_name);
    }
    if (mods.size() > 1)
    {
      LOG_WARN << "Warning (ModificationsDB::getModification): more than one modification with name '" + mod_name + "', residue '" + residue + "', specificity '" + String(Int(term_spec)) << "' found, picking the first one of:";
      for (set<const ResidueModification*>::const_iterator it = mods.begin();
           it != mods.end(); ++it)
      {
        LOG_WARN << " " << (*it)->getFullId();
      }
      LOG_WARN << "\n";
    }
    return **mods.begin();
  }


  bool ModificationsDB::has(String modification) const
  {
    OPENMS_PRECONDITION(!modification_names_.has(modification) || (int)findModificationIndex(modification) >= 0,
        "The modification being present implies that it can be found."); // NOTE: some very smart compilers may remove this statement ...
    return modification_names_.has(modification);
  }

  Size ModificationsDB::findModificationIndex(const String & mod_name) const
  {
    Int idx(-1);
    if (modification_names_.has(mod_name))
    {
      if (modification_names_[mod_name].size() > 1)
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "more than one element of name '" + mod_name + "' found!");
      }
    }
    else
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mod_name);
    }

    const ResidueModification* mod = *modification_names_[mod_name].begin();
    for (Size i = 0; i != mods_.size(); ++i)
    {
      if (mods_[i] == mod)
      {
        idx = i;
        break;
      }
    }

    // throw if we did not find the modification
    if (idx < 0)
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mod_name);
    }

    return idx;
  }


  void ModificationsDB::searchModificationsByDiffMonoMass(vector<String>& mods, double mass, double max_error, const String& residue, ResidueModification::TermSpecificity term_spec)
  {
    mods.clear();
    for (vector<ResidueModification*>::const_iterator it = mods_.begin();
         it != mods_.end(); ++it)
    {
      if ((fabs((*it)->getDiffMonoMass() - mass) <= max_error) &&
          residuesMatch_(residue, (*it)->getOrigin()) &&
          ((term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY) ||
           (term_spec == (*it)->getTermSpecificity())))
      {
        mods.push_back((*it)->getFullId());
      }
    }
  }


  const ResidueModification* ModificationsDB::getBestModificationByMonoMass(double mass, double max_error, const String& residue,
                                                                            ResidueModification::TermSpecificity term_spec)
  {
    double min_error = max_error;
    const ResidueModification* mod = nullptr;
    const Residue* residue_ = ResidueDB::getInstance()->getResidue(residue); // is NULL if not found
    for (vector<ResidueModification*>::const_iterator it = mods_.begin();
         it != mods_.end(); ++it)
    {
      double mono_mass = (*it)->getMonoMass();
      if ((mono_mass <= 0) && !residue.empty())
      {
        // Since not all modifications have a monoisotopic mass stored (they may
        // map to multiple residues), we calculate a monoisotopic mass from the
        // delta mass.
        // First the internal (inside an AA chain) weight of the residue:
        if (residue_ != nullptr) continue; // @TODO: throw an exception here?
        double internal_weight = residue_->getMonoWeight() -
          residue_->getInternalToFull().getMonoWeight();
        mono_mass = (*it)->getDiffMonoMass() + internal_weight;
      }
      // using less instead of less-or-equal will pick the first matching
      // modification of equally heavy modifications (in our case this is the
      // first matching UniMod entry)
      double mass_error = fabs(mono_mass - mass);
      if ((mass_error < min_error) &&
          residuesMatch_(residue, (*it)->getOrigin()) &&
          ((term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY) ||
           (term_spec == (*it)->getTermSpecificity())))
      {
        min_error = mass_error;
        mod = *it;
      }
    }
    return mod;
  }


  const ResidueModification* ModificationsDB::getBestModificationByDiffMonoMass(double mass, double max_error, const String& residue, ResidueModification::TermSpecificity term_spec)
  {
    double min_error = max_error;
    const ResidueModification* mod = nullptr;
    for (vector<ResidueModification*>::const_iterator it = mods_.begin();
         it != mods_.end(); ++it)
    {
      // using less instead of less-or-equal will pick the first matching
      // modification of equally heavy modifications (in our case this is the
      // first matching UniMod entry)
      double mass_error = fabs((*it)->getDiffMonoMass() - mass);
      if ((mass_error < min_error) &&
          residuesMatch_(residue, (*it)->getOrigin()) &&
          ((term_spec == ResidueModification::NUMBER_OF_TERM_SPECIFICITY) ||
           (term_spec == (*it)->getTermSpecificity())))
      {
        min_error = mass_error;
        mod = *it;
      }
    }
    return mod;
  }

  void ModificationsDB::readFromUnimodXMLFile(const String& filename)
  {
    vector<ResidueModification*> new_mods;
    UnimodXMLFile().load(filename, new_mods);

    for (vector<ResidueModification*>::iterator it = new_mods.begin(); it != new_mods.end(); ++it)
    {
      // create full ID based on other information:
      (*it)->setFullId();

      // e.g. Oxidation (M)
      modification_names_[(*it)->getFullId()].insert(*it);
      // e.g. Oxidation
      modification_names_[(*it)->getId()].insert(*it);
      // e.g. Oxidized
      modification_names_[(*it)->getFullName()].insert(*it);
      // e.g. UniMod:312
      modification_names_[(*it)->getUniModAccession()].insert(*it);
      mods_.push_back(*it);
    }
  }

  void ModificationsDB::addModification(ResidueModification* new_mod)
  {
    if (has(new_mod->getFullId()))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Modification already exists in ModificationsDB.", String(new_mod->getFullId()));
    }
    modification_names_[new_mod->getFullId()].insert(new_mod);
    modification_names_[new_mod->getId()].insert(new_mod);
    modification_names_[new_mod->getFullName()].insert(new_mod);
    modification_names_[new_mod->getUniModAccession()].insert(new_mod);

    mods_.push_back(new_mod); // we probably want that
  }

  void ModificationsDB::readFromOBOFile(const String& filename)
  {
    ResidueModification mod;
    // add multiple mods for multiple specificities
    //Map<String, ResidueModification> all_mods;
    multimap<String, ResidueModification> all_mods;

    ifstream is(File::find(filename).c_str());
    String line, line_wo_spaces, id;
    String origin = "";

    bool reading_cross_link = false;

    //parse file
    while (getline(is, line, '\n'))
    {
      line.trim();
      line_wo_spaces = line;
      line_wo_spaces.removeWhitespaces();

      if (line == "" || line[0] == '!') //skip empty lines and comments
      {
        continue;
      }

      if (line_wo_spaces == "[Term]")       //new term
      {
        // if the last [Term] was a moon-link, then it does not belong in CrossLinksDB
        if (id != "" && !reading_cross_link) //store last term
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

          // for mono-links from XLMOD.obo:
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
        else if (reading_cross_link) // re-initialize before reading next [Term]
        {
          id = "";
          origin = "";
          mod = ResidueModification();
          reading_cross_link = false;
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
          if (val_split[1] == "2")
          {
            reading_cross_link = true;
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

    if (id != "") //store last term
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

      // for mono-links from XLMOD.obo:
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

  void ModificationsDB::getAllSearchModifications(vector<String>& modifications) const
  {
    modifications.clear();

    for (vector<ResidueModification*>::const_iterator it = mods_.begin(); it != mods_.end(); ++it)
    {
      if ((*it)->getUniModRecordId() > 0)
      {
        modifications.push_back((*it)->getFullId());
      }
    }
    sort(modifications.begin(), modifications.end());
  }


  bool ModificationsDB::residuesMatch_(const String& residue, char origin) const
  {
    return (residue.empty() || (origin == residue[0]) || (residue == "X") || (origin == 'X') || (residue == "."));
  }

} // namespace OpenMS
