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

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
  ResidueDB::ResidueDB()
  {
    readResiduesFromFile_("CHEMISTRY/Residues.xml");
    buildResidueNames_();
  }

  ResidueDB::~ResidueDB()
  {
    clear_();
  }

  const Residue* ResidueDB::getResidue(const String& name) const
  {
    if (residue_names_.find(name) != residue_names_.end())
    {
      return residue_names_.at(name);
    }
    return nullptr;
  }

  const Residue* ResidueDB::getResidue(const unsigned char& one_letter_code) const
  {
    return residue_by_one_letter_code_[one_letter_code];
  }

  Size ResidueDB::getNumberOfResidues() const
  {
    return residues_.size();
  }

  Size ResidueDB::getNumberOfModifiedResidues() const
  {
    return modified_residues_.size();
  }

  const set<const Residue*> ResidueDB::getResidues(const String& residue_set) const
  {
    if (!residues_by_set_.has(residue_set))
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Residue set cannot be found: '" + residue_set + "'");
    }

    return residues_by_set_[residue_set];
  }

  void ResidueDB::setResidues(const String& file_name)
  {
    clearResidues_();
    readResiduesFromFile_(file_name);
    buildResidueNames_();
  }

  void ResidueDB::addResidue(const Residue& residue)
  {
    Residue* r = new Residue(residue);
    addResidue_(r);
  }

  void ResidueDB::addResidue_(Residue* r)
  {
    vector<String> names;
    if (r->getName() != "")
    {
      names.push_back(r->getName());
    }
    if (r->getShortName() != "")
    {
      names.push_back(r->getShortName());
    }
    set<String> synonyms = r->getSynonyms();
    for (set<String>::iterator it = synonyms.begin(); it != synonyms.end(); ++it)
    {
      names.push_back(*it);
    }

    if (!r->isModified())
    {
      for (vector<String>::const_iterator it = names.begin(); it != names.end(); ++it)
      {
        residue_names_[*it] = r;
        residue_by_one_letter_code_[(unsigned char)(*it)[0]] = r;
      }
      residues_.insert(r);
      const_residues_.insert(r);
    }
    else
    {
      modified_residues_.insert(r);
      const_modified_residues_.insert(r);

      // get all modification names
      vector<String> mod_names;
      const ResidueModification* mod = r->getModification();

      mod_names.push_back(mod->getId());
      mod_names.push_back(mod->getFullName());
      mod_names.push_back(mod->getFullId());
      const set<String>& mod_synonyms = mod->getSynonyms();
      for (set<String>::const_iterator it = mod_synonyms.begin(); it != mod_synonyms.end(); ++it)
      {
        mod_names.push_back(*it);
      }

      for (vector<String>::const_iterator it = names.begin(); it != names.end(); ++it)
      {
        for (vector<String>::const_iterator mod_it = mod_names.begin(); mod_it != mod_names.end(); ++mod_it)
        {
          if ( mod_it->empty() || it->empty() ) continue;
          residue_mod_names_[*it][*mod_it] = r;
        }
      }
    }
    buildResidueNames_();
    return;
  }

  bool ResidueDB::hasResidue(const String& res_name) const
  {
    if (residue_names_.find(res_name) != residue_names_.end())
    {
      return true;
    }
    return false;
  }

  bool ResidueDB::hasResidue(const Residue* residue) const
  {
    if (const_residues_.find(residue) != const_residues_.end() ||
        const_modified_residues_.find(residue) != const_modified_residues_.end())
    {
      return true;
    }
    return false;
  }

  void ResidueDB::readResiduesFromFile_(const String& file_name)
  {
    String file = File::find(file_name);

    Param param;
    ParamXMLFile paramFile;
    paramFile.load(file, param);

    if (!param.begin().getName().hasPrefix("Residues"))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "");
    }

    try
    {
      vector<String> split;
      param.begin().getName().split(':', split);
      String prefix = split[0] + split[1];
      Residue* res_ptr = nullptr;

      Map<String, String> values;

      for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
      {
        it.getName().split(':', split);
        if (prefix != split[0] + split[1])
        {
          // add residue
          res_ptr = parseResidue_(values);
          values.clear();
          residues_.insert(res_ptr);
          const_residues_.insert(res_ptr);
          prefix = split[0] + split[1];
        }

        String value = it->value;
        String key = it.getName();
        values[key] = value;

      }

      // add last residue
      res_ptr = parseResidue_(values);
      residues_.insert(res_ptr);
      const_residues_.insert(res_ptr);
    }
    catch (Exception::BaseException& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, e.what(), "");
    }
  }

  void ResidueDB::clear_()
  {
    clearResidues_();
    //clearResidueModifications_();
  }

  void ResidueDB::clearResidues_()
  {
    set<Residue*>::iterator it;
    for (it = residues_.begin(); it != residues_.end(); ++it)
    {
      delete *it;
    }

    residues_.clear();
    residue_names_.clear();
    const_residues_.clear();
  }

  Residue* ResidueDB::parseResidue_(Map<String, String>& values)
  {
    vector<EmpiricalFormula> low_mass_ions;
    Residue* res_ptr = new Residue();

    for (Map<String, String>::iterator it = values.begin(); it != values.end(); ++it)
    {
      String key(it->first);
      String value(it->second);

      if (key.hasSuffix(":Name"))
      {
        res_ptr->setName(value);
        continue;
      }
      if (key.hasSuffix(":ShortName"))
      {
        res_ptr->setShortName(value);
        continue;
      }
      if (key.hasSuffix(":ThreeLetterCode"))
      {
        res_ptr->setThreeLetterCode(value);
        continue;
      }
      if (key.hasSuffix(":OneLetterCode"))
      {
        res_ptr->setOneLetterCode(value);
        continue;
      }
      if (key.hasSuffix(":Formula"))
      {
        EmpiricalFormula formula(value);
        res_ptr->setFormula(EmpiricalFormula(value));
        res_ptr->setAverageWeight(formula.getAverageWeight());
        res_ptr->setMonoWeight(formula.getMonoWeight());
        continue;
      }

      if (key.hasSubstring(":Losses:LossName"))
      {
        res_ptr->addLossName(value);
        continue;
      }
      if (key.hasSubstring(":Losses:LossFormula"))
      {
        EmpiricalFormula loss(value);
        res_ptr->addLossFormula(loss);
        continue;
      }

      if (key.hasSubstring("NTermLosses:LossName"))
      {
        res_ptr->addNTermLossName(value);
        continue;
      }

      if (key.hasSubstring("NTermLosses:LossFormula"))
      {
        EmpiricalFormula loss(value);
        res_ptr->addNTermLossFormula(loss);
        continue;
      }

      if (key.hasSubstring("LowMassIons"))
      {
        // no markers defined?
        if (!key.hasSuffix(":"))
        {
          low_mass_ions.push_back(EmpiricalFormula(value));
        }
        continue;
      }
      if (key.hasSubstring("Synonyms"))
      {
        // no synonyms defined?
        if (!key.hasSuffix(":"))
        {
          res_ptr->addSynonym(value);
        }
        continue;
      }
      if (key.hasSubstring("pka"))
      {
        // no pka defined?
        if (!key.hasSuffix(":"))
        {
          res_ptr->setPka(value.toDouble());
        }
        continue;
      }
      if (key.hasSubstring("pkb"))
      {
        // no pkb defined?
        if (!key.hasSuffix(":"))
        {
          res_ptr->setPkb(value.toDouble());
        }
        continue;
      }
      if (key.hasSubstring("pkc"))
      {
        // no pkc defined?
        if (!key.hasSuffix(":"))
        {
          res_ptr->setPkc(value.toDouble());
        }
        continue;
      }
      if (key.hasSubstring("GB_SC"))
      {
        res_ptr->setSideChainBasicity(value.toDouble());
        continue;
      }
      if (key.hasSubstring("GB_BB_L"))
      {
        res_ptr->setBackboneBasicityLeft(value.toDouble());
        continue;
      }
      if (key.hasSubstring("GB_BB_R"))
      {
        res_ptr->setBackboneBasicityRight(value.toDouble());
        continue;
      }
      if (key.hasSubstring("ResidueSets"))
      {
        StringList residue_sets = ListUtils::create<String>(value);
        for (StringList::const_iterator local_it = residue_sets.begin(); local_it != residue_sets.end(); ++local_it)
        {
          res_ptr->addResidueSet(*local_it);
          residue_sets_.insert(*local_it);
        }
        continue;
      }
      cerr << "unknown key: " << key << ", with value: " << value << endl;
    }

    if (!low_mass_ions.empty())
    {
      res_ptr->setLowMassIons(low_mass_ions);
    }

    for (set<String>::const_iterator it = res_ptr->getResidueSets().begin(); it != res_ptr->getResidueSets().end(); ++it)
    {
      residues_by_set_[*it].insert(res_ptr);
    }

    return res_ptr;
  }

  const set<String>& ResidueDB::getResidueSets() const
  {
    return residue_sets_;
  }

  void ResidueDB::buildResidueNames_()
  {
    // initialize lookup table to null pointer
    for (Size i = 0; i != sizeof(residue_by_one_letter_code_)/sizeof(residue_by_one_letter_code_[0]); ++i)
    {
      residue_by_one_letter_code_[i] = nullptr;
    }

    set<Residue*>::iterator it;
    for (it = residues_.begin(); it != residues_.end(); ++it)
    {
      residue_names_[(*it)->getName()] = *it;
      if ((*it)->getThreeLetterCode() != "")
      {
        residue_names_[(*it)->getThreeLetterCode()] = *it;
      }
      if ((*it)->getOneLetterCode() != "")
      {
        residue_names_[(*it)->getOneLetterCode()] = *it;
        const unsigned char l = (*it)->getOneLetterCode()[0];
        residue_by_one_letter_code_[l] = *it;
      }
      if ((*it)->getShortName() != "")
      {
        residue_names_[(*it)->getShortName()] = *it;
      }
      set<String>::iterator sit;
      set<String> syn = (*it)->getSynonyms();
      for (sit = syn.begin(); sit != syn.end(); ++sit)
      {
        if (*sit != "")
        {
          residue_names_[*sit] = *it;
        }
      }
    }
  }

  const Residue* ResidueDB::getModifiedResidue(const String& modification)
  {
    // terminal mods. don't apply to residue (side chain), so don't consider them:
    const ResidueModification& mod = ModificationsDB::getInstance()->getModification(modification, "", ResidueModification::ANYWHERE);
    return getModifiedResidue(getResidue(mod.getOrigin()), mod.getFullId());
  }

  const Residue* ResidueDB::getModifiedResidue(const Residue* residue, const String& modification)
  {
    OPENMS_PRECONDITION(!modification.empty(), "Modification cannot be empty")
    // search if the mod already exists
    String res_name = residue->getName();

    if (residue_names_.find(res_name) == residue_names_.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       String("Residue with name " + res_name + " was not registered in residue DB, register first!").c_str());
    }

    // terminal mods. don't apply to residue (side chain), so don't consider them:
    const ResidueModification& mod = ModificationsDB::getInstance()->getModification(modification, residue->getOneLetterCode(), ResidueModification::ANYWHERE);
    String id = mod.getId();
    if (id.empty()) id = mod.getFullId();

    if (residue_mod_names_.has(res_name) && residue_mod_names_[res_name].has(id))
    {
      return residue_mod_names_[res_name][id];
    }

    Residue* res = new Residue(*residue_names_[res_name]);
    res->setModification_(mod);
    //res->setLossFormulas(vector<EmpiricalFormula>());
    //res->setLossNames(vector<String>());

    // now register this modified residue
    addResidue_(res);
    return res;
  }

}
