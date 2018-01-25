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

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationDefinition.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>

#include <set>
#include <iostream>

using namespace std;

namespace OpenMS
{
  MassDecompositionAlgorithm::MassDecompositionAlgorithm() :
    DefaultParamHandler("MassDecompositionAlgorithm"),
    alphabet_(nullptr),
    decomposer_(nullptr)
  {
    defaults_.setValue("decomp_weights_precision", 0.01, "precision used to calculate the decompositions, this only affects cache usage!", ListUtils::create<String>("advanced"));
    defaults_.setValue("tolerance", 0.3, "tolerance which is allowed for the decompositions");

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    defaults_.setValue("fixed_modifications", ListUtils::create<String>(""), "fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'");
    defaults_.setValidStrings("fixed_modifications", all_mods);
    defaults_.setValue("variable_modifications", ListUtils::create<String>(""), "variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'");
    defaults_.setValidStrings("variable_modifications", all_mods);
    defaults_.setValue("residue_set", "Natural19WithoutI", "The predefined amino acid set that should be used, see doc of ResidueDB for possible residue sets", ListUtils::create<String>("advanced"));
    set<String> residue_sets = ResidueDB::getInstance()->getResidueSets();
    vector<String> valid_strings;
    for (set<String>::const_iterator it = residue_sets.begin(); it != residue_sets.end(); ++it)
    {
      valid_strings.push_back(*it);
    }
    defaults_.setValidStrings("residue_set", valid_strings);
    defaultsToParam_();
  }

  MassDecompositionAlgorithm::~MassDecompositionAlgorithm()
  {
    delete alphabet_;
    delete decomposer_;
  }

  void MassDecompositionAlgorithm::getDecompositions(vector<MassDecomposition> & decomps, double mass)
  {
    double tolerance((double) param_.getValue("tolerance"));
    ims::RealMassDecomposer::decompositions_type decompositions = decomposer_->getDecompositions(mass, tolerance);

    for (ims::RealMassDecomposer::decompositions_type::const_iterator pos = decompositions.begin(); pos != decompositions.end(); ++pos)
    {
      String d;
      for (ims::IMSAlphabet::size_type i = 0; i < alphabet_->size(); ++i)
      {
        if ((*pos)[i] > 0)
        {
          d += alphabet_->getName(i) + String((*pos)[i]) + " ";
        }
      }
      d.trim();
      MassDecomposition decomp(d);
      decomps.push_back(decomp);
    }

    return;
  }

  void MassDecompositionAlgorithm::updateMembers_()
  {
    // todo add accessor to tolerance, it is called very often in CID mode

    Map<char, double> aa_to_weight;

    set<const Residue *> residues = ResidueDB::getInstance()->getResidues((String)param_.getValue("residue_set"));

    for (set<const Residue *>::const_iterator it = residues.begin(); it != residues.end(); ++it)
    {
      aa_to_weight[(*it)->getOneLetterCode()[0]] = (*it)->getMonoWeight(Residue::Internal);
    }

    // now handle the modifications
    ModificationDefinitionsSet mod_set(param_.getValue("fixed_modifications"), param_.getValue("variable_modifications"));
    set<ModificationDefinition> fixed_mods = mod_set.getFixedModifications();
    for (set<ModificationDefinition>::const_iterator it = fixed_mods.begin(); it != fixed_mods.end(); ++it)
    {
      const ResidueModification& mod = it->getModification();
      char aa = ' ';
      if (mod.getOrigin() == 'X')
      {
        cerr << "MassDecompositionAlgorithm: Warning: cannot handle modification " << mod.getName() << ", because aa is ambiguous (" << mod.getOrigin() << "), ignoring modification!" << endl;
        continue;
      }
      else
      {
        aa = mod.getOrigin();
      }

      if (mod.getMonoMass() != 0)
      {
        aa_to_weight[aa] = mod.getMonoMass();
      }
      else
      {
        if (mod.getDiffMonoMass() != 0)
        {
          aa_to_weight[aa] += mod.getDiffMonoMass();
        }
        else
        {
          cerr << "MassDecompositionAlgorithm: Warning: cannot handle modification " << mod.getName() << ", because no monoisotopic mass value was found! Ignoring modification!" << endl;
          continue;
        }
      }
    }

    const StringList mod_names(ListUtils::create<String>("a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z"));
    vector<String>::const_iterator actual_mod_name = mod_names.begin();
    set<ModificationDefinition> var_mods = mod_set.getVariableModifications();
    for (set<ModificationDefinition>::const_iterator it = var_mods.begin(); it != var_mods.end(); ++it)
    {
      ResidueModification mod = it->getModification();
      //cerr << it->getModification() << " " << mod.getOrigin() << " " << mod.getId() << " " << mod.getFullId() << " " << mod.getUniModAccession() << " " << mod.getPSIMODAccession() << endl;
      char aa = (*actual_mod_name)[0];
      char origin_aa = ' ';
      ++actual_mod_name;

      if (mod.getOrigin() == 'X')
      {
        cerr << "MassDecompositionAlgorithm: Warning: cannot handle modification " << mod.getName() << ", because aa is ambiguous (" << mod.getOrigin() << "), ignoring modification!" << endl;
        continue;
      }
      else
      {
        origin_aa = mod.getOrigin();
      }

      if (mod.getMonoMass() != 0)
      {
        aa_to_weight[aa] = mod.getMonoMass();
      }
      else
      {
        if (mod.getDiffMonoMass() != 0)
        {
          aa_to_weight[aa] = aa_to_weight[origin_aa] + mod.getDiffMonoMass();
        }
        else
        {
          cerr << "Warning: cannot handle modification " << mod.getName() << ", because no monoisotopic mass value was found! Ignoring modification!" << endl;
          continue;
        }
      }
    }

    if (alphabet_ != nullptr)
    {
      delete alphabet_;
    }
    if (decomposer_ != nullptr)
    {
      delete decomposer_;
    }

    // init mass decomposer
    alphabet_ = new ims::IMSAlphabet();
    for (Map<char, double>::ConstIterator it = aa_to_weight.begin(); it != aa_to_weight.end(); ++it)
    {
      alphabet_->push_back(String(it->first), it->second);
    }

    // initializes weights
    ims::Weights weights(alphabet_->getMasses(), (double) param_.getValue("decomp_weights_precision"));

    // optimize alphabet by dividing by gcd
    weights.divideByGCD();

    // decomposes real values
    decomposer_ = new ims::RealMassDecomposer(weights);

    return;
  }

} // namespace OpenMS
