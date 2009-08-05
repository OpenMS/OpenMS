// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/DENOVO/MassDecompositionAlgorithm.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationDefinition.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <set>

using namespace std;

namespace OpenMS
{
	MassDecompositionAlgorithm::MassDecompositionAlgorithm()
	: DefaultParamHandler("MassDecompositionAlgorithm"),
		alphabet_(0),
		decomposer_(0)
	{
		defaults_.setValue("decomp_weights_precision", 0.01, "precision used to calculate the decompositions, this only affects cache usage!");
		defaults_.setValue("tolerance", 0.3, "tolerance which is allowed for the decompositions");
		defaults_.setValue("fixed_modifications", "", "fixed modifications, specified using PSI-MOD terms, e.g. MOD:01214,MOD:00048");
		defaults_.setValue("variable_modifications", "MOD:00719,MOD:01329", "variable modifications, specified using PSI-MOD terms, e.g. MOD:01214,MOD:00048");
		defaultsToParam_();
	}

	MassDecompositionAlgorithm::~MassDecompositionAlgorithm()
	{
		delete alphabet_;
		delete decomposer_;
	}

	void MassDecompositionAlgorithm::getDecompositions(vector<MassDecomposition>& decomps, double mass)
	{
		double tolerance((double)param_.getValue("tolerance"));
    ims::RealMassDecomposer::decompositions_type decompositions = decomposer_->getDecompositions(mass, tolerance);

    for (ims::RealMassDecomposer::decompositions_type::const_iterator pos = decompositions.begin(); pos != decompositions.end(); ++pos)
    {
      String d;
      for (ims::Alphabet::size_type i = 0; i < alphabet_->size(); ++i)
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
    Map<char, double> aa_to_weight;

    aa_to_weight['K'] = 128.095;
    aa_to_weight['M'] = 131.04;
    aa_to_weight['F'] = 147.068;
    aa_to_weight['P'] = 97.0528;
    aa_to_weight['S'] = 87.032;
    aa_to_weight['T'] = 101.048;
    aa_to_weight['W'] = 186.079;
    aa_to_weight['Y'] = 163.063;
    aa_to_weight['V'] = 99.0684;
    aa_to_weight['A'] = 71.0371;
    aa_to_weight['R'] = 156.101;
    aa_to_weight['N'] = 114.043;
    aa_to_weight['D'] = 115.027;
    aa_to_weight['C'] = 103.00919;
    aa_to_weight['E'] = 129.043;
    aa_to_weight['Q'] = 128.059;
    aa_to_weight['G'] = 57.0215;
    aa_to_weight['H'] = 137.059;
    aa_to_weight['L'] = 113.084;

		// now handle the modifications
		ModificationDefinitionsSet mod_set((String)param_.getValue("fixed_modifications"), (String)param_.getValue("variable_modifications"));
		set<ModificationDefinition> fixed_mods = mod_set.getFixedModifications();
    for (set<ModificationDefinition>::const_iterator it = fixed_mods.begin(); it != fixed_mods.end(); ++it)
    {
      ResidueModification mod = ModificationsDB::getInstance()->getModification(it->getModification());
      char aa=' ';
      if (mod.getOrigin().size() != 1 || mod.getOrigin() == "X")
      {
        cerr << "Warning: cannot handle modification " << it->getModification() << ", because aa is ambiguous (" << mod.getOrigin() << "), ignoring modification!" << endl;
        continue;
      }
      else
      {
        aa = mod.getOrigin()[0];
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
          cerr << "Warning: cannot handle modification " << it->getModification() << ", because no monoisotopic mass value was found! Ignoringmodification!" << endl;
          continue;
        }
      }
    }

		const StringList mod_names(StringList::create("a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z"));
    vector<String>::const_iterator actual_mod_name = mod_names.begin();
    set<ModificationDefinition> var_mods = mod_set.getVariableModifications();
    for (set<ModificationDefinition>::const_iterator it = var_mods.begin(); it != var_mods.end(); ++it)
    {
      ResidueModification mod = ModificationsDB::getInstance()->getModification(it->getModification());
      char aa = (*actual_mod_name)[0];
      char origin_aa = ' ';
      ++actual_mod_name;

      if (mod.getOrigin().size() != 1 || mod.getOrigin() == "X")
      {
        cerr << "Warning: cannot handle modification " << it->getModification() << ", because aa is ambiguous (" << mod.getOrigin() << "), ignoring modification!" << endl;
        continue;
      }
      else
      {
        origin_aa = mod.getOrigin()[0];
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
          cerr << "Warning: cannot handle modification " << it->getModification() << ", because no monoisotopic mass value was found! Ignoringmodification!" << endl;
          continue;
        }
      }
    }

		if (alphabet_ != 0)
		{
			delete alphabet_;
		}
		if (decomposer_ != 0)
		{
			delete decomposer_;
		}

    // init mass decomposer
    alphabet_ = new ims::Alphabet();
		for (Map<char, double>::ConstIterator it = aa_to_weight.begin(); it != aa_to_weight.end(); ++it)
		{
			alphabet_->push_back(String(it->first), it->second);
		}

    // initializes weights
    ims::Weights weights(alphabet_->getMasses(), (double)param_.getValue("decomp_weights_precision"));

    // optimize alphabet by dividing by gcd
    weights.divideByGCD();

    // decomposes real values
    decomposer_ = new ims::RealMassDecomposer(weights);

		return;
	}

	

} // namespace OpenMS

