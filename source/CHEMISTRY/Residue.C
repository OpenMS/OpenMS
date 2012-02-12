// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
//

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <cstdlib>

using namespace std;

namespace OpenMS
{
	// residue
	Residue::Residue()
		: name_("unknown"),
			average_weight_(0.0f),
			mono_weight_(0.0f),
			is_modified_(false),
			modification_(""),
			loss_average_weight_(0.0f),
			loss_mono_weight_(0.0f),
			pka_(0.0),
			pkb_(0.0),
			pkc_(-1.0)
	{
	}

	Residue::Residue(	const String& name,
										const String& three_letter_code,
										const String& one_letter_code,
										const EmpiricalFormula& formula)
		:	name_(name),
			three_letter_code_(three_letter_code),
			one_letter_code_(one_letter_code),
			formula_(formula),
			average_weight_(0),
			mono_weight_(0),
			is_modified_(false),
			modification_(""),
			loss_average_weight_(0.0f),
			loss_mono_weight_(0.0f),
			pka_(0.0),
			pkb_(0.0),
			pkc_(-1.0),
			gb_sc_(0.0),
			gb_bb_l_(0.0),
			gb_bb_r_(0.0)
	{
		if (formula_ != "")
		{
			internal_formula_ = formula_ - getInternalToFull();
		}
	}

	Residue::Residue(const Residue& residue)
		: name_(residue.name_),
			short_name_(residue.short_name_),
			synonyms_(residue.synonyms_),
			three_letter_code_(residue.three_letter_code_),
			one_letter_code_(residue.one_letter_code_),
			formula_(residue.formula_),
			internal_formula_(residue.internal_formula_),
			average_weight_(residue.average_weight_),
			mono_weight_(residue.mono_weight_),
			is_modified_(residue.is_modified_),
			pre_mod_name_(residue.pre_mod_name_),
			modification_(residue.modification_),
			loss_names_(residue.loss_names_),
			loss_formulas_(residue.loss_formulas_),
			NTerm_loss_names_(residue.NTerm_loss_names_),
			NTerm_loss_formulas_(residue.NTerm_loss_formulas_),
			loss_average_weight_(residue.loss_average_weight_),
			loss_mono_weight_(residue.loss_mono_weight_),
			low_mass_ions_(residue.low_mass_ions_),
			pka_(residue.pka_),
			pkb_(residue.pkb_),
			pkc_(residue.pkc_),
			gb_sc_(residue.gb_sc_),
			gb_bb_l_(residue.gb_bb_l_),
			gb_bb_r_(residue.gb_bb_r_),
			residue_sets_(residue.residue_sets_)
	{
	}
	
	Residue::~Residue()
	{
	}

	Residue& Residue::operator = (const Residue& residue)
	{
		if (this != &residue)
		{
			name_ = residue.name_;
			short_name_ = residue.short_name_;
			synonyms_ = residue.synonyms_;
			three_letter_code_ = residue.three_letter_code_;
			one_letter_code_ = residue.one_letter_code_;
			formula_ = residue.formula_;
			internal_formula_ = residue.internal_formula_;
			average_weight_ = residue.average_weight_;
			mono_weight_ = residue.mono_weight_;
			is_modified_ = residue.is_modified_;
			pre_mod_name_ = residue.pre_mod_name_;
			modification_ = residue.modification_;
			loss_names_ = residue.loss_names_;
			loss_formulas_ = residue.loss_formulas_;
			NTerm_loss_names_ = residue.NTerm_loss_names_;
			NTerm_loss_formulas_ = residue.NTerm_loss_formulas_;
			loss_average_weight_ = residue.loss_average_weight_;
			loss_mono_weight_ = residue.loss_mono_weight_;
			low_mass_ions_ = residue.low_mass_ions_;
			pka_ = residue.pka_;
			pkb_ = residue.pkb_;
			pkc_ = residue.pkc_;
			gb_sc_ = residue.gb_sc_;
			gb_bb_l_ = residue.gb_bb_l_;
			gb_bb_r_ = residue.gb_bb_r_;
			residue_sets_ = residue.residue_sets_;
		}
		return *this;
	}

	void Residue::setName(const String& name)
	{
		name_ = name;
	}

	const String& Residue::getName() const
	{
		return name_;
	}

	String Residue::getResidueTypeName(const Residue::ResidueType res_type)
	{
		String ion("-ion");
		switch (res_type)
		{
  		case Residue::Full:
	  		return "full";
		  case Residue::Internal:
			  return "internal";
		  case Residue::NTerminal:
			  return "N-terminal";
		  case Residue::CTerminal:
			  return "C-terminal";
			case Residue::AIon:
				return "a" + ion;
			case Residue::BIon:
				return "b" + ion;
			case Residue::CIon:
				return "c" + ion;
			case Residue::CIonMinusOne:
				return "c-1" + ion;
			case Residue::CIonPlusOne:
				return "c+1" + ion;
			case Residue::CIonPlusTwo:
				return "c+2" + ion;
			case Residue::XIon:
				return "x" + ion;
			case Residue::YIon:
				return "y" + ion;
			case Residue::ZIon:
				return "z" + ion;
			case Residue::ZIonMinusOne:
				return "z-1" + ion;
			case Residue::ZIonPlusOne:
				return "z+1" + ion;
			case Residue::ZIonPlusTwo:
				return "z+2" + ion;
			default:
				cerr << "Residue::getResidueTypeName: residue type has no name" << endl;
		}
		return "";
	}

	void Residue::setShortName(const String& short_name)
	{
		short_name_ = short_name;
	}

	const String& Residue::getShortName() const
	{
		return short_name_;
	}

	void Residue::setSynonyms(const set<String>& synonyms)
	{
		synonyms_ = synonyms;
	}

	void Residue::addSynonym(const String& synonym)
	{
		synonyms_.insert(synonym);
	}
	
	const set<String>& Residue::getSynonyms() const
	{
		return synonyms_;
	}
	
	void Residue::setThreeLetterCode(const String& three_letter_code)
	{
		three_letter_code_ = three_letter_code;
	}

	const String& Residue::getThreeLetterCode() const
	{
		return three_letter_code_;
	}

	void Residue::setOneLetterCode(const String& one_letter_code)
	{
		one_letter_code_ = one_letter_code;
	}

	const String& Residue::getOneLetterCode() const
	{
		return one_letter_code_;
	}

	DoubleReal Residue::getPka() const
	{
		return pka_;
	}

	DoubleReal Residue::getPkb() const
	{
		return pkb_;
	}

	DoubleReal Residue::getPkc() const
	{
		return pkc_;
	}

	DoubleReal Residue::getPiValue() const
	{
		DoubleReal temp1 = 0.0;
		DoubleReal temp2 = 0.0;
		DoubleReal temp3 = 0.0;
		DoubleReal pi = 0;
		
		temp1 = getPka();
		temp2 = getPkb();
		temp3 = getPkc();
				
		if (temp3 >= 0 && temp3 < temp1)
		{			
			pi = (temp3 + temp2) / 2;
		}
		else if (temp3 >= temp2)
		{			
			pi = (temp1 + temp3) / 2;
		}
		else
		{
			pi = (temp1 + temp2) / 2;
		}
		
		return pi;
	}

	void Residue::setPka(DoubleReal value)
	{
		pka_ = value;
	}

	void Residue::setPkb(DoubleReal value)
	{
		pkb_ = value;
	}

	void Residue::setPkc(DoubleReal value)
	{
		pkc_ = value;
	}

	void Residue::setLossFormulas(const vector<EmpiricalFormula>& loss_formulas)
	{
		loss_formulas_ = loss_formulas;
	}

	void Residue::addLossFormula(const EmpiricalFormula& loss_formula)
	{
		loss_formulas_.push_back(loss_formula);
	}
	
	const vector<EmpiricalFormula>& Residue::getLossFormulas() const
	{
		return loss_formulas_;
	}

	void Residue::addLossName(const String& name)
	{
		loss_names_.push_back(name);
	}
	
	void Residue::setLossNames(const vector<String>& names)
	{
		loss_names_ = names;
	}

	const vector<String>& Residue::getLossNames() const
	{
		return loss_names_;
	}


  void Residue::setNTermLossFormulas(const vector<EmpiricalFormula>& NTerm_loss_formulas)
  {
    NTerm_loss_formulas_ = NTerm_loss_formulas;
  }

  void Residue::addNTermLossFormula(const EmpiricalFormula& NTerm_loss_formula)
  {
    NTerm_loss_formulas_.push_back(NTerm_loss_formula);
  }

  const vector<EmpiricalFormula>& Residue::getNTermLossFormulas() const
  {
    return NTerm_loss_formulas_;
  }

  void Residue::addNTermLossName(const String& name)
  {
    NTerm_loss_names_.push_back(name);
  }

  void Residue::setNTermLossNames(const vector<String>& names)
  {
    NTerm_loss_names_ = names;
  }

  const vector<String>& Residue::getNTermLossNames() const
  {
    return NTerm_loss_names_;
  }
	
	
	void Residue::setFormula(const EmpiricalFormula& formula)
	{
		formula_ = formula; 
		internal_formula_ = formula_ - getInternalToFull();
	}

	EmpiricalFormula Residue::getFormula(ResidueType res_type) const
	{
		switch (res_type)
		{
			case Full:
				return formula_;
			case Internal:
				return internal_formula_;
			case NTerminal:
				return formula_ - getNTerminalToFull();
			case CTerminal:
				return formula_ - getCTerminalToFull();
			case BIon:
				return formula_ - getBIonToFull();
			case AIon:
				return formula_ - getAIonToFull();
			case CIonMinusOne:
				return formula_ - getCIonMinusOneToFull();
			case CIon:
				return formula_ - EmpiricalFormula("OH") + EmpiricalFormula("NH");
			case XIon:
				return formula_ + getXIonToFull();
			case YIon:
				return formula_ + getYIonToFull();
			case ZIonMinusOne:
					return formula_ - getZIonMinusOneToFull();
			case ZIon:
				return formula_ - getZIonToFull();
			case ZIonPlusOne:
				return formula_ - getZIonPlusOneToFull();
			case ZIonPlusTwo:
				return formula_ - getZIonPlusTwoToFull();
			default:
				cerr << "Residue::getFormula: unknown ResidueType" << endl;
				return formula_;
		}
	}

	void Residue::setAverageWeight(DoubleReal weight) 
	{
		average_weight_ = weight;
		return;
	}

	DoubleReal Residue::getAverageWeight(ResidueType res_type) const
	{
		switch (res_type)
		{
			case Full:
				return average_weight_;
			case Internal:
				return average_weight_ - getInternalToFullAverageWeight();
			case NTerminal:
				return average_weight_ - getNTerminalToFullAverageWeight();
			case CTerminal:
				return average_weight_ - getCTerminalToFullAverageWeight();
			case BIon:
				return average_weight_ - getBIonToFullAverageWeight();
			case AIon:
				return average_weight_ - getAIonToFullAverageWeight();
			case CIonMinusOne:
				return average_weight_ - getCIonMinusOneToFullAverageWeight();
			case CIon:
				return average_weight_ - EmpiricalFormula("OH").getAverageWeight() + EmpiricalFormula("NH").getAverageWeight();
			case CIonPlusOne:
				return average_weight_ - getCIonPlusOneToFullAverageWeight();
			case CIonPlusTwo:
				return average_weight_ - getCIonPlusTwoToFullAverageWeight();
			case XIon:
				return average_weight_ + getXIonToFullAverageWeight();
			case YIon:
				return average_weight_ + getYIonToFullAverageWeight();
			case ZIonMinusOne:
				return average_weight_ - getZIonMinusOneToFullAverageWeight();
			case ZIon:
				return average_weight_ - getZIonToFullAverageWeight();
			case ZIonPlusOne:
				return average_weight_ - getZIonPlusOneToFullAverageWeight();
			case ZIonPlusTwo:
				return average_weight_ - getZIonPlusTwoToFullAverageWeight();
			default:
				cerr << "Residue::getAverageWeight: unknown ResidueType" << endl;
				return average_weight_;
		}
	}

	void Residue::setMonoWeight(DoubleReal weight)
	{
		mono_weight_ = weight;
		return;
	}

	DoubleReal Residue::getMonoWeight(ResidueType res_type) const
	{
		switch (res_type)
		{
			case Full:
				return mono_weight_;
			case Internal:
				return mono_weight_ - getInternalToFullMonoWeight();
			case NTerminal:
				return mono_weight_ - getNTerminalToFullMonoWeight();
			case CTerminal:
				return mono_weight_ - getCTerminalToFullMonoWeight();
			case BIon:
				return mono_weight_ - getBIonToFullMonoWeight();
			case AIon:
				return mono_weight_ - getAIonToFullMonoWeight();
			case CIonMinusOne:
				return mono_weight_ - getCIonMinusOneToFullMonoWeight();
			case CIon:
				return mono_weight_ - EmpiricalFormula("OH").getMonoWeight() + EmpiricalFormula("NH").getMonoWeight();
			case CIonPlusOne:
				return mono_weight_ - getCIonPlusOneToFullMonoWeight();
			case CIonPlusTwo:
				return mono_weight_ - getCIonPlusTwoToFullMonoWeight();
			case XIon:
				return mono_weight_ + getXIonToFullMonoWeight();
			case YIon:
				return mono_weight_ + getYIonToFullMonoWeight();
     	case ZIonMinusOne:
        return mono_weight_ - getZIonMinusOneToFullMonoWeight();
			case ZIon:
				return mono_weight_ - getZIonToFullMonoWeight();
      case ZIonPlusOne:
        return mono_weight_ - getZIonPlusOneToFullMonoWeight();
      case ZIonPlusTwo:
        return mono_weight_ - getZIonPlusTwoToFullMonoWeight();
			default:
				cerr << "Residue::getMonoWeight: unknown ResidueType" << endl;
				return mono_weight_;
		}
	}

	void Residue::setModification(const String& modification)
	{
		//modification_ = modification;

		ModificationsDB* mod_db = ModificationsDB::getInstance();
		ResidueModification mod = mod_db->getModification(one_letter_code_, modification, ResidueModification::ANYWHERE);

		modification_ = mod.getId();
		// update all the members
		if (mod.getAverageMass() != 0)
		{
			average_weight_ = mod.getAverageMass();
		}
		if (mod.getMonoMass() != 0)
		{
			mono_weight_ = mod.getMonoMass();
		}

		bool updated_formula(false);
		if (mod.getDiffFormula() != "")
		{
			updated_formula = true;
			setFormula(getFormula() + mod.getDiffFormula());
		}
		if (mod.getFormula() != "" && !updated_formula)
		{
			updated_formula = true;
			String formula = mod.getFormula();
			formula.removeWhitespaces();
			formula_ = formula;
		}
		
		if (updated_formula)
		{
			average_weight_ = formula_.getAverageWeight();
			mono_weight_ = formula_.getMonoWeight();
		}
		else
		{
			if (mod.getAverageMass() != 0)
			{
				average_weight_ = mod.getAverageMass();
			}
			if (mod.getMonoMass() != 0)
			{
				mono_weight_ = mod.getMonoMass();
			}
		}
	
		// neutral losses
		loss_formulas_.clear();
		loss_names_.clear();
		if (mod.hasNeutralLoss())
		{
			loss_formulas_.push_back(mod.getNeutralLossDiffFormula());
			loss_names_.push_back(mod.getNeutralLossDiffFormula().getString());
		}

		is_modified_ = true;
	}
	
	const String& Residue::getModification() const
	{
		return modification_;
	}

	void Residue::setLowMassIons(const vector<EmpiricalFormula>& low_mass_ions)
	{
		low_mass_ions_ = low_mass_ions;
	}

	const vector<EmpiricalFormula>& Residue::getLowMassIons() const
	{
		return low_mass_ions_;
	}

	DoubleReal Residue::getBackboneBasicityRight() const
	{
		return gb_bb_r_;
	}

	void Residue::setBackboneBasicityRight(DoubleReal gb_bb_r)
	{
		gb_bb_r_ = gb_bb_r;
	}

	DoubleReal Residue::getBackboneBasicityLeft() const
	{
		return gb_bb_l_;
	}

	void Residue::setBackboneBasicityLeft(DoubleReal gb_bb_l)
	{
		gb_bb_l_ = gb_bb_l;
	}

	DoubleReal Residue::getSideChainBasicity() const
	{
		return gb_sc_;
	}

	void Residue::setSideChainBasicity(DoubleReal gb_sc)
	{
		gb_sc_ = gb_sc;
	}

	void Residue::setResidueSets(const set<String>& residue_sets)
	{
		residue_sets_ = residue_sets;
	}

	const set<String>& Residue::getResidueSets() const
	{
		return residue_sets_;
	}

	void Residue::addResidueSet(const String& residue_set)
	{
		residue_sets_.insert(residue_set);
	}

	bool Residue::isModified() const
	{
		return is_modified_;
	}
	
	bool Residue::hasNeutralLoss() const
	{
    return ( !loss_formulas_.empty() );
	}

	bool Residue::hasNTermNeutralLosses() const
  {
    return ( !NTerm_loss_formulas_.empty() );
	}
	
	bool Residue::operator == (const Residue& residue) const
	{
		return (name_ == residue.name_ &&
						short_name_ == residue.short_name_ &&
						synonyms_ == residue.synonyms_ &&
						three_letter_code_ == residue.three_letter_code_ &&
						one_letter_code_ == residue.one_letter_code_ &&
						formula_ == residue.formula_ &&
						average_weight_ == residue.average_weight_ &&
						mono_weight_ == residue.mono_weight_ &&
						is_modified_ == residue.is_modified_ &&
						pre_mod_name_ == residue.pre_mod_name_ &&
						modification_ == residue.modification_ &&
						loss_names_ == residue.loss_names_ &&
						loss_formulas_ == residue.loss_formulas_ &&
						NTerm_loss_names_ == residue.NTerm_loss_names_ &&
						NTerm_loss_formulas_ == residue.NTerm_loss_formulas_ &&
						loss_average_weight_ == residue.loss_average_weight_ &&
						loss_mono_weight_ == residue.loss_mono_weight_ &&
						low_mass_ions_ == residue.low_mass_ions_ &&
						pka_ == residue.pka_ &&
						pkb_ == residue.pkb_ &&
						pkc_ == residue.pkc_ &&
						gb_sc_ == residue.gb_sc_ &&
						gb_bb_l_ == residue.gb_bb_l_ &&
						gb_bb_r_ == residue.gb_bb_r_ &&
						residue_sets_ == residue.residue_sets_);
	}

	bool Residue::operator == (char one_letter_code) const
	{
		return one_letter_code_[0] == one_letter_code;
	}

	bool Residue::operator != (char one_letter_code) const
	{
		return one_letter_code_[0] != one_letter_code;
	}

	bool Residue::operator != (const Residue& residue) const
	{
		return !(*this == residue);
	}

	bool Residue::isInResidueSet(const String& residue_set)
	{
		return residue_sets_.find(residue_set) != residue_sets_.end();
	}

	ostream& operator << (ostream& os, const Residue& residue) 
	{
		os << residue.name_ << " "
			 << residue.three_letter_code_ << " "
			 << residue.one_letter_code_ << " "
			 << residue.formula_;
		return os;
	}

}

