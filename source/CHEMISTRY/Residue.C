// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/Residue.h>

using namespace std;

namespace OpenMS
{
	// modification
	Residue::Modification::Modification()
		: add_average_weight_(0.0f),
			add_mono_weight_(0.0f),
			del_average_weight_(0.0f),
			del_mono_weight_(0.0f)
	{
	}

	Residue::Modification::Modification(const Modification& modification)
		:	name_(modification.name_),
			short_name_(modification.short_name_),
			name_prefix_(modification.name_prefix_),
			synonyms_(modification.synonyms_),
			add_formula_(modification.add_formula_),
			add_average_weight_(modification.add_average_weight_),
			add_mono_weight_(modification.add_mono_weight_),
			del_formula_(modification.del_formula_),
			del_average_weight_(modification.del_average_weight_),
			del_mono_weight_(modification.del_mono_weight_),
			valid_residues_(modification.valid_residues_)
	{
	}

	Residue::Modification::~Modification()
	{
	}

	Residue::Modification& Residue::Modification::operator = (const Modification& modification)
	{
		if (this != &modification)
		{
			name_ = modification.name_;
			name_prefix_ = modification.name_prefix_;
			short_name_ = modification.short_name_;
			synonyms_ = modification.synonyms_;
			add_formula_ = modification.add_formula_;
			add_average_weight_ = modification.add_average_weight_;
			add_mono_weight_ = modification.add_mono_weight_;
			del_formula_ = modification.del_formula_;
			del_average_weight_ = modification.del_average_weight_;
			del_mono_weight_ = modification.del_mono_weight_;
			valid_residues_ = modification.valid_residues_;
		}
		return *this;
	}
	
	void Residue::Modification::setName(const String& name)
	{
		name_ = name;
	}

	const String& Residue::Modification::getName() const
	{
		return name_;
	}

	void Residue::Modification::setShortName(const String& short_name)
	{
		short_name_ = short_name;
	}

	const String& Residue::Modification::getShortName() const
	{
		return short_name_;
	}

	void Residue::Modification::setNamePrefix(const String& prefix)
	{
		name_prefix_ = prefix;
	}

	const String& Residue::Modification::getNamePrefix() const
	{
		return name_prefix_;
	}

	void Residue::Modification::setSynonyms(const std::set<String>& synonyms)
	{
		synonyms_ = synonyms;
	}

	void Residue::Modification::addSynonym(const String& synonym)
	{
		synonyms_.insert(synonym);
	}

	const std::set<String>& Residue::Modification::getSynonyms() const
	{
		return synonyms_;
	}
	
	void Residue::Modification::setAddFormula(const EmpiricalFormula& formula)
	{
		add_formula_ = formula;
	}

	const EmpiricalFormula& Residue::Modification::getAddFormula() const
	{
		return add_formula_;
	}

	void Residue::Modification::setAddAverageWeight(Real weight)
	{
		add_average_weight_ = weight;
	}

	Real Residue::Modification::getAddAverageWeight() const
	{
		return add_average_weight_;
	}

	void Residue::Modification::setAddMonoWeight(Real weight)
	{
		add_mono_weight_ = weight;
	}

	Real Residue::Modification::getAddMonoWeight() const
	{
		return add_mono_weight_;
	}

	void Residue::Modification::setDelFormula(const EmpiricalFormula& formula)
	{
		del_formula_ = formula;
	}

	const EmpiricalFormula& Residue::Modification::getDelFormula() const
	{
		return del_formula_;
	}

	void Residue::Modification::setDelAverageWeight(Real weight)
	{
		del_average_weight_ = weight;
	}

	Real Residue::Modification::getDelAverageWeight() const
	{
		return del_average_weight_;
	}
	
	void Residue::Modification::setDelMonoWeight(Real weight)
	{
		del_mono_weight_ = weight;
	}

	Real Residue::Modification::getDelMonoWeight() const
	{
		return del_mono_weight_;
	}

	void Residue::Modification::setValidResidues(const set<Residue*>& valid_residues)
	{
		valid_residues_ = valid_residues;
	}

	void Residue::Modification::addValidResidue(Residue* valid_residue)
	{
		valid_residues_.insert(valid_residue);
	}
	
	const set<Residue*>& Residue::Modification::getValidResidues() const
	{
		return valid_residues_;
	}

	bool Residue::Modification::operator == (const Residue::Modification& modification) const
	{
		return 	name_ == modification.name_ &&
						name_prefix_ == modification.name_prefix_ &&
						synonyms_ == modification.synonyms_ &&
						add_formula_ == modification.add_formula_ &&
						add_average_weight_ == modification.add_average_weight_ &&
						add_mono_weight_ == modification.add_mono_weight_ &&
						del_formula_ == modification.del_formula_ &&
						del_average_weight_ == modification.del_average_weight_ &&
						del_mono_weight_ == modification.del_mono_weight_ &&
						valid_residues_ == modification.valid_residues_;
	}
	
	bool Residue::Modification::operator != (const Residue::Modification& modification) const
	{
		return !(*this == modification);
	}

	// residue
	Residue::Residue()
		: name_("unknown"),
			average_weight_(0.0f),
			mono_weight_(0.0f),
			is_modified_(false),
			modification_(0),
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
										const EmpiricalFormula& formula,
										const EmpiricalFormula& neutral_loss)
		:	name_(name),
			three_letter_code_(three_letter_code),
			one_letter_code_(one_letter_code),
			formula_(formula),
			internal_formula_(formula_ - getInternalToFull()),
			average_weight_(0),
			mono_weight_(0),
			is_modified_(false),
			modification_(0),
			loss_formula_(neutral_loss),
			loss_average_weight_(0.0f),
			loss_mono_weight_(0.0f),
			pka_(0.0),
			pkb_(0.0),
			pkc_(-1.0)
	{
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
			loss_name_(residue.loss_name_),
			loss_formula_(residue.loss_formula_),
			loss_average_weight_(residue.loss_average_weight_),
			loss_mono_weight_(residue.loss_mono_weight_),
			low_mass_ions_(residue.low_mass_ions_),
			pka_(residue.pka_),
			pkb_(residue.pkb_),
			pkc_(residue.pkc_)
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
			loss_name_ = residue.loss_name_;
			loss_formula_ = residue.loss_formula_;
			loss_average_weight_ = residue.loss_average_weight_;
			loss_mono_weight_ = residue.loss_mono_weight_;
			low_mass_ions_ = residue.low_mass_ions_;
			pka_ = residue.pka_;
			pkb_ = residue.pkb_;
			pkc_ = residue.pkc_;
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

	String getResidueTypeName(Residue::ResidueType res_type)
	{
		String ion("-ion");
		switch (res_type)
		{
			case Residue::AIon:
				return "a"+ion;
			case Residue::BIon:
				return "b"+ion;
			case Residue::CIon:
				return "c"+ion;
			case Residue::XIon:
				return "x"+ion;
			case Residue::YIon:
				return "y"+ion;
			case Residue::ZIon:
				return "z"+ion;
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

	void Residue::setLossFormula(const EmpiricalFormula& loss_formula)
	{
		loss_formula_ = loss_formula;
	}

	const EmpiricalFormula& Residue::getLossFormula() const
	{
		return loss_formula_;
	}

	void Residue::setLossAverageWeight(Real weight)
	{
		loss_average_weight_ = weight;
	}

	Real Residue::getLossAverageWeight() const
	{
		return loss_average_weight_;
	}

	void Residue::setLossMonoWeight(Real weight)
	{
		loss_mono_weight_ = weight;
	}

	Real Residue::getLossMonoWeight() const
	{
		return loss_mono_weight_;
	}

	void Residue::setLossName(const String& name)
	{
		loss_name_ = name;
	}

	const String& Residue::getLossName() const
	{
		return loss_name_;
	}
	
	void Residue::setFormula(const EmpiricalFormula& formula, ResidueType res_type)
	{
		switch (res_type)
		{
			case Full: 
				formula_ = formula; 
				internal_formula_ = formula_ - getInternalToFull();
				return;
			case Internal:
				formula_ = formula + getInternalToFull();
				internal_formula_ = formula_ - getInternalToFull();
				return;
			case NTerminal:
				formula_ = formula + getNTerminalToFull();
				internal_formula_ = formula_ - getInternalToFull();
				return;
			case CTerminal:
				formula_ = formula + getCTerminalToFull();
				internal_formula_ = formula_ - getInternalToFull();
				return;
			case BIon:
				formula_ = formula + getBIonToFull();
				internal_formula_ = formula_ - getInternalToFull();
				return;
			case YIon:
				formula_ = formula - getYIonToFull();
				internal_formula_ = formula_ - getInternalToFull();
				return;
			case AIon:
				formula_ = formula + getAIonToFull();
				internal_formula_ = formula_ - getInternalToFull();
				return;
			default: 
				cerr << "Residue::setFormula: unknown ResidueType" << endl;
				return;
		}
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
			case CIon:
				return formula_ - EmpiricalFormula("OH") + EmpiricalFormula("NH");
			case XIon:
				return formula_ + getXIonToFull();
			case YIon:
				return formula_ + getYIonToFull();
			case ZIon:
				return formula_ - getZIonToFull();
			default:
				cerr << "Residue::getFormula: unknown ResidueType" << endl;
				return formula_;
		}
	}

	void Residue::setAverageWeight(Real weight, ResidueType res_type) 
	{
		switch (res_type)
		{
			case Full:
				average_weight_ = weight;
				return;
			case Internal:
				average_weight_ = weight + getInternalToFullAverageWeight();
				return;
			case NTerminal:
				average_weight_ = weight + getNTerminalToFullAverageWeight();
				return;
			case CTerminal:
				average_weight_ = weight + getCTerminalToFullAverageWeight();
				return;
			case BIon:
				average_weight_ = weight + getBIonToFullAverageWeight();
				return;
			case AIon:
				average_weight_ = weight + getAIonToFullAverageWeight();
				return;
			case YIon:
				average_weight_ = weight - getYIonToFullAverageWeight();
				return;
			default:
				cerr << "Residue::setAverageWeight: unknown ResidueType" << endl;
				average_weight_ = weight;
				return;
		}
	}

	Real Residue::getAverageWeight(ResidueType res_type) const
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
			case CIon:
				return average_weight_ - EmpiricalFormula("OH").getAverageWeight() + EmpiricalFormula("NH").getAverageWeight();
			case XIon:
				return average_weight_ + getXIonToFullAverageWeight();
			case YIon:
				return average_weight_ + getYIonToFullAverageWeight();
			case ZIon:
				return average_weight_ - getZIonToFullAverageWeight();
			default:
				cerr << "Residue::getAverageWeight: unknown ResidueType" << endl;
				return average_weight_;
		}
	}

	void Residue::setMonoWeight(Real weight, ResidueType res_type)
	{
		switch (res_type)
		{
			case Full:
				mono_weight_ = weight;
				return;
			case Internal:
				mono_weight_ = weight + getInternalToFullMonoWeight();
				return;
			case NTerminal:
				mono_weight_ = weight + getNTerminalToFullMonoWeight();
				return;
			case CTerminal:
				mono_weight_ = weight + getCTerminalToFullMonoWeight();
				return;
			case BIon:
				mono_weight_ = weight + getBIonToFullMonoWeight();
				return;
			case AIon:
				mono_weight_ = weight + getAIonToFullMonoWeight();
				return;
			case YIon:
				mono_weight_ = weight - getYIonToFullMonoWeight();
				return;
			default:
				cerr << "Residue::setMonoWeight: unknown ResidueType" << endl;
				mono_weight_ = weight;
				return;
		}
	}

	Real Residue::getMonoWeight(ResidueType res_type) const
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
			case CIon:
				return mono_weight_ - EmpiricalFormula("OH").getAverageWeight() + EmpiricalFormula("NH").getAverageWeight();
			case XIon:
				return mono_weight_ + getXIonToFullMonoWeight();
			case YIon:
				return mono_weight_ + getYIonToFullMonoWeight();
			case ZIon:
				return mono_weight_ - getZIonToFullMonoWeight();
			default:
				cerr << "Residue::getMonoWeight: unknown ResidueType" << endl;
				return mono_weight_;
		}
	}

	void Residue::setModification(Modification* modification)
	{
		modification_ = modification;
		is_modified_ = true;
	}
	
	const Residue::Modification* Residue::getModification() const
	{
		return modification_;
	}

	void Residue::setUnmodifiedName(const String& name)
	{
		pre_mod_name_ = name;
	}
	
	const String& Residue::getUnmodifiedName() const
	{
		return pre_mod_name_;
	}

	void Residue::setLowMassIons(const vector<EmpiricalFormula>& low_mass_ions)
	{
		low_mass_ions_ = low_mass_ions;
	}

	const vector<EmpiricalFormula>& Residue::getLowMassIons() const
	{
		return low_mass_ions_;
	}
	
	bool Residue::isModified() const
	{
		return is_modified_;
	}
	
	bool Residue::hasNeutralLoss() const
	{
		return !loss_formula_.isEmpty();
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
						loss_name_ == residue.loss_name_ &&
						loss_formula_ == residue.loss_formula_ &&
						loss_average_weight_ == residue.loss_average_weight_ &&
						loss_mono_weight_ == residue.loss_mono_weight_ &&
						low_mass_ions_ == residue.low_mass_ions_);
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

	ostream& operator << (ostream& os, const Residue& residue) 
	{
		os << residue.name_ << " "
			 << residue.three_letter_code_ << " "
			 << residue.one_letter_code_ << " "
			 << residue.formula_;
		return os;
	}

}

