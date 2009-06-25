// -*- mode: C++; tab-width: 2; -*-
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
//

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <algorithm>

using namespace std;

namespace OpenMS 
{
	//ResidueDB * AASequence::custom_res_db_ = 0;

	// AASequence
	AASequence::AASequence()
		:	valid_(true)
	{
	}

	AASequence::AASequence(const AASequence& rhs)
		:	peptide_(rhs.peptide_),
			sequence_string_(rhs.sequence_string_),
			valid_(rhs.valid_),
			n_term_mod_(rhs.n_term_mod_),
			c_term_mod_(rhs.c_term_mod_)
	{
	}

	AASequence::AASequence(const String& peptide)
		:	valid_(true)
	{
		parseString_(peptide_, peptide);
	}

	AASequence::AASequence(const char* peptide)
		: valid_(true)
	{
		parseString_(peptide_, String(peptide));
	}
	
	AASequence::~AASequence()
	{
	}

	AASequence& AASequence::operator = (const AASequence& rhs)
	{
		if (this != &rhs)
		{
			peptide_ = rhs.peptide_;
			sequence_string_ = rhs.sequence_string_;
			valid_ = rhs.valid_;
			n_term_mod_ = rhs.n_term_mod_;
			c_term_mod_ = rhs.c_term_mod_;
		}
		return *this;
	}
	
	const Residue& AASequence::getResidue(SignedSize index) const
	{
		if (index >= 0 && Size(index) >= peptide_.size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peptide_.size());
		}
		if (index < 0)
		{
			throw Exception::IndexUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, 0);
		}
		return *peptide_[index];
	}

	const Residue& AASequence::getResidue(Size index) const
	{
		if (index >= peptide_.size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peptide_.size());
		}
		return *peptide_[index];
	}

	bool AASequence::isValid() const
	{
		return valid_;
	}
	
	String AASequence::toString() const
	{
		stringstream ss;
		ss << *this;
		return String(ss.str());
	}
	
	String AASequence::toUnmodifiedString() const
	{
		if (valid_)
		{
			String tmp;
			for (ConstIterator it = begin(); it != end(); ++it)
			{
				tmp += it->getOneLetterCode();
			}
			return tmp;
		}
		return sequence_string_;
	}


	bool AASequence::operator < (const AASequence& rhs) const
	{
		if (!valid_)
		{
		if (!rhs.valid_)
			{
				return sequence_string_ < rhs.sequence_string_;
			}
			else
			{
				return sequence_string_ < rhs.toString();
			}
		}
		else
		{
			if (!rhs.valid_)
			{
				return toString() < rhs.sequence_string_;
			}
			else
			{
				return toString() < rhs.toString();
			}
		}
		return false;
	}
	
	EmpiricalFormula AASequence::getFormula(Residue::ResidueType type, Int charge) const
	{
		EmpiricalFormula ef;
		static EmpiricalFormula H("H");
		static EmpiricalFormula OH("OH");
		static EmpiricalFormula NH("NH");
		for (Int i=0; i<charge; ++i)
		{
			ef += H;
		}


    // terminal modifications
    if (n_term_mod_ != "" && 
				(type == Residue::Full || type == Residue::AIon || type == Residue::BIon || type == Residue::CIon || type == Residue::NTerminal)
			 )
    {
      ef += ModificationsDB::getInstance()->getModification(n_term_mod_).getDiffFormula();
    }


    if (c_term_mod_ != "" &&
				(type == Residue::Full || type == Residue::XIon || type == Residue::YIon || type == Residue::ZIon || type == Residue::CTerminal)
			 )
    {
      ef += ModificationsDB::getInstance()->getModification(c_term_mod_).getDiffFormula();
    }



		if (peptide_.size() > 0)
		{
			if (peptide_.size() == 1)
			{	
				ef += peptide_[0]->getFormula(type);
			}
			else
			{
				for (Size i=0;i!=peptide_.size();++i)
				{
					ef += peptide_[i]->getFormula(Residue::Internal);
				}
				
				// add the missing formula part
				switch(type)
				{
					case Residue::Full:
						return ef+Residue::getInternalToFull();
					case Residue::Internal: 
						return ef/* + add_protons*/;
					case Residue::NTerminal:
						return ef + Residue::getInternalToFull() - Residue::getNTerminalToFull();
					case Residue::CTerminal:
						return ef + Residue::getInternalToFull() - Residue::getCTerminalToFull();
					case Residue::BIon:
						return ef + Residue::getInternalToFull() - Residue::getBIonToFull() - H;
					case Residue::AIon:
						return ef + Residue::getInternalToFull() - Residue::getAIonToFull() - H;
					case Residue::CIon:
						return ef + Residue::getInternalToFull() - OH + NH;
					case Residue::XIon:
						return ef + Residue::getInternalToFull() + Residue::getXIonToFull();
					case Residue::YIon:
						return ef + Residue::getInternalToFull() + Residue::getYIonToFull();
					case Residue::ZIon:
						return ef + Residue::getInternalToFull() - Residue::getZIonToFull();
					default:
						cerr << "AASequence::getFormula: unknown ResidueType" << endl;
				}
			}			
		}

		return ef;
	}
	
	DoubleReal AASequence::getAverageWeight(Residue::ResidueType type, Int charge) const
	{
		return getFormula(type, charge).getAverageWeight();
	}

	DoubleReal AASequence::getMonoWeight(Residue::ResidueType type, Int charge) const
	{
		return getFormula(type, charge).getMonoWeight();
	}

	/*void AASequence::getNeutralLosses(Map<const EmpiricalFormula, UInt) const
	{
		// the following losses are from the Zhang paper (AC, 76, 14, 2004)
		// charge directed*/
		/*
		static const EmpiricalFormula R_44("NH2CHNH"); 
		static const EmpiricalFormula R_59("CN3H5"); // guanidium
		static const EmpiricalFormula R_61("N2H4CH");
		// charge remote
		static const EmpiricalFormula R_60("N2H4CO"); // combination of NH=C=NH + C-terminal H2O
		static const EmpiricalFormula H2O("H2O"); // loss from the C-terminus
		static const EmpiricalFormula NH3("NH3");
		Map<const EmpiricalFormula*, UInt> losses;
		
		for (Size i=0;i!=peptide_.size();++i)
		{
			if (peptide_[i]->hasNeutralLoss())
			{
				const EmpiricalFormula* loss = peptide_[i]->getLossFormulas();
				if (losses.find(loss) != losses.end())
				{
					losses[loss]++;
				}
				else
				{
					losses[loss] = 1;
				}
			}
	
			
			// TODO: hack this should be in the data file
			if (peptide_[i]->getOneLetterCode() == "R")
			{
				losses[&R_44] = 1;
				losses[&R_59] = 1;
				losses[&R_61] = 1;
				losses[&R_60] = 1;
			}
			losses[&H2O] = 1;
			losses[&NH3] = 1;
		}	
		return losses;
	}*/

	const Residue& AASequence::operator [] (SignedSize index) const
	{
		if (index < 0)
		{
			throw Exception::IndexUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, 0);
		}
		else
		{
			if (Size(index) >= size())
			{
				throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
			}
		}
		return *peptide_[Size(index)];
	}
	
	const Residue& AASequence::operator [] (Size index) const
	{
		if (index >= size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
		}
		return *peptide_[index];
	}
	
	AASequence AASequence::operator + (const AASequence& sequence) const
	{
		AASequence seq;
		seq.peptide_ = peptide_;
		for (Size i=0;i!=sequence.peptide_.size();++i)
		{
			seq.peptide_.push_back(sequence.peptide_[i]);
		}
		return seq;
	}
	
	AASequence AASequence::operator + (const String& peptide) const
	{
		AASequence seq(peptide);
		return *this + seq;
	}

	AASequence AASequence::operator + (const char* peptide) const
	{
		return *this + String(peptide);
	}

	AASequence AASequence::operator + (const Residue* residue) const
	{
		if (!ResidueDB::getInstance()->hasResidue(residue))
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, "given residue");
		}
		AASequence seq = *this;
		seq += residue;
		return seq;
	}
	
	AASequence& AASequence::operator += (const AASequence& sequence)
	{
		for (Size i=0;i!=sequence.peptide_.size();++i)
		{
			peptide_.push_back(sequence.peptide_[i]);
		}
		return *this;
	}
	
	AASequence& AASequence::operator += (const String& peptide)
	{
		vector<const Residue*> vec;
		parseString_(vec, peptide);
		for (Size i=0;i!=vec.size();++i)
		{
			peptide_.push_back(vec[i]);
		}
		return *this;
	}

	AASequence& AASequence::operator += (const char* peptide)
	{
		*this += String(peptide);
		return *this;
	}


	AASequence& AASequence::operator += (const Residue* residue)
	{
		if (!ResidueDB::getInstance()->hasResidue(residue))
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, "given residue");
		}
		peptide_.push_back(residue);
		return *this;
	}
	
	Size AASequence::size() const
	{
		return peptide_.size();
	}

	AASequence AASequence::getPrefix(Size index) const
	{
		if (index > size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
		}
		AASequence seq;
		for (Size i=0;i<index;++i)
		{
			seq.peptide_.push_back(peptide_[i]);
		}
		return seq;
	}

	AASequence AASequence::getSuffix(Size index) const
	{
		if (index > size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
		}
		AASequence seq;
		for (Size i=size()-index;i!=size();++i)
		{
			seq.peptide_.push_back(peptide_[i]);
		}
		return seq;
	}

	AASequence AASequence::getSubsequence(Size index, UInt num) const
	{
		if (index >= size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
		}
		if (index + num > size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index+num, size());
		}
		AASequence seq;
		for (Size i = index; i != index + num; ++i)
		{
			seq.peptide_.push_back(peptide_[i]);
		}
		return seq;
	}
	
	bool AASequence::has(const Residue& residue) const 
	{
		for (Size i = 0; i != peptide_.size(); ++i)
		{
			if (*peptide_[i] == residue)
			{
				return true;
			}
		}
		return false;
	}

	bool AASequence::has(const String& residue) const
	{
		if (!getResidueDB_()->hasResidue(residue))
		{
			return false;
		}
		return has(*getResidueDB_()->getResidue(residue));
	}

	
	bool AASequence::hasSubsequence(const AASequence& sequence) const
	{
		if (sequence.size() == 0)
		{
			return true;
		}
		else
		{
			if (sequence.size() <= peptide_.size())
			{
				for (Size i=0;i!=peptide_.size();++i)
				{
					if (peptide_[i] == sequence.peptide_[0])
					{
						Size j=0;
						for (;j+i!=peptide_.size() && j!=sequence.peptide_.size();++j)
						{
							if (peptide_[j+i] == sequence.peptide_[j])
							{
								if (j == sequence.peptide_.size() -1)
								{
									return true;
								}
							}
							else
							{
								break;
							}
						}
					}
				}
			}
		}
		return false;
	}

	bool AASequence::hasSubsequence(const String& sequence) const
	{
		AASequence aa_seq(sequence);
		vector<const Residue*> seq = aa_seq.peptide_;
		if (seq.size() == 0)
		{
			return true;
		}
		else
		{
			if (seq.size() <= peptide_.size())
			{
				for (Size i=0;i!=peptide_.size();++i)
				{
					Size j=0;
					for (;j+i!=peptide_.size() && j!=seq.size();++j)
					{
						if (peptide_[j+i] == seq[j])
						{
							if (j == seq.size()-1)
							{
								return true;
							}
						}
						else
						{
							break;
						}
					}	
				}
			}
		}
		return false;
	}
	
	bool AASequence::hasPrefix(const AASequence& sequence) const
	{
		if (sequence.size() == 0)
		{
			return true;
		}
		if (sequence.size() > peptide_.size())
		{
			return false;
		}
		for (Size i=0;i!=sequence.size();++i)
		{
			if (sequence.peptide_[i] != peptide_[i])
			{
				return false;
			}
		}
		return true;
	}

	bool AASequence::hasPrefix(const String& sequence) const
	{
		AASequence seq(sequence);
		return hasPrefix(seq);
	}
	
	bool AASequence::hasSuffix(const AASequence& sequence) const
	{
		if (sequence.size() == 0)
		{
			return true;
		}
		if (sequence.size() > peptide_.size())
		{
			return false;
		}
		for (Size i=0;i!=sequence.size();++i)
		{
			if (sequence.peptide_[sequence.size()-1-i] != peptide_[size()-1-i])
			{
				return false;
			}
		}
		return true;
	}

	bool AASequence::hasSuffix(const String& sequence) const
	{
		AASequence seq(sequence);
		return hasSuffix(seq);
	}
	
	bool AASequence::operator == (const AASequence& peptide) const
	{
		if (!valid_)
		{
			if (peptide.valid_)
			{
				return false;
			}
			else
			{
				return sequence_string_ == peptide.sequence_string_ && 
							 n_term_mod_ == peptide.n_term_mod_ &&
							 c_term_mod_ == peptide.c_term_mod_;
			}
		}
		
		if (size() != peptide.size())
		{
			return false;
		}
		for (Size i = 0; i != size(); ++i)
		{
			if (peptide_[i] != peptide.peptide_[i])
			{
				return false;
			}
		}
		
		if (n_term_mod_ != peptide.n_term_mod_)
		{
			return false;
		}
		if (c_term_mod_ != peptide.c_term_mod_)
		{
			return false;
		}
		return true;
	}
	
	bool AASequence::operator == (const String& peptide) const
	{
		AASequence sequence(peptide);
		return *this == sequence;
	}

	bool AASequence::operator == (const char* peptide) const
	{
		return *this == String(peptide);
	}
	
	bool AASequence::operator != (const AASequence& peptide) const
	{
		return !(*this == peptide);
	}

	bool AASequence::operator != (const String& sequence) const
	{
		return !(*this == sequence);
	}

	bool AASequence::operator != (const char* sequence) const
	{
		return *this != String(sequence);
	}

	bool AASequence::isModified() const
	{
		if (n_term_mod_ != "" || c_term_mod_ != "")
		{
			return true;
		}
		for (vector<const Residue*>::const_iterator it = peptide_.begin(); it != peptide_.end(); ++it)
		{
			if ((*it)->isModified())
			{
				return true;
			}
		}
		return false;
	}
	
	bool AASequence::isModified(Size position) const
	{
		if (position >= peptide_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide_.size(), position);
    }

		return peptide_[position]->isModified();
	}
	
	ostream& operator << (ostream& os, const AASequence& peptide)
	{
		if (!peptide.valid_)
		{
			os << peptide.sequence_string_;
			return os;
		}
    if (peptide.n_term_mod_ != "")
    {
      os << "(" << peptide.n_term_mod_ << ")";
    }

		for (Size i = 0; i != peptide.size(); ++i)
		{
			if (peptide.peptide_[i]->isModified())
			{
				if (peptide.peptide_[i]->getOneLetterCode() != "")
				{
					os << peptide.peptide_[i]->getOneLetterCode();
				}
				else
				{
					os << "[" << peptide.peptide_[i]->getMonoWeight() << "]";
				}
				
				String id = ModificationsDB::getInstance()->getModification(peptide.peptide_[i]->getOneLetterCode(), peptide.peptide_[i]->getModification()).getId();
				if (id != "")
				{
					os << "(" << id << ")";
				}
				else
				{
					os << "([" << ModificationsDB::getInstance()->getModification(peptide.peptide_[i]->getOneLetterCode(), peptide.peptide_[i]->getModification()).getDiffMonoMass() << "])";
				}
			}
			else
			{
				if (peptide.peptide_[i]->getOneLetterCode() != "")
				{
					os << peptide.peptide_[i]->getOneLetterCode();
				}
				else
				{
					if (peptide.peptide_[i]->getShortName() != "")
					{
						os << peptide.peptide_[i]->getShortName();
					}
					else
					{
						os << "[" << peptide.peptide_[i]->getMonoWeight() << "]";
					}
				}
			}
		}
		
    if (peptide.c_term_mod_ != "")
    {
      os << "(" << peptide.c_term_mod_ << ")";
    }
		return os;
	}

	
	void AASequence::parseString_(vector<const Residue*>& sequence, const String& pep)
	{
		sequence.clear();
		String peptide(pep);
		peptide.trim();

		if (peptide.size() == 0)
		{
			return;
		}
		
		// split the peptide in its residues
		vector<String> split;
		Size pos(0);
		bool mod_open(false), tag_open(false);
		if (peptide[0] == '[')
		{
			tag_open = true;
		}
		if (peptide[0] == '(')
		{
			mod_open = true;
		}
		for (Size i = 1; i < peptide.size(); ++i)
		{
			if ((isalpha(peptide[i]) && isupper(peptide[i]) && !mod_open) ||
					(peptide[i] == '[' && !mod_open))
			{
				split.push_back(peptide.substr(pos, i-pos));
				pos = i;
				if (mod_open)
				{
					mod_open = false;
				}
			}
			if (peptide[i] == '(')
			{
				mod_open = true;
				continue;
			}
			if (peptide[i] == ')')
			{
				mod_open = false;
				continue;
			}
			if (peptide[i] == '[')
			{
				tag_open = true;
				continue;
			}
			if (peptide[i] == ']')
			{
				tag_open = false;
				continue;
			}
		}
		
		// push_back last residue
		split.push_back(peptide.substr(pos, peptide.size()-pos));

		if (split.size() > 0 && split[0].size() > 0 && split[0][0] == '(')
		{
			String mod = split[0];
			mod.remove('(');
			mod.remove(')');
			n_term_mod_ = ModificationsDB::getInstance()->getModification(mod).getId();

			split.erase(split.begin());
		}

		// test the last split if there is a C-terminal modification
		if (split.size() == 0)
		{
			return;
		}

		String c_term = *(split.end() - 1);
		Size c_term_mods = count(c_term.begin(), c_term.end(), '(');
		if (c_term_mods > 0)
		{
			// now we have found a potential C-term modification
			String mod;
			mod = c_term.suffix('(');
			mod = mod.prefix(')');
			const ResidueModification* potential_mod = &ModificationsDB::getInstance()->getModification(mod);
			if (potential_mod->getTermSpecificity() == ResidueModification::C_TERM)
			{
				c_term_mod_ = potential_mod->getId();			
				split[split.size() - 1] = c_term.substr(0, c_term.size() - mod.size() - 2);
			}
		}
		
		// parse the residues
		for (Size i = 0; i != split.size(); ++i)
		{
			String res = split[i];
			String name, mod, tag;
			for (Size j = 0; j != res.size(); ++j)
			{
				if (isalpha(res[j]))
				{
					name += res[j];
				}
				else
				{
					if (res[j] == '(')
					{
						for (Size k = j + 1; res[k] != ')';++k)
						{
							if (k == res.size())
							{
								valid_ = false;
								sequence_string_.concatenate(split.begin(), split.end());
								sequence.clear();
								cerr << "AASequence: cannot convert string '" << peptide << "' into meaningful amino acid sequence, missing ')'!" << endl;
								return;
							}
							mod += res[k];
						}
						break;
					}
					else
					{
						if (res[j] == '[')
						{
							for (Size k = j + 1; res[k] != ']'; ++k)
							{
								if (k == res.size())
								{
									valid_ = false;
									sequence_string_.concatenate(split.begin(), split.end());
									sequence.clear();
									cerr << "AASequence: cannot convert string '" << peptide << "' into meaningful amino acid sequence, missing ']'!" << endl;
									return;
								}
								tag += res[k];
							}
							break;
						}
						else
						{
							valid_ = false;
							sequence_string_.concatenate(split.begin(), split.end());
							sequence.clear();
							cerr << "AASequence: cannot convert string '" << peptide << "' into meaningful amino acid sequence, residue '" << res << "' unknown at position " << j << "!" << endl;
							return;
						}
					}
				}
			}

			// now we have the name and the modification name (if there is one), or the tag
			const Residue* res_ptr = getResidueDB_()->getResidue(name);
			if (res_ptr == 0 && tag == "")
			{
				valid_ = false;
				sequence_string_.concatenate(split.begin(), split.end());
				sequence.clear();
				cerr << "AASequence: cannot parse residue with name: '" << name << "' from sequence '" << peptide << "'" << endl;
				return;
			}
			if (mod != "")
			{
				sequence.push_back(ResidueDB::getInstance()->getModifiedResidue(res_ptr, mod));
			}
			else
			{
				if (tag != "")
				{
					// if the residue db does not have this tag-residue, we add one
					if (res_ptr == 0)
					{
						Residue res(tag, String(""), String(""), EmpiricalFormula(""));
						res.setMonoWeight(tag.toDouble(), Residue::Internal);
						res.setAverageWeight(tag.toDouble(), Residue::Internal);
						getResidueDB_()->addResidue(res);
						sequence.push_back(getResidueDB_()->getResidue(tag));
					}
					else
					{
						sequence.push_back(res_ptr);
					}
				}
				else
				{
					sequence.push_back(res_ptr);
				}
				if (sequence.size() < i)
				{
					valid_ = false;
					sequence_string_.concatenate(split.begin(), split.end());
					sequence.clear();
					return;
				}
			}
		}
	}

	ResidueDB* AASequence::getResidueDB_() const
	{
		return ResidueDB::getInstance();
	}	

	Size AASequence::getNumberOf(const String& residue) const
	{
		Size count(0);
		const Residue* res = getResidueDB_()->getResidue(residue);
		if (valid_)
		{
			for (vector<const Residue*>::const_iterator it = peptide_.begin(); it != peptide_.end(); ++it)
			{
				if (*it == res)
				{
					++count;
				}
			}
		}
		else
		{
			for (String::ConstIterator it = sequence_string_.begin(); it != sequence_string_.end(); ++it)
			{
				if (String(*it) == res->getOneLetterCode())
				{
					++count;
				}
			}
		}

		return count;
	}

	void AASequence::getAAFrequencies(Map<String, Size>& frequency_table) const
	{
		frequency_table.clear();
	
		if (valid_)
		{
			for (vector<const Residue*>::const_iterator it = peptide_.begin(); it != peptide_.end(); ++it)
			{
				frequency_table[(*it)->getOneLetterCode()] +=1;
			}
		}
		else
		{
			for (String::ConstIterator it = sequence_string_.begin(); it != sequence_string_.end(); ++it)
			{
				frequency_table[String(*it)] += 1;
			}
		}
	}

	AASequence::AASequence(ConstIterator begin, ConstIterator end)
		: valid_(true)
	{
		for (ConstIterator it = begin; it != end; ++it)
		{
			peptide_.push_back(&*it);
		}
	}

	void AASequence::setModification(Size index, const String& modification)
	{
		if (index >= peptide_.size())
		{
			if (valid_)
			{
				throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peptide_.size());
			}
			else
			{
				cerr << "AASequence: Warning: cannot set modification '" << modification 
						 << "' at position '" << index << "' because sequence '" << sequence_string_ 
						 << "' could not be parsed!" << endl;
				return;
			}
		}
		peptide_[index] = getResidueDB_()->getModifiedResidue(peptide_[index], modification);
	}

	void AASequence::setNTerminalModification(const String& modification)
	{
		if (modification == "")
		{
			n_term_mod_ = "";
			return;
		}
		set<String> mods(ModificationsDB::getInstance()->searchModifications(modification));
		if (mods.size() > 1)
		{
			cerr << "AASequence::setNTerminalModification: Error more than one modification with name '" 
					 << modification << "' found; using first one: " 
					 << *mods.begin() << "!" << endl;
		}
		if (mods.size() == 0)
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, modification);
		}
		// TODO check term specificity
		n_term_mod_ = ModificationsDB::getInstance()->getModification(*mods.begin()).getId();
	}

	void AASequence::setCTerminalModification(const String& modification)
	{
		if (modification == "")
		{
			c_term_mod_ = "";
			return;
		}
		set<String> mods = ModificationsDB::getInstance()->searchModifications(modification);
    if (mods.size() > 1)
    {
      cerr << "AASequence::setCTerminalModification: Error more than one modification with name '" 
					 << modification << "' found; using first one: "
					 << *mods.begin() << "!" << endl;
    }
    if (mods.size() == 0)
    {
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, modification);
    }
    // TODO check term specificity
		c_term_mod_ = ModificationsDB::getInstance()->getModification(*mods.begin()).getId();
	}

	const String& AASequence::getNTerminalModification() const
	{
		return n_term_mod_;
	}

	const String& AASequence::getCTerminalModification() const
	{
		return c_term_mod_;
	}

	bool AASequence::hasNTerminalModification() const
	{
		if (n_term_mod_ != "")
		{
			return true;
		}
		return false;
	}

	bool AASequence::hasCTerminalModification() const
	{
		if (c_term_mod_ != "")
		{
			return true;
		}
		return false;
	}

	bool AASequence::setStringSequence(const String& sequence)
	{
		c_term_mod_ = "";
		n_term_mod_ = "";
		parseString_(peptide_, sequence);
		return valid_;
	}
}
