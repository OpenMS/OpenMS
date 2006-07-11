// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CHEMISTRY/PeptideSequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>

using namespace std;

namespace OpenMS 
{
	ResidueDB * PeptideSequence::custom_res_db_ = 0;

	PeptideSequence::PeptideSequence()
	{
	}

	PeptideSequence::PeptideSequence(const PeptideSequence& peptide_string)
		:	peptide_(peptide_string.peptide_)
	{
	}

	PeptideSequence::PeptideSequence(const String& peptide) throw(Exception::ParseError)
	{
		parseString_(peptide_, peptide);
	}

	PeptideSequence::PeptideSequence(ResidueDB* res_db_ptr)
	{
		custom_res_db_ = res_db_ptr;
	}

	PeptideSequence::~PeptideSequence()
	{
	}

	const Residue* PeptideSequence::getResidue(SignedInt index) const
		throw(Exception::IndexUnderflow, Exception::IndexOverflow)
	{
		if (index >= 0 && UnsignedInt(index) <= peptide_.size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peptide_.size());
		}
		if (index < 0)
		{
			throw Exception::IndexUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, 0);
		}
		return peptide_[index];
	}

	const Residue* PeptideSequence::getResidue(UnsignedInt index) const
		throw(Exception::IndexOverflow)
	{
		if (index <= peptide_.size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peptide_.size());
		}
		return peptide_[index];
	}
	
	EmpiricalFormula PeptideSequence::getFormula(Residue::ResidueType type, SignedInt charge) const
	{
		EmpiricalFormula ef;
		for (SignedInt i=0; i<charge; ++i)
		{
			ef += Formulas::H;
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
						return ef + Residue::getInternalToFull() - Residue::getBIonToFull() - Formulas::H;
					case Residue::AIon:
						return ef + Residue::getInternalToFull() - Residue::getAIonToFull() - Formulas::H;
					case Residue::CIon:
						return ef + Residue::getInternalToFull() - Formulas::OH + Formulas::NH - Formulas::H;
					case Residue::XIon:
						return ef + Residue::getInternalToFull() + Residue::getXIonToFull();
					case Residue::YIon:
						return ef + Residue::getInternalToFull() + Residue::getYIonToFull();
					case Residue::ZIon:
						return ef + Residue::getInternalToFull() - Residue::getZIonToFull();
					default:
						cerr << "PeptideSequence::getFormula: unknown ResidueType" << endl;
				}
			}			
		}
		return ef;
	}
	
	Real PeptideSequence::getAverageWeight(Residue::ResidueType type, SignedInt charge) const
	{
		return getFormula(type, charge).getAverageWeight();
	}

	Real PeptideSequence::getMonoWeight(Residue::ResidueType type, SignedInt charge) const
	{
		return getFormula(type, charge).getMonoWeight();
	}

	HashMap<const EmpiricalFormula*, Size> PeptideSequence::getNeutralLosses() const
	{
		// the following losses are from the Zhang paper (AC, 76, 14, 2004)
		// charge directed
		static const EmpiricalFormula R_44("NH2CHNH"); 
		static const EmpiricalFormula R_59("CN3H5"); // guanidium
		static const EmpiricalFormula R_61("N2H4CH");
		// charge remote
		static const EmpiricalFormula R_60("N2H4CO"); // combination of NH=C=NH + C-terminal H2O
		static const EmpiricalFormula H2O("H2O"); // loss from the C-terminus
		static const EmpiricalFormula NH3("NH3");
		HashMap<const EmpiricalFormula*, Size> losses;
		for (Size i=0;i!=peptide_.size();++i)
		{
			if (peptide_[i]->hasNeutralLoss())
			{
				const EmpiricalFormula* loss = &peptide_[i]->getLossFormula();
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
	}

	const Residue* PeptideSequence::operator [] (SignedInt index) const
		throw(Exception::IndexUnderflow, Exception::IndexOverflow)
	{
		if (index < 0)
		{
			throw Exception::IndexUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, 0);
		}
		else
		{
			if (UnsignedInt(index) >= size())
			{
				throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
			}
		}
		return peptide_[UnsignedInt(index)];
	}
	
	const Residue* PeptideSequence::operator [] (UnsignedInt index) const
		throw(Exception::IndexOverflow)
	{
		if (index >= size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
		}
		return peptide_[index];
	}
	
	PeptideSequence PeptideSequence::operator + (const PeptideSequence& sequence) const
	{
		PeptideSequence seq;
		seq.peptide_ = peptide_;
		for (Size i=0;i!=sequence.peptide_.size();++i)
		{
			seq.peptide_.push_back(sequence.peptide_[i]);
		}
		return seq;
	}
	
	PeptideSequence PeptideSequence::operator + (const String& peptide) const
		throw(Exception::ParseError)
	{
		PeptideSequence seq;
		seq.peptide_ = peptide_;
		vector<const Residue*> vec;
		parseString_(vec, peptide);
		for (Size i=0;i!=vec.size();++i)
		{
			seq.peptide_.push_back(vec[i]);
		}
		return seq;
	}
	
	PeptideSequence& PeptideSequence::operator += (const PeptideSequence& sequence)
	{
		for (Size i=0;i!=sequence.peptide_.size();++i)
		{
			peptide_.push_back(sequence.peptide_[i]);
		}
		return *this;
	}
	
	PeptideSequence& PeptideSequence::operator += (const String& peptide)
		throw(Exception::ParseError)
	{
		vector<const Residue*> vec;
		parseString_(vec, peptide);
		for (Size i=0;i!=vec.size();++i)
		{
			peptide_.push_back(vec[i]);
		}
		return *this;
	}

	void PeptideSequence::setResidueDB(ResidueDB* res_db)
	{
		custom_res_db_ = res_db;
	}

	Size PeptideSequence::size() const
	{
		return peptide_.size();
	}

	PeptideSequence PeptideSequence::getPrefix(Size index) const
		throw(Exception::IndexOverflow)
	{
		if (index > size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
		}
		PeptideSequence seq;
		for (Size i=0;i<index;++i)
		{
			seq.peptide_.push_back(peptide_[i]);
		}
		return seq;
	}

	PeptideSequence PeptideSequence::getSuffix(Size index) const
		throw(Exception::IndexOverflow)
	{
		if (index > size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
		}
		PeptideSequence seq;
		for (Size i=size()-index;i!=size();++i)
		{
			seq.peptide_.push_back(peptide_[i]);
		}
		return seq;
	}

	PeptideSequence PeptideSequence::getSubsequence(Size index, Size num) const
		throw(Exception::IndexOverflow)
	{
		if (index > size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
		}
		if (index+num >= size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index+num, size());
		}
		PeptideSequence seq;
		for (Size i=index;i!=num;++i)
		{
			seq.peptide_.push_back(peptide_[i]);
		}
		return seq;
	}
	
	bool PeptideSequence::has(const Residue* residue) const 
	{
		for (Size i=0;i!=peptide_.size();++i)
		{
			if (peptide_[i] == residue)
			{
				return true;
			}
		}
		return false;
	}

	bool PeptideSequence::has(const String& residue) const
	{
		if (!getResidueDB_()->hasResidue(residue))
		{
			return false;
		}
		return has(getResidueDB_()->getResidue(residue));
	}

	
	bool PeptideSequence::hasSubsequence(const PeptideSequence& sequence) const
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

	bool PeptideSequence::hasSubsequence(const String& sequence) const
		throw(Exception::ParseError)
	{
		vector<const Residue*> seq;
		parseString_(seq, sequence);
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
	
	bool PeptideSequence::hasPrefix(const PeptideSequence& sequence) const
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

	bool PeptideSequence::hasPrefix(const String& sequence) const
		throw(Exception::ParseError)
	{
		PeptideSequence seq(sequence);
		return hasPrefix(seq);
	}
	
	bool PeptideSequence::hasSuffix(const PeptideSequence& sequence) const
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

	bool PeptideSequence::hasSuffix(const String& sequence) const
		throw(Exception::ParseError)
	{
		PeptideSequence seq(sequence);
		return hasSuffix(seq);
	}
	
	bool PeptideSequence::operator == (const PeptideSequence& peptide) const
	{
		if (size() != peptide.size())
		{
			return false;
		}
		for (Size i=0;i!=size();++i)
		{
			if (peptide_[i] != peptide.peptide_[i])
			{
				return false;
			}
		}
		return true;
	}
	
	bool PeptideSequence::operator == (const String& peptide) const
		throw(Exception::ParseError)
	{
		vector<const Residue*> sequence;
		parseString_(sequence, peptide);
		if (size() != sequence.size())
		{
			return false;
		}
		for (Size i=0;i!=size();++i)
		{
			if (peptide_[i] != sequence[i])
			{
				return false;
			}	
		}
		return true;
	}

	bool PeptideSequence::operator != (const PeptideSequence& peptide) const
	{
		return !(*this == peptide);
	}

	bool PeptideSequence::operator != (const String& sequence) const
		throw(Exception::ParseError)
	{
		return !(*this == sequence);
	}

	ostream& operator << (ostream& os, const PeptideSequence& peptide)
	{
		for (Size i=0;i!=peptide.size();++i)
		{
			if (peptide.peptide_[i]->isModified())
			{
				if (peptide.peptide_[i]->getUnmodifiedName() != "")
				{
					os << peptide.peptide_[i]->getUnmodifiedName();
				}
				else
				{
					os << "[" << peptide.peptide_[i]->getMonoWeight() << "]";
				}
				if (peptide.peptide_[i]->getModification()->getShortName() != "")
				{
					os << "(" << peptide.peptide_[i]->getModification()->getShortName() << ")";
				}
				else
				{
					os << "([" << peptide.peptide_[i]->getModification()->getAddMonoWeight()  -
												peptide.peptide_[i]->getModification()->getDelMonoWeight() << "])";
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
		return os;
	}

	
	void PeptideSequence::parseString_(vector<const Residue*>& sequence, const String& peptide) const
		throw(Exception::ParseError)
	{
		if (peptide.size() > 0)
		{
			// split the peptide in its residues
			vector<String> split;
			if (!isalpha(peptide[0]) || !isupper(peptide[0]) || peptide[0] == '[')
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide, String(0));
			}
			else
			{
				Size pos(0);
				bool mod_open(false), tag_open(false);
				if (peptide[0] == '[')
				{
					tag_open = true;
				}
				for (Size i=1;i!=peptide.size();++i)
				{
					if (isalpha(peptide[i]) && isupper(peptide[i]) && !mod_open ||
							peptide[i] == '[' && !mod_open)
					{
						split.push_back(peptide.substr(pos, i-pos));
						pos = i;
					}
					if (peptide[i] == '(')
					{
						mod_open = true;
					}
					else
					{
						if (peptide[i] == ')')
						{
							mod_open = false;
						}
					}
					if (peptide[i] == '[')
					{
						tag_open = true;
					}
					else
					{
						if (peptide[i] == ']')
						{
							mod_open = false;
						}
					}
				}
				// push_back last residue
				split.push_back(peptide.substr(pos, peptide.size()-pos));
			}

			// parse the residues
			for (Size i=0;i!=split.size();++i)
			{
				String res = split[i];
				String name;
				String mod;
				String tag;
				for (Size j=0;j!=res.size();++j)
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
									throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, res, " '(' found but no ')'");
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
										throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, res, " '[' found but no ']'");
									}
									tag += res[k];
								}
								break;
							}
							else
							{
								throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, res, String(j));
							}
						}
					}
				}

				// now we have the name and the modification name (if there is one), or the tag
				const Residue * res_ptr = getResidueDB_()->getResidue(name);
				if (res_ptr == 0 && tag == "")
				{
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, name, "");
				}
				if (mod != "")
				{
					Size seq_size = sequence.size();
					set<const Residue*> mod_res = getResidueDB_()->getResidues(mod);
					for (set<const Residue*>::const_iterator it=mod_res.begin();it!=mod_res.end();++it)
					{
						if (getResidueDB_()->getResidue((*it)->getUnmodifiedName()) == res_ptr)
						{
							sequence.push_back(*it);
							break;
						}
					}
					if (seq_size == sequence.size())
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, res, mod);
					}
				}
				else
				{
					if (tag != "")
					{
						// if the residue db does not have this tag-residue, we add one
						if (res_ptr == 0)
						{
							Residue res(tag, String(""), String(""), EmpiricalFormula(""), EmpiricalFormula(""));
							res.setMonoWeight(tag.toFloat());
							res.setAverageWeight(tag.toFloat());
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
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, res, "");
					}
				}
			}
		}
	}

	ResidueDB* PeptideSequence::getResidueDB_() const
	{
		static ResidueDB * res_db = 0;
		if (res_db == 0)
		{
			res_db = new ResidueDB;
		}
		if (custom_res_db_ == 0)
		{
			return res_db;
		}
		else
		{
			return custom_res_db_;
		}
	}	
}
