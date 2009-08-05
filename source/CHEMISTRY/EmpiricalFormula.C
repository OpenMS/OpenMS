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

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <iostream>

using namespace std;

namespace OpenMS
{
	EmpiricalFormula::EmpiricalFormula()
		:	charge_(0),
			element_db_(ElementDB::getInstance())
	{
	}

	EmpiricalFormula::EmpiricalFormula(const EmpiricalFormula& formula)
		:	formula_(formula.formula_),
			charge_(formula.charge_),
			element_db_(formula.element_db_)
	{
	}

	EmpiricalFormula::EmpiricalFormula(const String& formula) 
		:	element_db_(ElementDB::getInstance())
	{
		charge_ = parseFormula_(formula_, formula);
	}

	EmpiricalFormula::EmpiricalFormula(SignedSize number, const Element* element, Size charge)
	{
		formula_[element] = number;
		charge_ = charge;
	}

	EmpiricalFormula::~EmpiricalFormula()
	{
	}
	
	DoubleReal EmpiricalFormula::getMonoWeight() const
	{
		DoubleReal weight(0);
		weight += Constants::PROTON_MASS_U * charge_;
		Map<const Element*, SignedSize>::ConstIterator it=formula_.begin();
		for (; it != formula_.end(); ++it)
		{
			weight += it->first->getMonoWeight() * (DoubleReal)it->second;
		}
		return weight;
	}

	DoubleReal EmpiricalFormula::getAverageWeight() const 
	{
		DoubleReal weight(0);
		weight += Constants::PROTON_MASS_U * charge_;
		Map<const Element*, SignedSize>::ConstIterator it=formula_.begin();
		for (; it != formula_.end(); ++it)
		{
			weight += it->first->getAverageWeight() * (DoubleReal)it->second;
		}
		return weight;
	}

	IsotopeDistribution EmpiricalFormula::getIsotopeDistribution(UInt max_depth) const
	{
		IsotopeDistribution result(max_depth);
		Map<const Element*, SignedSize>::ConstIterator it=formula_.begin();
		for (; it!=formula_.end(); ++it)
		{
			IsotopeDistribution tmp = it->first->getIsotopeDistribution(); 
			tmp.setMaxIsotope(max_depth);
			result += tmp * it->second;
		}
		result.renormalize();
		return result;
	}
	
	const Element* EmpiricalFormula::getElement(const String& name) const
	{
		return element_db_->getElement(name);
	}

	const Element* EmpiricalFormula::getElement(UInt atomic_number) const
	{
		return element_db_->getElement(atomic_number);
	}
	
	const ElementDB* EmpiricalFormula::getElementDB() const
	{
		return element_db_;
	}

	Size EmpiricalFormula::getNumberOf(UInt atomic_number) const
	{
		if (element_db_->hasElement(atomic_number))
		{
			if (formula_.has(element_db_->getElement(atomic_number)))
			{
				return formula_[element_db_->getElement(atomic_number)];
			}
		}
		return 0;
	}

	Size EmpiricalFormula::getNumberOf(const String& name) const
	{
		if (element_db_->hasElement(name))
		{
			if (formula_.has(element_db_->getElement(name)))
			{
				return formula_[element_db_->getElement(name)];
			}
		}
		return 0;
	}

	Size EmpiricalFormula::getNumberOf(const Element* element) const
	{
		if (formula_.has(element))
		{
			return formula_[element];
		}
		return 0;
	}
	
	Size EmpiricalFormula::getNumberOfAtoms() const
	{
		Size num_atoms(0);
		Map<const Element*, SignedSize>::ConstIterator it = formula_.begin();
		for (; it!=formula_.end(); ++it)
		{
			num_atoms += it->second;
		}
		return num_atoms;
	}
	
	void EmpiricalFormula::setCharge(Size charge)
	{
		charge_ = charge;
	}
	
	Size EmpiricalFormula::getCharge() const
	{
		return charge_;
	}

	String EmpiricalFormula::getString() const
	{
		String formula;
		Map<String, SignedSize> new_formula;
		
		for (Map<const Element*, SignedSize>::ConstIterator it = formula_.begin(); it!=formula_.end(); ++it)
		{
			new_formula[it->first->getSymbol()] = it->second;
		}

		for (Map<String, SignedSize>::ConstIterator it = new_formula.begin(); it != new_formula.end(); ++it)
		{						
			formula += it->first + String(it->second);
		}
		return formula;
	}

	EmpiricalFormula& EmpiricalFormula::operator = (const EmpiricalFormula& formula)
	{
		if (this != &formula)
		{
			formula_ = formula.formula_;
			charge_ = formula.charge_;
		}
		return *this;
	}

	EmpiricalFormula& EmpiricalFormula::operator = (const String& formula)
	{
		charge_ = 0;
		formula_.clear();
		charge_ = parseFormula_(formula_, formula);
		return *this;
	}

	EmpiricalFormula EmpiricalFormula::operator * (const SignedSize& times) const
	{
		EmpiricalFormula ef(*this);
		Map<const Element*, SignedSize>::ConstIterator it=formula_.begin();
		for (;it!=formula_.end();++it)
		{
				ef.formula_[it->first] *= times;
		}
		ef.charge_ *= times;
		return ef;
	}	
	
	EmpiricalFormula EmpiricalFormula::operator + (const EmpiricalFormula& formula) const
	{
		EmpiricalFormula ef;
		ef.formula_ = formula.formula_;
		Map<const Element*, SignedSize>::ConstIterator it=formula_.begin();
		for (;it!=formula_.end();++it)
		{
			if (ef.formula_.has(it->first))
			{
				ef.formula_[it->first] += it->second;
			}
			else
			{
				ef.formula_[it->first] = it->second;
			}
		}
		ef.charge_ = charge_ + formula.charge_;
		return ef;
	}

	EmpiricalFormula EmpiricalFormula::operator + (const String& formula) const
	{
		EmpiricalFormula ef;
		Size charge = parseFormula_(ef.formula_, formula);
		Map<const Element*, SignedSize>::ConstIterator it=formula_.begin();
		for (;it!=formula_.end();++it)
		{
			if (ef.formula_.has(it->first))
			{
				ef.formula_[it->first] += it->second;
			}
			else
			{
				ef.formula_[it->first] = it->second;
			}
		}
		ef.charge_ = charge_ + charge;
		return ef;
	}
	
	EmpiricalFormula& EmpiricalFormula::operator += (const EmpiricalFormula& formula)
	{
		Map<const Element*, SignedSize>::ConstIterator it=formula.formula_.begin();
		for (;it!=formula.formula_.end();++it)
		{
			if (formula_.has(it->first))
			{			
				formula_[it->first] += it->second;
			}
			else
			{
				formula_[it->first] = it->second;
			}
		}
		charge_ += formula.charge_;
		return *this;
	}

	EmpiricalFormula& EmpiricalFormula::operator += (const String& formula) 
	{
		Map<const Element*, SignedSize> str_formula;
		Size charge = parseFormula_(str_formula, formula);
		charge_ += charge;
		Map<const Element*, SignedSize>::ConstIterator it;
		for (it=str_formula.begin(); it!=str_formula.end(); ++it)
		{
			if (formula_.has(it->first))
			{
				formula_[it->first] += it->second;
			}
			else
			{
				formula_[it->first] = it->second;
			}
		}
		return *this;
	}
	
	EmpiricalFormula EmpiricalFormula::operator - (const EmpiricalFormula& formula) const
	{
		EmpiricalFormula ef;
		Map<const Element*, SignedSize>::ConstIterator it=formula.formula_.begin();
		for (; it!=formula.formula_.end(); ++it)
		{
			const Element * e = it->first;
			SignedSize num = it->second;
			if (formula_.has(e))
			{
				if (formula_[e] - num != 0)
				{
					ef.formula_[e] = formula_[e] - num;
				}
			}
			else
			{
				ef.formula_[e] = -num;
			}
		}

		for (it=formula_.begin();it!=formula_.end();++it)
		{
			if (!ef.hasElement(it->first) && !formula.formula_.has(it->first))
			{
				ef.formula_[it->first] = it->second;
			}
		}
		
		ef.charge_ = charge_ - formula.charge_;
		return ef;
	}

	EmpiricalFormula EmpiricalFormula::operator - (const String& formula) const
	{
		EmpiricalFormula ef;
		Map<const Element*, SignedSize> str_formula;
		Size charge = parseFormula_(str_formula, formula);
		ef.charge_ = charge_ - charge;
		Map<const Element*, SignedSize>::ConstIterator it=str_formula.begin();
		for (; it!=str_formula.end(); ++it)
		{
			const Element * e = it->first;
			SignedSize num = it->second;
			if (formula_.has(e))
			{
				if (formula_[e] - num != 0)
				{
					ef.formula_[e] = formula_[e] - num;
				}
			}
			else
			{
				ef.formula_[e] = -num;
			}
		}

		for (it=formula_.begin();it!=formula_.end();++it)
		{
			if (!ef.hasElement(it->first) && !str_formula.has(it->first))
			{
				ef.formula_[it->first] = it->second;
			}
		}
		
		return ef;
	}
	
	EmpiricalFormula& EmpiricalFormula::operator -= (const EmpiricalFormula& formula) 
	{
		Map<const Element*, SignedSize>::ConstIterator it=formula.formula_.begin();
		for (; it!=formula.formula_.end(); ++it)
		{
			if (formula_.has(it->first))
			{
				if (formula_[it->first] != it->second)
				{
					formula_[it->first] -= it->second;
				}
				else
				{
					formula_.erase(it->first);
				}
			}
			else
			{
				formula_[it->first] = -it->second;
			}
		}
		charge_ -= formula.charge_;
		return *this;
	}

	EmpiricalFormula& EmpiricalFormula::operator -= (const String& formula) 
	{
		Map<const Element*, SignedSize> str_formula;
		Size charge = parseFormula_(str_formula, formula);
		charge_ -= charge;
		Map<const Element*, SignedSize>::ConstIterator it = str_formula.begin();
		for (; it!=str_formula.end(); ++it)
		{
			if (formula_.has(it->first))
			{
				if (formula_[it->first] != it->second)
				{
					formula_[it->first] -= it->second;
				}
				else
				{
					formula_.erase(it->first);
				}
			}
			else
			{
				formula_[it->first] = -it->second;
			}
		}
		return *this;
	}

	bool EmpiricalFormula::isCharged() const
	{
		return charge_!=0;
	}
	
	bool EmpiricalFormula::isEmpty() const
	{
		return (formula_.size() == 0);
	}

	bool EmpiricalFormula::hasElement(const Element* element) const
	{
		return formula_.has(element);
	}

	bool EmpiricalFormula::hasElement(const String& element) const
	{
		if (!element_db_->hasElement(element))
		{
			return false;
		}
		else
		{
			if (formula_.has(element_db_->getElement(element)))
			{
				return true;
			}
		}
		return false;
	}

	bool EmpiricalFormula::hasElement(UInt atomic_number) const
	{
		if (!element_db_->hasElement(atomic_number))
		{
			return false;
		}
		else
		{
			if (formula_.has(element_db_->getElement(atomic_number)))
			{
				return true;
			}
		}
		return false;
	}
	
	bool EmpiricalFormula::operator == (const EmpiricalFormula& formula) const
	{
		return (formula_ == formula.formula_ && charge_ == formula.charge_);
	}
	
	bool EmpiricalFormula::operator == (const String& formula) const
	{
		Map<const Element*, SignedSize> str_formula;
		Size charge = parseFormula_(str_formula, formula);
		return (formula_ == str_formula && charge_ == charge);
	}
	
	bool EmpiricalFormula::operator != (const EmpiricalFormula& formula) const
	{
		return (formula_ != formula.formula_ || charge_ != formula.charge_);
	}
	
	bool EmpiricalFormula::operator != (const String& formula) const
	{
		Map<const Element *, SignedSize> str_formula;
		Size charge = parseFormula_(str_formula, formula);
		return (formula_ != str_formula || charge_ != charge);
	}
	
	ostream& operator << (ostream& os, const EmpiricalFormula& formula)
	{
		Map<String, SignedSize> new_formula;
		for (Map<const Element*, SignedSize>::ConstIterator it = formula.formula_.begin(); it != formula.formula_.end(); ++it)
		{
			new_formula[it->first->getSymbol()] = it->second;
		}

		for (Map<String, SignedSize>::ConstIterator it = new_formula.begin(); it != new_formula.end(); ++it)
		{
			os << it->first;
			if (it->second > 1)
			{
				os << it->second;
			}
		}
		if (formula.charge_ != 0)
		{
			if (formula.charge_ == 1)
			{
				os << "+";
			}
			else
			{
				if (formula.charge_ == 1)
				{
					os << "+";
				}
				else
				{
					os << "+" << formula.charge_;
				}
			}
		}
		return os;
	}

	Size EmpiricalFormula::parseFormula_(Map<const Element*, SignedSize>& ef, const String& formula) const 
	{
		Size charge = 0;
		
		String symbol;
		String number;

		// split the formula
		vector<String> splitter;	
		if (formula.size() > 0)
		{
			if (!isdigit(formula[0]) || formula[0] == '(')
			{
				bool is_isotope(false), is_symbol(false);
				String split;
				for (Size i=0; i<formula.size(); ++i)
				{
					if ((isupper(formula[i]) && (!is_isotope || is_symbol))
							|| formula[i] == '(' 
							|| formula[i] == '+')
					{
						if (split != "")
						{
							splitter.push_back(split);
							is_isotope = false;
							is_symbol = false;
						}
						split = String(1, formula[i]);
					}
					else
					{
						split += String(1, formula[i]);
					}
					if (formula[i] == '(')
          {
            is_isotope = true;
          }
					if (isupper(formula[i]))
					{
						is_symbol = true;
					}
				}
				splitter.push_back(split);
			}
			else
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, formula, "This formula does not begin with an element!");
			}
		}

		// add up the elements and charge
		for (Size i=0;i!=splitter.size();++i)
		{
			//cerr << "splitter[" << i << "]='" << splitter[i] << "'" << endl;
			String split = splitter[i];
			String number;
			String symbol;
			bool had_symbol(false);
			for (SignedSize j = split.size() - 1; j >= 0; --j)
			{
				if (!had_symbol && (isdigit(split[j]) || split[j] == '-'))
				{
					number = split[j] + number;
				}
				else
				{
					symbol = split[j] + symbol;
					had_symbol = true;
				}
			}

			SignedSize num(1);
			if (number != "")
			{
				num = number.toInt();
			}
			
			if (element_db_->hasElement(symbol))
			{
				if (num != 0)
				{
					const Element* e = element_db_->getElement(symbol);
					if (ef.has(e))
					{
						ef[e] += num;
					}
					else
					{
						ef[e] = num;
					}
				}
			}
			else
			{
				if (symbol == "+")
				{
					charge += num;
				}
				else
				{
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, split, symbol);
				}
			}
		}
			
		return charge;
	}
}

