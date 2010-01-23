// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/Element.h>

using namespace std;

namespace OpenMS
{
	Element::Element()
		:	name_(OPENMS_CHEMISTRY_ELEMENT_NAME_DEFAULT),
			symbol_(OPENMS_CHEMISTRY_ELEMENT_SYMBOL_DEFAULT),
			atomic_number_(OPENMS_CHEMISTRY_ELEMENT_ATOMICNUMBER_DEFAULT),
			average_weight_(OPENMS_CHEMISTRY_ELEMENT_WEIGHT_DEFAULT),
			mono_weight_(OPENMS_CHEMISTRY_ELEMENT_WEIGHT_DEFAULT)
	{
	}

	Element::Element(const Element& e)
		:	name_(e.name_),
			symbol_(e.symbol_),
			atomic_number_(e.atomic_number_),
			average_weight_(e.average_weight_),
			mono_weight_(e.mono_weight_),
			isotopes_(e.isotopes_)
	{
	}
	
	Element::Element( const String& name,
										const String& symbol,
										UInt atomic_number,
										DoubleReal average_weight,
										DoubleReal mono_weight,
										const IsotopeDistribution& isotopes)
		:	name_(name),
			symbol_(symbol),
			atomic_number_(atomic_number),
			average_weight_(average_weight),
			mono_weight_(mono_weight),
			isotopes_(isotopes)
	{
	}

	Element::~Element()
	{
	}

	void Element::setAtomicNumber(UInt atomic_number)
	{
		atomic_number_ = atomic_number;
	}	
	
	UInt Element::getAtomicNumber() const
	{
		return atomic_number_;
	}
	
	void Element::setAverageWeight(DoubleReal weight)
	{
		average_weight_ = weight;
	}
	
	DoubleReal Element::getAverageWeight() const
	{
		return average_weight_;
	}
	
	void Element::setMonoWeight(DoubleReal weight)
	{
		mono_weight_ = weight;
	}
	
	DoubleReal Element::getMonoWeight() const
	{
		return mono_weight_;
	}
	
	void Element::setIsotopeDistribution(const IsotopeDistribution& distribution)
	{
		isotopes_ = distribution;
	}
	
	const IsotopeDistribution& Element::getIsotopeDistribution() const
	{
		return isotopes_;
	}

	void Element::setName(const String& name)
	{
		name_ = name;
	}
	
	const String& Element::getName() const
	{
		return name_;
	}

	void Element::setSymbol(const String& symbol)
	{
		symbol_ = symbol;
	}
	
	const String& Element::getSymbol() const
	{
		return symbol_;
	}

	Element& Element::operator = (const Element& element)
	{
		name_ = element.name_;
		symbol_ = element.symbol_;
		atomic_number_ = element.atomic_number_;
		average_weight_ = element.average_weight_;
		mono_weight_ = element.mono_weight_;
		isotopes_ = element.isotopes_;
		return *this;
	}

	bool Element::operator == (const Element& element) const
	{
		return 	name_ == element.name_ &&
						symbol_ == element.symbol_ &&
						atomic_number_ == element.atomic_number_ &&
						average_weight_ == element.average_weight_ &&
						mono_weight_ == element.mono_weight_ &&
						isotopes_ == element.isotopes_;
	}

	bool Element::operator != (const Element& element) const
	{
		return !(*this == element);
	}
	
	std::ostream& operator << (std::ostream& os, const Element& element)
	{
		os 	<< element.name_ << " " 
				<< element.symbol_ << " " 
				<< element.atomic_number_ << " "
				<< element.average_weight_ << " "
				<< element.mono_weight_;

		for (IsotopeDistribution::ConstIterator it=element.isotopes_.begin(); it!=element.isotopes_.end(); ++it)
		{
			if (it->second > 0.0f)
			{
				os << " " << it->first << "=" << it->second*100 << "%";
			}
		}
		return os;
	}
} // namespace OpenMS

