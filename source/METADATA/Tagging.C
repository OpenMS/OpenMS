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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Tagging.h>

using namespace std;

namespace OpenMS
{
	const std::string Tagging::NamesOfIsotopeVariant[] = {"LIGHT","HEAVY"};

	Tagging::Tagging() :
		Modification(),
		mass_shift_(0.0),
		variant_(LIGHT)
	{
		type_= "Tagging";
	}
	
	Tagging::Tagging(const Tagging& source):
			Modification(source),
			mass_shift_(source.mass_shift_),
			variant_(source.variant_)
	{
	  //????
	}
	
	Tagging::~Tagging()
	{
	  //????
	}
	
	Tagging& Tagging::operator = (const Tagging& source)
	{
	  if (&source == this) return *this;
	  
	  Modification::operator=(source);
		mass_shift_ = source.mass_shift_;
		variant_ = source.variant_;
	  
	  return *this;
	}
	
  bool Tagging::operator== (const SampleTreatment& rhs) const
  {
  	if (type_!=rhs.getType()) return false;
  	
  	const Tagging* tmp = dynamic_cast<const Tagging*>(&rhs);
  	return
  		Modification::operator==(rhs) &&
			mass_shift_ == tmp->mass_shift_ &&
			variant_ == tmp->variant_
  		;
  }	
	
	SampleTreatment* Tagging::clone() const
	{
		SampleTreatment* tmp = new Tagging(*this);
		return tmp;
	}
	
	
	DoubleReal Tagging::getMassShift() const
	{
	  return mass_shift_;
	}
	
	void Tagging::setMassShift(DoubleReal mass_shift)
	{
	  mass_shift_ = mass_shift;
	}
	
	
	const Tagging::IsotopeVariant& Tagging::getVariant() const
	{
	  return variant_;
	}
	
	void Tagging::setVariant(const Tagging::IsotopeVariant& variant)
	{
	  variant_ = variant;
	}
	
}


