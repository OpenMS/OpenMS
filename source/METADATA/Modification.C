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

#include <OpenMS/METADATA/Modification.h>

using namespace std;

namespace OpenMS
{

	const std::string Modification::NamesOfSpecificityType[] = {"AA","AA_AT_CTERM","AA_AT_NTERM","CTERM","NTERM"};
	
	Modification::Modification() :
		SampleTreatment("Modification"),
		reagent_name_(""),
		mass_(0.0),
		specificity_type_(AA),
		affected_amino_acids_("")
	{
		
	}
	
	Modification::Modification(const Modification& source):
			SampleTreatment(source),
			reagent_name_(source.reagent_name_),
			mass_(source.mass_),
			specificity_type_(source.specificity_type_),
			affected_amino_acids_(source.affected_amino_acids_)
	{

	}
	
	Modification::~Modification()
	{

	}
	
	Modification& Modification::operator = (const Modification& source)
	{
	  if (&source == this) return *this;
	  
	  SampleTreatment::operator=(source);
		reagent_name_ = source.reagent_name_;
		mass_ = source.mass_;
		specificity_type_ = source.specificity_type_;
		affected_amino_acids_ = source.affected_amino_acids_;
	  
	  return *this;
	}

  bool Modification::operator== (const SampleTreatment& rhs) const
  {
  	if (type_!=rhs.getType()) return false;
  	
  	const Modification* tmp = dynamic_cast<const Modification*>(&rhs);
  	return
  		SampleTreatment::operator==(*tmp) &&
			reagent_name_ == tmp->reagent_name_ &&
			mass_ == tmp->mass_ &&
			specificity_type_ == tmp->specificity_type_ &&
			affected_amino_acids_ == tmp->affected_amino_acids_
  		;
  }

	SampleTreatment* Modification::clone() const
	{
		SampleTreatment* tmp = new Modification(*this);
		return tmp;
	}
	
	const String& Modification::getReagentName() const
	{
	  return reagent_name_;
	}
	
	void Modification::setReagentName(const String& reagent_name)
	{
	  reagent_name_ = reagent_name;
	}
	
	
	DoubleReal Modification::getMass() const
	{
	  return mass_;
	}
	
	void Modification::setMass(DoubleReal mass)
	{
	  mass_ = mass;
	}
	
	
	const Modification::SpecificityType& Modification::getSpecificityType() const
	{
	  return specificity_type_;
	}
	
	void Modification::setSpecificityType(const Modification::SpecificityType& specificity_type)
	{
	  specificity_type_ = specificity_type;
	}
	
	
	const String& Modification::getAffectedAminoAcids() const
	{
	  return affected_amino_acids_;
	}
	
	void Modification::setAffectedAminoAcids(const String& affected_amino_acids)
	{
	  affected_amino_acids_ = affected_amino_acids;
	}
	
}


