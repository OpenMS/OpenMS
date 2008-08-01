// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Precursor.h>

using namespace std;

namespace OpenMS
{

	const std::string Precursor::NamesOfActivationMethod[] = {"Unknown","Collision-induced dissociation","Post-source decay","Plasma desorption","Surface-induced dissociation","Blackbody infrared radiative dissociation","Electron capture dissociation","Infrared multiphoton dissociation","Sustained off-resonance irradiation","High-energy collision-induced dissociation","Low-energy collision-induced dissociation","Photodissociation","Electron transfer dissociation","Pulsed q dissociation"};
	const std::string Precursor::NamesOfEnergyUnits[] = {"Unknown","Electron volt","Precent"};

	Precursor::Precursor():
		MetaInfoInterface(),
		activation_method_(ACTMETHNULL),
		activation_energy_(0.0),
		activation_energy_unit_(UNITSNULL),
		window_size_(0.0)
	{
		
	}
	
	Precursor::Precursor(const Precursor& source):
		MetaInfoInterface(source),
	  activation_method_(source.activation_method_),
	  activation_energy_(source.activation_energy_),
	  activation_energy_unit_(source.activation_energy_unit_),
	  window_size_(source.window_size_)
	{
	  
	}
	
	Precursor::~Precursor()
	{
	  
	}
	
	Precursor& Precursor::operator = (const Precursor& source)
	{
	  if (&source == this) return *this;
	  
	  MetaInfoInterface::operator=(source);
	  activation_method_ = source.activation_method_;
	  activation_energy_ = source.activation_energy_;
	  activation_energy_unit_ = source.activation_energy_unit_;
	  window_size_ = source.window_size_;

	  return *this;
	}

  bool Precursor::operator== (const Precursor& rhs) const
  {
  	return 
	    activation_method_ == rhs.activation_method_ &&
	    activation_energy_ == rhs.activation_energy_ &&
	    activation_energy_unit_ == rhs.activation_energy_unit_ &&
	    window_size_ == rhs.window_size_ &&
  		MetaInfoInterface::operator==(rhs)
  		;
  }
  
  bool Precursor::operator!= (const Precursor& rhs) const
  {
  	return !(operator==(rhs));
 	}
	
	Precursor::ActivationMethod Precursor::getActivationMethod() const 
	{
	  return activation_method_; 
	}
	
	void Precursor::setActivationMethod(Precursor::ActivationMethod activation_method)
	{
	  activation_method_ = activation_method; 
	}
	
	float Precursor::getActivationEnergy() const 
	{
	  return activation_energy_; 
	}
	
	void Precursor::setActivationEnergy(float activation_energy)
	{
	  activation_energy_ = activation_energy; 
	}
	
	Precursor::EnergyUnits Precursor::getActivationEnergyUnit() const 
	{
	  return activation_energy_unit_; 
	}
	
	void Precursor::setActivationEnergyUnit(Precursor::EnergyUnits activation_energy_unit)
	{
	  activation_energy_unit_ = activation_energy_unit; 
	}
	
	float Precursor::getWindowSize() const
	{
		return window_size_;
	}
	
	void Precursor::setWindowSize(float size)
	{
		window_size_ = size;
	}

}



