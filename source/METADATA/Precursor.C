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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Precursor.h>

using namespace std;

namespace OpenMS
{

	const std::string Precursor::NamesOfActivationMethod[] = {"Unknown","Collision-induced dissociation","Post-source decay","Plasma desorption","Surface-induced dissociation","Blackbody infrared radiative dissociation","Electron capture dissociation","Infrared multiphoton dissociation","Sustained off-resonance irradiation","High-energy collision-induced dissociation","Low-energy collision-induced dissociation","Photodissociation","Electron transfer dissociation","Pulsed q dissociation"};

	Precursor::Precursor():
		RichPeak1D(),
		activation_method_(ACTMETHNULL),
		activation_energy_(0.0),
		window_low_(0.0),
		window_up_(0.0),
		charge_(0),
		possible_charge_states_()
	{
	}
	
	Precursor::Precursor(const Precursor& source):
		RichPeak1D(source),
	  activation_method_(source.activation_method_),
	  activation_energy_(source.activation_energy_),
	  window_low_(source.window_low_),
	  window_up_(source.window_up_),
		charge_(source.charge_),
		possible_charge_states_(source.possible_charge_states_)
	{
	}
	
	Precursor::~Precursor()
	{
	}
	
	Precursor& Precursor::operator = (const Precursor& source)
	{
	  if (&source == this) return *this;
	  
	  RichPeak1D::operator=(source);
	  activation_method_ = source.activation_method_;
	  activation_energy_ = source.activation_energy_;
	  window_low_ = source.window_low_;
	  window_up_ = source.window_up_;
		charge_ = source.charge_;
		possible_charge_states_ = source.possible_charge_states_;
		
	  return *this;
	}

  bool Precursor::operator== (const Precursor& rhs) const
  {
  	return 
	    activation_method_ == rhs.activation_method_ &&
	    activation_energy_ == rhs.activation_energy_ &&
	    window_low_ == rhs.window_low_ &&
	    window_up_ == rhs.window_up_ &&
			charge_ == rhs.charge_ &&
			possible_charge_states_ == rhs.possible_charge_states_ &&
			RichPeak1D::operator==(rhs)
 
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
	
	DoubleReal Precursor::getActivationEnergy() const 
	{
	  return activation_energy_; 
	}
	
	void Precursor::setActivationEnergy(DoubleReal activation_energy)
	{
	  activation_energy_ = activation_energy; 
	}
	
	DoubleReal Precursor::getIsolationWindowLowerOffset() const
	{
		return window_low_;
	}
	
	void Precursor::setIsolationWindowLowerOffset(DoubleReal bound)
	{
		window_low_ = bound;
	}
	
	DoubleReal Precursor::getIsolationWindowUpperOffset() const
	{
		return window_up_;
	}

	void Precursor::setIsolationWindowUpperOffset(DoubleReal bound)
	{
		window_up_ = bound;
	}

	Int Precursor::getCharge() const
	{
		return charge_;
	}

	void Precursor::setCharge( Int charge )
	{
		charge_ = charge;
		return;
	}

	std::vector<Int>& Precursor::getPossibleChargeStates()
	{
		return possible_charge_states_;
	}

	const std::vector<Int>& Precursor::getPossibleChargeStates() const
	{
		return possible_charge_states_;
	}

	void Precursor::setPossibleChargeStates(const std::vector<Int>& possible_charge_states)
	{
		possible_charge_states_ = possible_charge_states;
	}
	
}



