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
// $Maintainer: Mathias Walzer $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Precursor.h>

using namespace std;

namespace OpenMS
{

	const std::string Precursor::NamesOfActivationMethod[] = {"Collision-induced dissociation","Post-source decay","Plasma desorption","Surface-induced dissociation","Blackbody infrared radiative dissociation","Electron capture dissociation","Infrared multiphoton dissociation","Sustained off-resonance irradiation","High-energy collision-induced dissociation","Low-energy collision-induced dissociation","Photodissociation","Electron transfer dissociation","Pulsed q dissociation"};

	Precursor::Precursor():
		CVTermList(),
		Peak1D(),
		activation_methods_(),
		activation_energy_(0.0),
		window_low_(0.0),
		window_up_(0.0),
		charge_(0),
		possible_charge_states_()
	{
	}
	
	Precursor::Precursor(const Precursor& source):
		CVTermList(source),
		Peak1D(source),
	  activation_methods_(source.activation_methods_),
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
	 
		CVTermList::operator=(source); 
	  Peak1D::operator=(source);
	  activation_methods_ = source.activation_methods_;
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
	    activation_methods_ == rhs.activation_methods_ &&
	    activation_energy_ == rhs.activation_energy_ &&
	    window_low_ == rhs.window_low_ &&
	    window_up_ == rhs.window_up_ &&
			charge_ == rhs.charge_ &&
			possible_charge_states_ == rhs.possible_charge_states_ &&
			Peak1D::operator==(rhs) &&
			CVTermList::operator == (rhs)
 
 		;
  }
  
  bool Precursor::operator!= (const Precursor& rhs) const
  {
  	return !(operator==(rhs));
 	}
	
	const set<Precursor::ActivationMethod>& Precursor::getActivationMethods() const 
	{
	  return activation_methods_;
	}

	set<Precursor::ActivationMethod>& Precursor::getActivationMethods() 
	{
	  return activation_methods_;
	}
	
	void Precursor::setActivationMethods(const set<Precursor::ActivationMethod>& activation_methods)
	{
	  activation_methods_ = activation_methods; 
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



