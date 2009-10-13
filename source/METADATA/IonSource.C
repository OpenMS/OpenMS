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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/IonSource.h>

using namespace std;

namespace OpenMS
{

	const std::string IonSource::NamesOfInletType[] = {"Unknown","Direct","Batch","Chromatography","Particle beam","Membrane sparator","Open split","Jet separator","Septum","Reservoir","Moving belt","Moving wire","Flow injection analysis","Electro spray","Thermo spray","Infusion","Continuous flow fast atom bombardment","Inductively coupled plasma","Membrane inlet","Nanospray inlet"};
	
	const std::string IonSource::NamesOfIonizationMethod[] = {"Unknown","Electrospray ionisation","Electron ionization","Chemical ionisation","Fast atom bombardment","Thermospray","Laser desorption","Field desorption","Flame ionization","Plasma desorption","Secondary ion MS","Thermal ionization","Atmospheric pressure ionisation","ISI","Collsion induced decomposition","Collsiona activated decomposition","HN","Atmospheric pressure chemical ionization","Atmospheric pressure photo ionization","Inductively coupled plasma","Nano electrospray ionization","Micro electrospray ionization","Surface enhanced laser desorption ionization","Surface enhanced neat desorption","Fast ion bombardment","Matrix-assisted laser desorption ionization","Multiphoton ionization","Desorption ionization","Flowing afterglow","Field ionization","Glow discharge ionization","Negative ion chemical ionization","Neutralization reionization mass spectrometry","Photoionization","Pyrolysis mass spectrometry","Resonance enhanced multiphoton ionization","Adiabatic ionization","Associative ionization","Autodetachment","Autoionization","Charge exchange ionization","Chemi-ionization","Dissociative ionization","Liquid secondary ionization","Penning ionization","Soft ionization","Spark ionization","Surface ionization","Vertical ionization","Atmospheric pressure matrix-assisted laser desorption ionization","Desorption/ionization on silicon","Surface-assisted laser desorption ionization"};

	const std::string IonSource::NamesOfPolarity[] = {"Unknown","Positive","Negative"};

	IonSource::IonSource():
		MetaInfoInterface(),
		inlet_type_(INLETNULL),
		ionization_method_(IONMETHODNULL),
		polarity_(POLNULL),
		order_(0)
	{
	}
	
	IonSource::IonSource(const IonSource& source):
		MetaInfoInterface(source),
	  inlet_type_(source.inlet_type_),
	  ionization_method_(source.ionization_method_),
	  polarity_(source.polarity_),
		order_(source.order_)
	{
	}
	
	IonSource::~IonSource()
	{
	}
	
	IonSource& IonSource::operator = (const IonSource& source)
	{
	  if (&source == this) return *this;
	  
	  order_ = source.order_;
	  inlet_type_ = source.inlet_type_;
	  ionization_method_ = source.ionization_method_;
	  polarity_ = source.polarity_;
	  MetaInfoInterface::operator=(source);
	  
	  return *this;
	}

  bool IonSource::operator== (const IonSource& rhs) const
  {
  	return
	 		order_ == rhs.order_ &&
		  inlet_type_ == rhs.inlet_type_ &&
		  ionization_method_ == rhs.ionization_method_ &&
		  polarity_ == rhs.polarity_ &&
  		MetaInfoInterface::operator==(rhs)
  		;
  }
  
  bool IonSource::operator!= (const IonSource& rhs) const
  {
  	return !(operator==(rhs));
 	}
	
	IonSource::InletType IonSource::getInletType() const 
	{
	  return inlet_type_; 
	}
	
	void IonSource::setInletType(IonSource::InletType inlet_type)
	{
	  inlet_type_ = inlet_type; 
	}
	
	IonSource::IonizationMethod IonSource::getIonizationMethod() const 
	{
	  return ionization_method_; 
	}
	
	void IonSource::setIonizationMethod(IonSource::IonizationMethod ionization_type)
	{
	  ionization_method_ = ionization_type; 
	}
	
	IonSource::Polarity IonSource::getPolarity() const 
	{
	  return polarity_; 
	}
	
	void IonSource::setPolarity(IonSource::Polarity polarity)
	{
	  polarity_ = polarity; 
	}
	
	Int IonSource::getOrder() const
  {
  	return order_;
  }
  
  void IonSource::setOrder(Int order)
  {
  	order_ = order;
  }
  
}

