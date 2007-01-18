// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/METADATA/IonSource.h>

using namespace std;

namespace OpenMS
{

	const std::string IonSource::NamesOfInletType[] = {"Unknown","DIRECT","BATCH","CHROMATOGRAPHY","PARTICLEBEAM","MEMBRANESEPARATOR","OPENSPLIT","JETSEPARATOR","SEPTUM","RESERVOIR","MOVINGBELT","MOVINGWIRE","FLOWINJECTIONANALYSIS","ELECTROSPRAYINLET","THERMOSPRAYINLET","INFUSION","CONTINUOUSFLOWFASTATOMBOMBARDMENT","INDUCTIVELYCOUPLEDPLASMA"};
	const std::string IonSource::NamesOfIonizationMethod[] = {"Unknown","ESI (electrospray ionisation)","EI (electron impact)","CI (chemical ionisation)","FAB (fast atom bombardment)","TSP (thermospray)","LD (laser desorption)","FD (field desorption)","FI","PD (plasma desorption)","SI (secondary ion MS)","TI","API (atmospheric pressure ionisation)","ISI","CID (collsion induced decomposition)","CAD (collsiona activated decomposition)","HN","APCI (atmospheric pressure chemical ionization)","APPI","ICP (inductively coupled plasma)"};
	const std::string IonSource::NamesOfPolarity[] = {"Unknown","POSITIVE","NEGATIVE"};

	IonSource::IonSource():
		MetaInfoInterface(),
		inlet_type_(INLETNULL),
		ionization_method_(IONMETHODNULL),
		polarity_(POLNULL)
	{
		
	}
	
	IonSource::IonSource(const IonSource& source):
		MetaInfoInterface(source),
	  inlet_type_(source.inlet_type_),
	  ionization_method_(source.ionization_method_),
	  polarity_(source.polarity_)
	{
	  
	}
	
	IonSource::~IonSource()
	{
	  
	}
	
	IonSource& IonSource::operator = (const IonSource& source)
	{
	  if (&source == this) return *this;
	  
	  inlet_type_ = source.inlet_type_;
	  ionization_method_ = source.ionization_method_;
	  polarity_ = source.polarity_;
	  MetaInfoInterface::operator=(source);
	  
	  return *this;
	}

  bool IonSource::operator== (const IonSource& rhs) const
  {
  	return 
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
	
}

