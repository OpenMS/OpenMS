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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/IonDetector.h>

using namespace std;

namespace OpenMS
{

	const std::string IonDetector::NamesOfType[] = {"Unknown","Electron multiplier","Photo multiplier","Focal plane array","Faraday cup","Conversion dynode electron multiplier","Conversion dynode photo multiplier","Multi-collector","Channel electron multiplier","channeltron","daly detector","microchannel plate detector","array detector","conversion dynode","dynode","focal plane collector","ion-to-photon detector","point collector","postacceleration detector","photodiode array detector","inductive detector","electron multiplier tube"};

	const std::string IonDetector::NamesOfAcquisitionMode[] = {"Unknown","Pulse counting","Analog-digital converter","Time-digital converter","Transient recorder"};
	
	IonDetector::IonDetector():
		MetaInfoInterface(),
		type_(TYPENULL),
		acquisition_mode_(ACQMODENULL),
		resolution_(0.0),
		ADC_sampling_frequency_(0.0),
		order_(0)
	{
		
	}
	
	IonDetector::IonDetector(const IonDetector& source):
		MetaInfoInterface(source),
	  type_(source.type_),
	  acquisition_mode_(source.acquisition_mode_),
	  resolution_(source.resolution_),
	  ADC_sampling_frequency_(source.ADC_sampling_frequency_),
		order_(source.order_)
	{
	  
	}
	
	IonDetector::~IonDetector()
	{
	  
	}
	
	IonDetector& IonDetector::operator = (const IonDetector& source)
	{
	  if (&source == this) return *this;
	  
	  order_ = source.order_;
	  type_ = source.type_;
	  acquisition_mode_ = source.acquisition_mode_;
	  resolution_ = source.resolution_;
	  ADC_sampling_frequency_ = source.ADC_sampling_frequency_;
	  MetaInfoInterface::operator=(source);
	  
	  return *this;
	}

  bool IonDetector::operator== (const IonDetector& rhs) const
  {
  	return
	 		order_ == rhs.order_ &&
		  type_ == rhs.type_ &&
		  acquisition_mode_ == rhs.acquisition_mode_ &&
		  resolution_ == rhs.resolution_ &&
		  ADC_sampling_frequency_ == rhs.ADC_sampling_frequency_ &&
  		MetaInfoInterface::operator==(rhs)
  		;
  }
  
  bool IonDetector::operator!= (const IonDetector& rhs) const
  {
  	return !(operator==(rhs));
 	}
	
	IonDetector::Type IonDetector::getType() const 
	{
	  return type_; 
	}
	
	void IonDetector::setType(IonDetector::Type type)
	{
	  type_ = type; 
	}
	
	IonDetector::AcquisitionMode IonDetector::getAcquisitionMode() const 
	{
	  return acquisition_mode_; 
	}
	
	void IonDetector::setAcquisitionMode(IonDetector::AcquisitionMode acquisition_mode)
	{
	  acquisition_mode_ = acquisition_mode; 
	}
	
	DoubleReal IonDetector::getResolution() const 
	{
	  return resolution_; 
	}
	
	void IonDetector::setResolution(DoubleReal resolution)
	{
	  resolution_ = resolution; 
	}
	
	DoubleReal IonDetector::getADCSamplingFrequency() const 
	{
	  return ADC_sampling_frequency_; 
	}
	
	void IonDetector::setADCSamplingFrequency(DoubleReal ADC_sampling_frequency)
	{
	  ADC_sampling_frequency_ = ADC_sampling_frequency; 
	}

  Int IonDetector::getOrder() const
  {
  	return order_;
  }
  
  void IonDetector::setOrder(Int order)
  {
  	order_ = order;
  }

}	
	
