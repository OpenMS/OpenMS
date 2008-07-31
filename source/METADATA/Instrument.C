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

#include <OpenMS/METADATA/Instrument.h>

using namespace std;

namespace OpenMS
{

	Instrument::Instrument():
		MetaInfoInterface()
	{
	
	}
	
	Instrument::Instrument(const Instrument& source):
		MetaInfoInterface(source),
	  name_(source.name_),
	  vendor_(source.vendor_),
	  model_(source.model_),
	  customizations_(source.customizations_),
	  ion_source_(source.ion_source_),
	  mass_analyzers_(source.mass_analyzers_),
	  ion_detector_(source.ion_detector_)
	{
	  
	}
	
	Instrument::~Instrument()
	{
	  
	}
	
	Instrument& Instrument::operator = (const Instrument& source)
	{
		if (&source == this) return *this;
	  
    name_ = source.name_;
    vendor_ = source.vendor_;
    model_ = source.model_;
    customizations_ = source.customizations_;
    ion_source_ = source.ion_source_;
    mass_analyzers_ = source.mass_analyzers_;
    ion_detector_ = source.ion_detector_;
  	MetaInfoInterface::operator=(source);
	  
	  return *this;
	}
	
	bool Instrument::operator== (const Instrument& rhs) const
	{
		return 
	    name_ == rhs.name_ &&
	    vendor_ == rhs.vendor_ &&
	    model_ == rhs.model_ &&
	    customizations_ == rhs.customizations_ &&
	    ion_source_ == rhs.ion_source_ &&
	    mass_analyzers_ == rhs.mass_analyzers_ &&
	    ion_detector_ == rhs.ion_detector_ &&
			MetaInfoInterface::operator==(rhs)
			;
	}
	
	bool Instrument::operator!= (const Instrument& rhs) const
	{
		return !(operator==(rhs));
	}
	
	const String& Instrument::getName() const 
	{
	  return name_; 
	}
	
	void Instrument::setName(const String& name)
	{
	  name_ = name; 
	}
	
	const String& Instrument::getVendor() const 
	{
	  return vendor_; 
	}
	
	void Instrument::setVendor(const String& vendor)
	{
	  vendor_ = vendor; 
	}
	
	const String& Instrument::getModel() const 
	{
	  return model_; 
	}
	
	void Instrument::setModel(const String& model)
	{
	  model_ = model; 
	}
	
	const String& Instrument::getCustomizations() const 
	{
	  return customizations_; 
	}
	
	void Instrument::setCustomizations(const String& customizations)
	{
	  customizations_ = customizations; 
	}
	
	const IonSource& Instrument::getIonSource() const 
	{
	  return ion_source_; 
	}
	
	IonSource&  Instrument::getIonSource()
	{
	  return ion_source_; 
	}
	
	void Instrument::setIonSource(const IonSource& ion_source)
	{
	  ion_source_ = ion_source; 
	}
	
	const std::vector<MassAnalyzer>& Instrument::getMassAnalyzers() const 
	{
	  return mass_analyzers_; 
	}
	
	std::vector<MassAnalyzer>&  Instrument::getMassAnalyzers()
	{
	  return mass_analyzers_; 
	}
	
	void Instrument::setMassAnalyzers(const std::vector<MassAnalyzer>& mass_analyzers)
	{
	  mass_analyzers_ = mass_analyzers; 
	}
	
	const IonDetector& Instrument::getIonDetector() const 
	{
	  return ion_detector_; 
	}
	
	IonDetector&  Instrument::getIonDetector()
	{
	  return ion_detector_; 
	}
	
	void Instrument::setIonDetector(const IonDetector& ion_detector)
	{
	  ion_detector_ = ion_detector; 
	}

}


