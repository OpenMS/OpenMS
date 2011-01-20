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

#include <OpenMS/METADATA/Instrument.h>

using namespace std;

namespace OpenMS
{

const std::string Instrument::NamesOfIonOpticsType[] = {"Unknown","magnetic deflection","delayed extraction","collision quadrupole","selected ion flow tube","time lag focusing","reflectron","einzel lens","first stability region","fringing field","kinetic energy analyzer","static field"};

	Instrument::Instrument():
		MetaInfoInterface(),
		ion_optics_(UNKNOWN)
	{
	
	}
	
	Instrument::Instrument(const Instrument& source):
		MetaInfoInterface(source),
	  name_(source.name_),
	  vendor_(source.vendor_),
	  model_(source.model_),
	  customizations_(source.customizations_),
	  ion_sources_(source.ion_sources_),
	  mass_analyzers_(source.mass_analyzers_),
	  ion_detectors_(source.ion_detectors_),
	  software_(source.software_),
		ion_optics_(source.ion_optics_)
	{
	  
	}
	
	Instrument::~Instrument()
	{
	  
	}
	
	Instrument& Instrument::operator = (const Instrument& source)
	{
		if (&source == this) return *this;
	  
  	MetaInfoInterface::operator=(source);
	  software_ = source.software_;
    name_ = source.name_;
    vendor_ = source.vendor_;
    model_ = source.model_;
    customizations_ = source.customizations_;
    ion_sources_ = source.ion_sources_;
    mass_analyzers_ = source.mass_analyzers_;
    ion_detectors_ = source.ion_detectors_;
		ion_optics_ = source.ion_optics_;
	  
	  return *this;
	}
	
	bool Instrument::operator== (const Instrument& rhs) const
	{
//		if (software_ != rhs.software_) cout << "Instrument - " << __LINE__ << endl;
//		if (name_ != rhs.name_) cout << "Instrument - " << __LINE__ << endl;
//		if (vendor_ != rhs.vendor_) cout << "Instrument - " << __LINE__ << endl;
//		if (model_ != rhs.model_) cout << "Instrument - " << __LINE__ << endl;
//		if (customizations_ != rhs.customizations_) cout << "Instrument - " << __LINE__ << endl;
//		if (ion_sources_ != rhs.ion_sources_) cout << "Instrument - " << __LINE__ << endl;
//		if (mass_analyzers_ != rhs.mass_analyzers_) cout << "Instrument - " << __LINE__ << endl;
//		if (ion_detectors_ != rhs.ion_detectors_) cout << "Instrument - " << __LINE__ << endl;
//		if (ion_optics_ != rhs.ion_optics_) cout << "Instrument - " << __LINE__ << endl;
//		if (MetaInfoInterface::operator!=(rhs)) cout << "Instrument - " << __LINE__ << endl;
	
		return
	 		software_ == rhs.software_ &&
	    name_ == rhs.name_ &&
	    vendor_ == rhs.vendor_ &&
	    model_ == rhs.model_ &&
	    customizations_ == rhs.customizations_ &&
	    ion_sources_ == rhs.ion_sources_ &&
	    mass_analyzers_ == rhs.mass_analyzers_ &&
	    ion_detectors_ == rhs.ion_detectors_ &&
			ion_optics_ == rhs.ion_optics_ &&
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
	
	const std::vector<IonSource>& Instrument::getIonSources() const 
	{
	  return ion_sources_; 
	}
	
	std::vector<IonSource>& Instrument::getIonSources()
	{
	  return ion_sources_; 
	}
	
	void Instrument::setIonSources(const std::vector<IonSource>& ion_sources)
	{
	  ion_sources_ = ion_sources; 
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
	
	const std::vector<IonDetector>& Instrument::getIonDetectors() const 
	{
	  return ion_detectors_; 
	}
	
	std::vector<IonDetector>&  Instrument::getIonDetectors()
	{
	  return ion_detectors_; 
	}
	
	void Instrument::setIonDetectors(const std::vector<IonDetector>& ion_detectors)
	{
	  ion_detectors_ = ion_detectors; 
	}
	
	const Software& Instrument::getSoftware() const
	{
	  return software_;
	}
	
	Software& Instrument::getSoftware()
	{
	  return software_;
	}
	
	void Instrument::setSoftware(const Software& software)
	{
	  software_ = software;
	}

  Instrument::IonOpticsType Instrument::getIonOptics() const
  {
  	return ion_optics_;
  }
  
  void Instrument::setIonOptics(Instrument::IonOpticsType ion_optics)
  {
  	ion_optics_ = ion_optics;
  }

}


