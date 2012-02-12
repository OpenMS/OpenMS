// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/METADATA/Digestion.h>

using namespace std;

namespace OpenMS
{
	
	
	Digestion::Digestion() :
		SampleTreatment("Digestion"),
		enzyme_(""), 
		digestion_time_(0.0), 
		temperature_(0.0), 
		ph_(0.0)
	{

	}
	
	Digestion::Digestion(const Digestion& source):
		SampleTreatment(source),
		enzyme_(source.enzyme_), 
		digestion_time_(source.digestion_time_), 
		temperature_(source.temperature_), 
		ph_(source.ph_)
	{

	}
	
	Digestion::~Digestion()
	{

	}
	
	SampleTreatment* Digestion::clone() const
	{
		SampleTreatment* tmp = new Digestion(*this);
		return tmp;
	}
	
	Digestion& Digestion::operator = (const Digestion& source)
	{
	  if (&source == this) return *this;
	  
	  SampleTreatment::operator=(source);
		enzyme_ = source.enzyme_;
		digestion_time_ = source.digestion_time_;
		temperature_ = source.temperature_;
		ph_ = source.ph_;
	  
	  return *this;
	}

  bool Digestion::operator== (const SampleTreatment& rhs) const
  {
  	if (type_!=rhs.getType()) return false;
  	
  	const Digestion* tmp = dynamic_cast<const Digestion*>(&rhs);
  	return
  		SampleTreatment::operator==(*tmp) &&
	 		enzyme_ == tmp->enzyme_ &&
			digestion_time_ == tmp->digestion_time_ &&
			temperature_ == tmp->temperature_ &&
			ph_ == tmp->ph_
  		;
  }
	
  const String& Digestion::getEnzyme() const
  {
    return enzyme_;
  }

  void Digestion::setEnzyme(const String& enzyme)
  {
    enzyme_ = enzyme;
  }


  DoubleReal Digestion::getDigestionTime() const
  {
    return digestion_time_;
  }

  void Digestion::setDigestionTime(DoubleReal digestion_time)
  {
    digestion_time_ = digestion_time;
  }


  DoubleReal Digestion::getTemperature() const
  {
    return temperature_;
  }

  void Digestion::setTemperature(DoubleReal temperature)
  {
    temperature_ = temperature;
  }


  DoubleReal Digestion::getPh() const
  {
    return ph_;
  }

  void Digestion::setPh(DoubleReal ph)
  {
    ph_ = ph;
  }

}

