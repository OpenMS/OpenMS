// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: Digestion.h,v 1.1 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_DIGESTION_H
#define OPENMS_METADATA_DIGESTION_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/SampleTreatment.h>

namespace OpenMS 
{
	/**
		@brief Meta information about digestion of a sample
		
		Representation of a digestion.
		
		@ingroup Metadata
	*/
  class Digestion: public SampleTreatment
  {
    public:
    	/// default constructor
    	Digestion();
    	/// copy constructor
    	Digestion(const Digestion&);
    	/// destructor
    	virtual ~Digestion();

			/// assignment operator
    	Digestion& operator=(const Digestion&);

    	/// Equality operator
      virtual bool operator== (const SampleTreatment& rhs) const;
			
			/// clone method. See SampleTreatment
			virtual SampleTreatment* clone() const;
			
			/// returns the enzyme name
		  const String& getEnzyme() const;
		  /// sets the enzyme name
		  void setEnzyme(const String& enzyme);
		
			/// returns the digestion time in minutes
		  float getDigestionTime() const;
		  /// sets the digestion time in minutes
		  void setDigestionTime(float digestion_time);
		
			/// return the temperature during digestion in °C
		  float getTemperature() const;
		  /// sets the temperature during digestion in °C
		  void setTemperature(float temperature);
		
			/// returns the pH value
		  float getPh() const;
		  /// sets the pH value
		  void setPh(float ph);			

    protected:
			String enzyme_;
			float digestion_time_;
			float temperature_;
			float ph_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_DIGESTION_H



