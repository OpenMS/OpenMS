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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
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
  class OPENMS_DLLAPI Digestion: public SampleTreatment
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

    	/**
    		@brief Equality operator
      	
      	Although this operator takes a reference to a SampleTreatment as argument
      	it tests for the equality of Tagging instances!
      */
      virtual bool operator== (const SampleTreatment& rhs) const;
			
			/// clone method. See SampleTreatment
			virtual SampleTreatment* clone() const;
			
			/// returns the enzyme name (default is "")
		  const String& getEnzyme() const;
		  /// sets the enzyme name
		  void setEnzyme(const String& enzyme);
		
			/// returns the digestion time in minutes (default is 0.0)
		  DoubleReal getDigestionTime() const;
		  /// sets the digestion time in minutes
		  void setDigestionTime(DoubleReal digestion_time);
		
			/// return the temperature during digestion in degree C (default is 0.0)
		  DoubleReal getTemperature() const;
		  /// sets the temperature during digestion in degree C
		  void setTemperature(DoubleReal temperature);
		
			/// returns the pH value (default is 0.0)
		  DoubleReal getPh() const;
		  /// sets the pH value
		  void setPh(DoubleReal ph);			

    protected:
			String enzyme_;
			DoubleReal digestion_time_;
			DoubleReal temperature_;
			DoubleReal ph_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_DIGESTION_H



