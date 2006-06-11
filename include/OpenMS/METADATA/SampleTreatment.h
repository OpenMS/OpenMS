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
// $Id: SampleTreatment.h,v 1.1 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SAMPLETREATMENT_H
#define OPENMS_METADATA_SAMPLETREATMENT_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/** 
		@brief Base class for sample treatments (Digestion, Modification, Tagging, ...)
		
		Virtual base class for all sample treatments.
		
		@ingroup Metadata
	*/
  class SampleTreatment: public MetaInfoInterface
  {
    public:
    	/**
    		@brief Constructor.
    		
    		Use a unique type string for each treatment type
    	*/
    	SampleTreatment(const String& type);
    	/**
    		@brief Copy constructor
    		
    		Do not forget to call it when you derive a class from SampleTreatment!
    	*/
    	SampleTreatment(const SampleTreatment&);
    	/// destructor
    	virtual ~SampleTreatment();
    	
    	/**
    		@brief Assignment operator
    		
    		Do not forget to call it when you derive a class from SampleTreatment!
    	*/
    	SampleTreatment& operator=(const SampleTreatment&);

    	/// Equality operator
      virtual bool operator== (const SampleTreatment& rhs) const =0;
      /// Equality operator
      bool operator!= (const SampleTreatment& rhs) const;
    	
    	/**
    		@brief return the treatment type
    		
    		The type_ has to be set in the default constructor. 
    		It is used to determine the kind of sample treatment, when only a pointer to this base class is available.
			*/
			const String& getType() const;
			
			/**
				@brief A clone methode
				
				clone method that creates a copy and retuns a pointer (base class pointer). 
				Used to copy sample treatments when only a pointer to this base class is available. 
			*/
			virtual SampleTreatment* clone() const =0;

    protected:
    	String type_;
    
    private:
    	/// Default constructor hidden to force setting of a type
    	SampleTreatment();

  };
} // namespace OpenMS

#endif // OPENMS_METADATA_SAMPLETREATMENT_H



