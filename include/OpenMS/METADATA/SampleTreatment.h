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

#ifndef OPENMS_METADATA_SAMPLETREATMENT_H
#define OPENMS_METADATA_SAMPLETREATMENT_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/** 
		@brief Base class for sample treatments (Digestion, Modification, Tagging, ...)
		
		Virtual base class for all sample treatments.
		
		The type of the treatment can be determined through the getType() method.
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI SampleTreatment: public MetaInfoInterface
  {
    public:
    	/**
    		@brief Constructor.
    		
    		Use a unique type string for each treatment type
    	*/
    	SampleTreatment(const String& type);
    	/**
    		@brief Copy constructor
    		
    		@note Do not forget to call it when you derive a class from SampleTreatment!
    	*/
    	SampleTreatment(const SampleTreatment&);
    	/// destructor
    	virtual ~SampleTreatment();
    	
    	/**
    		@brief Assignment operator
    		
    		@note Do not forget to call it when you derive a class from SampleTreatment!
    	*/
    	SampleTreatment& operator=(const SampleTreatment&);

    	/**
    		@brief Equality operator
    		
    		The equality operators of derived classes also take a SampleTreatment reference as argument.
    		They check the type and cast the reference to the right type if the type matches.
    		
      	@note Do not forget to call it when you derive a class from SampleTreatment!
      */
      virtual bool operator== (const SampleTreatment& rhs) const;
    	
    	/**
    		@brief return the treatment type
    		
    		The type_ has to be set in the default constructor. 
    		It is used to determine the kind of sample treatment, when only a pointer to this base class is available.
			*/
			const String& getType() const;
			
    	/// returns the description of the sample treatment
			const String& getComment() const;
			
      /// sets the description of the sample treatment
      void setComment(const String& comment);
			
			/**
				@brief A clone methode
				
				clone method that creates a copy and retuns a pointer (base class pointer). 
				Used to copy sample treatments when only a pointer to this base class is available. 
			*/
			virtual SampleTreatment* clone() const =0;

    protected:
    	String type_;
    	String comment_;
    
    private:
    	/// Default constructor hidden to force setting of a type
    	SampleTreatment();

  };
} // namespace OpenMS

#endif // OPENMS_METADATA_SAMPLETREATMENT_H



