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

#ifndef OPENMS_METADATA_MODIFICATION_H
#define OPENMS_METADATA_MODIFICATION_H

#include <OpenMS/METADATA/SampleTreatment.h>

#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS 
{
	/**
		@brief Meta information about chemical modification of a sample
		
		Representation of some kind of modification. 
		It hold information about what amino acids are modified and how much the mass changes.
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI Modification: public SampleTreatment
  {
    public:
    	/// specificity of the reagent.
    	enum SpecificityType {
    		AA ///< specified amino acids are modified
    		,AA_AT_CTERM ///< specified amino acids are modified, if they are at the C-terminus
    		,AA_AT_NTERM ///< specified amino acids are modified, if they are at the N-terminus
    		,CTERM ///< the C-termiuns is modified
    		,NTERM ///< the N-termiuns is modified
    		,SIZE_OF_SPECIFICITYTYPE
    		};
			/// Names of specifitiy types
			static const std::string NamesOfSpecificityType[SIZE_OF_SPECIFICITYTYPE];
    	
    	/// default constructor
    	Modification();
    	/// copy constructor
    	Modification(const Modification&);
    	/// destructor
    	virtual ~Modification();
			
			/// assignment operator
    	Modification& operator=(const Modification&);
    	
    	/**
    		@brief Equality operator
      	
      	Although this operator takes a reference to a SampleTreatment as argument
      	it tests for the equality of Tagging instances!
      */
      virtual bool operator== (const SampleTreatment& rhs) const;
				
			/// clone method. See SampleTreatment
			virtual SampleTreatment* clone() const;
			
			/// returns the name of the reagent that was used (default: "")
			const String& getReagentName() const;
		  /// sets the name of the reagent that was used
		  void setReagentName(const String& reagent_name);
			
			/// returns the mass change (default: 0.0)
		  DoubleReal getMass() const;
		  /// sets the mass change
		  void setMass(DoubleReal mass);
		
			/// returns the specificity of the the reagent (default: AA)
		  const SpecificityType& getSpecificityType() const;
		  /// sets the specificity of the the reagent
		  void setSpecificityType(const SpecificityType& specificity_type);
			
			/// returns a string containing the one letter code of the amino acids that are affected by the reagent. (default: "")
		  const String& getAffectedAminoAcids() const;
		  /// returns a string containing the one letter code of the amino acids that are affected by the reagent. Do not separate them by space, tab or comma!
		  void setAffectedAminoAcids(const String& affected_amino_acids);

    protected:
			String reagent_name_;
			DoubleReal mass_;
			SpecificityType specificity_type_;
			String affected_amino_acids_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_MODIFICATION_H



