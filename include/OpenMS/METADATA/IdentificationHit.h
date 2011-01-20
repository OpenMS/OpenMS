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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_IDENTIFICATIONHIT_H
#define OPENMS_METADATA_IDENTIFICATIONHIT_H

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{   	
  /**
    @brief Represents a object which can store the information of an analysisXML instance
    
	 	//@todo docu (Andreas)

		@ingroup Metadata
  */
  class OPENMS_DLLAPI IdentificationHit
  	: public MetaInfoInterface
  {
	  public:
			
	    /// @name constructors,destructors,assignment operator
	    //@{
	    /// default constructor
	    IdentificationHit();
	    /// destructor
	    virtual ~IdentificationHit();
	    /// copy constructor
	    IdentificationHit(const IdentificationHit& source);
	    /// assignment operator
	    IdentificationHit& operator=(const IdentificationHit& source);    
			/// Equality operator
			bool operator == (const IdentificationHit& rhs) const;
			/// Inequality operator
			bool operator != (const IdentificationHit& rhs) const;
	    //@}

			/// @name Accessors
			//@{
			/// sets the identifier
			void setId(const String& id);

			/// returns the id
			const String& getId() const;

			/// sets the charge state of the peptide
			void setCharge(Int charge);

			/// returns the charge state
			Int getCharge() const;

			/// sets the calculated mass to charge ratio
			void setCalculatedMassToCharge(DoubleReal mz);

			/// returns the calculated mass to charge ratio
			DoubleReal getCalculatedMassToCharge() const;

			/// sets the experimental mass to charge ratio
			void setExperimentalMassToCharge(DoubleReal mz);

			/// returns the experimental mass to charge
			DoubleReal getExperimentalMassToCharge() const;

			/// sets the name 
			void setName(const String& name);

			/// returns the name
			const String& getName() const;

			/// sets whether the peptide passed the threshold
			void setPassThreshold(bool pass);

			/// returns whether the peptide passed the threshold
			bool getPassThreshold() const;

			/// set the rank of the peptide
			void setRank(Int rank);

			/// returns the rank of the peptide
			Int getRank() const;
			//@}


	  protected:
			
			String id_;								///< identifier
			Int charge_; 							///< peptide charge
			DoubleReal calculated_mass_to_charge_; ///< calculated mass to charge ratio
			DoubleReal experimental_mass_to_charge_; ///< experimental mass to charge ratio
			String name_; 						///< name
			bool pass_threshold_; 		///< pass threshold
			Int rank_; 								///< rank of the peptide
  };

} //namespace OpenMS
#endif // OPENMS_METADATA_IDENTIFICATIONHIT_H
