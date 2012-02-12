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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SPECTRUMIDENTIFICATION_H
#define OPENMS_METADATA_SPECTRUMIDENTIFICATION_H

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/IdentificationHit.h>

namespace OpenMS
{   	
  /**
    @brief Represents a object which can store the information of an analysisXML instance
    
	 	//@todo docu (Andreas)

		@ingroup Metadata
  */
  class OPENMS_DLLAPI SpectrumIdentification
  	: public MetaInfoInterface
  {
	  public:
			
	    /// @name constructors,destructors,assignment operator
	    //@{
	    /// default constructor
	    SpectrumIdentification();
	    /// destructor
	    virtual ~SpectrumIdentification();
	    /// copy constructor
	    SpectrumIdentification(const SpectrumIdentification& source);
	    /// assignment operator
	    SpectrumIdentification& operator=(const SpectrumIdentification& source);    
			/// Equality operator
			bool operator == (const SpectrumIdentification& rhs) const;
			/// Inequality operator
			bool operator != (const SpectrumIdentification& rhs) const;
	    //@}

			// @name Accessors
			//@{
			/// sets the identification hits of this spectrum identification (corresponds to single peptide hit in the list)
			void setHits(const std::vector<IdentificationHit>& hits);

			/// adds a single identification hit to the hits
			void addHit(const IdentificationHit& hit);

			/// returns the identificatio hits of this spectrum identification
			const std::vector<IdentificationHit>& getHits() const;
			//@}

	  protected:
			
			String id_;																///< Identifier
			std::vector<IdentificationHit> hits_;  		///< Single peptide hits
  };

} //namespace OpenMS
#endif // OPENMS_METADATA_SPECTRUMIDENTIFICATION_H
