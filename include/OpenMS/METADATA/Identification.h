// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_METADATA_IDENTIFICATION_H
#define OPENMS_METADATA_IDENTIFICATION_H

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/METADATA/SpectrumIdentification.h>

#include <vector>

namespace OpenMS
{   	
  /**
    @brief Represents a object which can store the information of an analysisXML instance
    
	 	//@todo docu (Andreas)

		@ingroup Metadata
  */
  class OPENMS_DLLAPI Identification
  	: public MetaInfoInterface
  {
	  public:
			
	    /// @name constructors,destructors,assignment operator
	    //@{
	    /// default constructor
	    Identification();
	    /// destructor
	    virtual ~Identification();
	    /// copy constructor
	    Identification(const Identification& source);
	    /// assignment operator
	    Identification& operator=(const Identification& source);    
			/// Equality operator
			bool operator == (const Identification& rhs) const;
			/// Inequality operator
			bool operator != (const Identification& rhs) const;
	    //@}

			/// @name Accessors
			//@{
			/// sets the date and time the file was written
			void setCreationDate(const DateTime& date);

			/// returns the date and time the file was created
			const DateTime& getCreationDate() const;

			/// sets the spectrum identifications
			void setSpectrumIdentifications(const std::vector<SpectrumIdentification>& ids);

			/// adds a spectrum identification
			void addSpectrumIdentification(const SpectrumIdentification& id);

			/// returns the spectrum identifications stored
			const std::vector<SpectrumIdentification>& getSpectrumIdentifications() const;
			//@}
	  protected:
			
			String id_;								///< Identifier
			DateTime creation_date_;	///< Date and time the search was performed
			std::vector<SpectrumIdentification> spectrum_identifications_;

  };

} //namespace OpenMS
#endif // OPENMS_METADATA_IDENTIFICATION_H
