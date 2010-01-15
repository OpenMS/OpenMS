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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PROTEINHIT_H
#define OPENMS_METADATA_PROTEINHIT_H

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
  /**
    @brief Representation of a protein hit
    
    It contains the fields score, score_type, rank, accession, 
    accession_type and sequence.
		
		@ingroup Metadata
  */
  class OPENMS_DLLAPI ProteinHit
  	: public MetaInfoInterface
  {
  	public:
		
			/**	@name Constructors and Destructor */
			//@{
			
			/// default constructor
	    ProteinHit();
	    
			/// values constructor
	    ProteinHit(DoubleReal score, UInt rank, String accession, String sequence);
	
			/// copy constructor
	    ProteinHit(const ProteinHit& source);
					
			/// destructor
	    virtual ~ProteinHit();
	    //@}
	    
			/// assignment operator
	    ProteinHit& operator=(const ProteinHit& source);
			
			/// Equality operator
			bool operator == (const ProteinHit& rhs) const;
			
			/// Inequality operator
			bool operator != (const ProteinHit& rhs) const;
	
	
			/**	@name Accessors */
			//@{
			
	    /// returns the score of the protein hit 
	    Real getScore() const;
	    
			/// returns the rank of the protein hit
	    UInt getRank() const;
			
			/// returns the protein sequence
			const String& getSequence() const;
			
			/// returns the accession of the protein
			const String& getAccession() const;
	
	    /// sets the score of the protein hit 
	    void setScore(DoubleReal score);
			
			/// sets the rank
	    void setRank(UInt newrank);
			
			/// sets the protein sequence
			void setSequence(const String& sequence);
			
			/// sets the accession of the protein
			void setAccession(const String& accession);
	    //@}

	  protected:
	    Real score_;					///< the score of the protein hit
			UInt rank_;    				///< the position(rank) where the hit appeared in the hit list
	    String accession_;	 	///< the protein identifier
	    String sequence_;		 	///< the amino acid sequence of the protein hit

  };

} // namespace OpenMS

#endif // OPENMS_METADATA_PROTEINHIT_H
