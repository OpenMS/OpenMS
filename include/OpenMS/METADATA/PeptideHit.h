// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PEPTIDEHIT_H
#define OPENMS_METADATA_PEPTIDEHIT_H

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

namespace OpenMS
{
  /**
    @brief Representation of a peptide hit
    
    It contains the fields score, score_type, rank, and sequence.
		
		@ingroup Metadata
  */
  class PeptideHit
  {
  	public:
  
		/**	@name Constructors and Destructor */
		//@{
		
		/// default constructor
    PeptideHit();
    
		/// values constructor
    PeptideHit(double score, 
    					 std::string score_type, 
    					 uint rank, 
							 SignedInt charge,
    					 String sequence);

		/// copy constructor
    PeptideHit(const PeptideHit& source);
				
		/// destructor
    virtual ~PeptideHit();
    //@}
    
		/// assignment operator
    PeptideHit& operator=(const PeptideHit& source);
		/// Equality operator
		bool operator == (const PeptideHit& rhs) const;
		/// Inequality operator
		bool operator != (const PeptideHit& rhs) const;

		/**	@name Accessors */
		//@{
		
    /// returns the score of the peptide hit 
    float getScore() const;
    
		/// returns the type of the score
    const std::string& getScoreType() const;
		
		/// returns the rank of the peptide hit
    UnsignedInt getRank() const;
		
		/// returns the peptide sequence without trailing or following spaces
  	String getSequence() const;
		
		/// returns the carge of the peptide
		SignedInt getCharge() const;
		
		/// returns the corresponding protein indices
		const std::vector< std::pair<String, String> >& getProteinIndices() const;
		
		/// returns a mutable reference to the corresponding protein indices
		std::vector< std::pair<String, String> >& getProteinIndices();

		/**
			@brief Sets the references to the protein hits of this peptide hit
		 	
			The format of one reference is < DateTime, Accession >
		*/
	  void setProteinIndices(const std::vector< std::pair<String, String> >& indices);
    
		/// sets the score of the peptide hit 
    void setScore(const double& score);
    
		/// sets the type of the score
    void setScoreType(const std::string& score_type);

		/// sets the rank
    void setRank(UnsignedInt newrank);
    
		/// sets the peptide sequence
		void setSequence(const String& sequence);

		/// sets the charge of the peptide
		void setCharge(SignedInt charge);
		
		/**
			@brief Adds a references to a protein hit of this peptide hit
		 	
			The format of one reference is < DateTime, Accession >
		*/
		void addProteinIndex(const std::pair<String, String>& index); 
				
		/**
			@brief Adds a references to a protein hit of this peptide hit
			 	
			The format of one reference is < DateTime, Accession >
		*/
		void addProteinIndex(const DateTime& date, const String& accession); 

    //@}

		/// clears all information of the peptide hit
    void clear();
  protected:
    float score_;									///< the score of the peptide hit
    std::string score_type_;    	///< the score type of the peptide hit 
		UnsignedInt rank_;    				///< the position(rank) where the hit 
																	///< appeared in the hit list

		/// the charge of the peptide
		SignedInt charge_;
		
    String sequence_;							///< the amino acid sequence of the 
    															///< peptide hit
    std::vector< std::pair<String, String> > corresponding_protein_indices_; ///< the indices 
    															///< of the corresponding proteins

  };

} // namespace OpenMS

#endif // OPENMS_METADATA_PEPTIDEHIT_H
