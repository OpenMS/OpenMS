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

#ifndef OPENMS_METADATA_PEPTIDEHIT_H
#define OPENMS_METADATA_PEPTIDEHIT_H

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

namespace OpenMS
{
  /**
    @brief Representation of a peptide hit
    
    It contains the fields score, score_type, rank, and sequence.
		
		@ingroup Metadata
  */
  class OPENMS_DLLAPI PeptideHit
  	: public MetaInfoInterface
  {
  	public:

			/// @name Comparators for PeptideHit and ProteinHit
			//@{
			/// Greater predicate for scores of hits
			class OPENMS_DLLAPI ScoreMore
			{
			  public:
			  	template<typename Arg>
			    bool operator()(const Arg& a, const Arg& b)
			    {
			      return a.getScore() > b.getScore();
			    }
			};
			
			/// Lesser predicate for scores of hits
			class OPENMS_DLLAPI ScoreLess
			{
			  public:
			  	template<typename Arg>
			    bool operator()(const Arg& a, const Arg& b)
			    {
			      return a.getScore() < b.getScore();
			    }
			};
			//@}

			/**	@name Constructors and Destructor */
			//@{
			/// default constructor
	    PeptideHit();
	    
			/// values constructor
	    PeptideHit(DoubleReal score, 
	    					 UInt rank, 
								 Int charge,
	    					 const AASequence& sequence);
	
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
	
			/**	@name Accessors 
			*/
			//@{
	    /// returns the score of the peptide hit 
	    DoubleReal getScore() const;
	    
			/// returns the rank of the peptide hit
	    UInt getRank() const;
			
			/// returns the peptide sequence without trailing or following spaces
	  	const AASequence& getSequence() const;
			
			/// returns the carge of the peptide
			Int getCharge() const;
			
			/// returns the corresponding protein accessions
			const std::vector<String>& getProteinAccessions() const;
	
			/// sets the corresponding protein accessions
		  void setProteinAccessions(const std::vector<String>& accessions);
	    
			/// sets the score of the peptide hit 
	    void setScore(DoubleReal score);
	    
			/// sets the rank
	    void setRank(UInt newrank);
	    
			/// sets the peptide sequence
			void setSequence(const AASequence& sequence);
	
			/// sets the charge of the peptide
			void setCharge(Int charge);
			
			/// adds an accession of a protein which contains this peptide hit
			void addProteinAccession(const String& accession); 
			
			/// sets the amino acid before the sequence
			void setAABefore(char acid);
			/// returns the amino acid before the sequence
			char getAABefore() const;

			/// sets the amino acid after the sequence
			void setAAAfter(char acid);
			/// returns the amino acid after the sequence
			char getAAAfter() const;

			
	    //@}
	
	  
		protected:
	    DoubleReal score_;		///< the score of the peptide hit
			UInt rank_;    				///< the position(rank) where the hit appeared in the hit list
			Int charge_; ///< the charge of the peptide
	    AASequence sequence_;							///< the amino acid sequence of the peptide hit
	    char aa_before_; ///< Amino acid before the sequence
	    char aa_after_; ///< Amino acid after the sequence
	    std::vector<String> corresponding_protein_accessions_; ///< the accessions of the corresponding proteins
	
	};

} // namespace OpenMS

#endif // OPENMS_METADATA_PEPTIDEHIT_H
