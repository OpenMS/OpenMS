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

#ifndef OPENMS_METADATA_IDENTIFICATION_H
#define OPENMS_METADATA_IDENTIFICATION_H

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <string>
#include <map>

namespace OpenMS
{   	
  /**
    @brief Representation of a database search with one spectrum
    
	  Represents the peptide and protein hits with additional parameters.
		
		@ingroup Metadata
  */
  class Identification : public ProteinIdentification
  {
  public:

    /** @name constructors,destructors,assignment operator <br> */
    //@{
    
    /// default constructor
    Identification();
    
    /// destructor
    virtual ~Identification();
    
    /// copy constructor
    Identification(const Identification& source);
    
    /// assignment operator
    Identification& operator=(const Identification& source);
    //@}
    
		/// Equality operator
		bool operator == (const Identification& rhs) const;
		
		/// Inequality operator
		bool operator != (const Identification& rhs) const;

    /// read access to peptide hits
    const std::vector<PeptideHit>& getPeptideHits() const;

    /// mutable access to peptide hits
    std::vector<PeptideHit>& getPeptideHits();

		/// Inserts a peptide hit into a container
    void insertPeptideHit(const PeptideHit& input);

		/// Sets the peptide and protein hits
    void setPeptideAndProteinHits(const std::vector<PeptideHit>& peptide_hits, 
    															const std::vector<ProteinHit>& protein_hits);

		/// retrival of the peptide significance threshold value
    float getPeptideSignificanceThreshold() const;

		/// setting of the peptide significance threshold value
		void setPeptideSignificanceThreshold(float value);

		/// clears all information of this instance
    void clear();

		/// sorts the peptide and protein hits according to their score
		void sort();

		/// tests whether there is no information stored
		bool empty() const;

		/// sorts the peptide hits by sort() and assigns ranks according to the sorting
    void assignRanks();
    
		/**
			@brief returns all referencing hits
    	
			Returns a vector of peptide hits that reference a protein hit with @code date_time @endcode
			and @code accession @endcode.
    */
    std::vector<PeptideHit>* getReferencingHits(String date_time, String accession) const;

		/**
			@brief returns all non referencing peptide hits
    	
			Returns a vector of peptide hits that do not reference a protein hit in [protein_hits_begin, protein_hits_end)
			together with the date given by @code date_time @endcode
    */
    template <class iteratorT>
  	std::vector<PeptideHit>* getNonReferencingHits(iteratorT protein_hits_begin, iteratorT protein_hits_end, const String&	date_time) const
  	{
	  	std::vector<PeptideHit>* 	found_hits 	= new std::vector<PeptideHit>();
	  	bool 											referenced 	= false;
	  	String 										accession 	= "";
	  			
	  	// for every peptide hit
			for(UnsignedInt i = 0; i < peptide_hits_.size(); i++)
			{
				const std::vector< std::pair<String, String> >& references = 
						peptide_hits_[i].getProteinIndices();
				//for every reference
				for(UnsignedInt j = 0; j < references.size(); j++)
				{
					// for every protein hit
					for(iteratorT protein_it = protein_hits_begin;
							protein_it != protein_hits_end;
							protein_it++)
					{
						accession = protein_it->getAccession();
						if (references[j].first == date_time && references[j].second == accession)
						{
							referenced = true;
						}
					}
				}
				if (!referenced)
				{
					found_hits->push_back(peptide_hits_[i]);
				}
				referenced = false;
			}
	
	  	return found_hits;  			  		
  	}

		/**
			@brief returns all non referencing peptide hits
    	
			Returns a vector of peptide hits that do not reference a protein hit in @code protein_hits @endcode.
			The String argument in the map stands for the String representation of the date_time object of
			the ProteinIdentification the particular protein hit belongs to.
    */
  	std::vector<PeptideHit>* getNonReferencingHits(const std::multimap< String, ProteinHit >& protein_hits) const;
		/// Inserts a protein hit into a container
	  void insertProteinHit(const ProteinHit& input);

  protected:

    std::vector<PeptideHit> peptide_hits_;						///< a list containing the peptide hits
    float peptide_significance_threshold_;						///< the peptide significance threshold
  };

	///Struct that contains an Identification and its corresponding RT and m/z value
	struct IdentificationData
	{
		double rt;
		double mz;
		Identification id;
	
		bool operator==(const IdentificationData& rhs) const
		{
			return rt==rhs.rt && mz==rhs.mz &&  id==rhs.id;
		}
	};
} //namespace OpenMS
#endif // OPENMS_METADATA_IDENTIFICATION_H
