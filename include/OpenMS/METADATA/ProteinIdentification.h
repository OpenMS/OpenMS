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

#ifndef OPENMS_METADATA_PROTEINIDENTIFICATION_H
#define OPENMS_METADATA_PROTEINIDENTIFICATION_H

#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/METADATA/ProteinHit.h>

#include <string>

namespace OpenMS
{
  /**
    @brief Representation of a protein identification derived from serveral Identification instances
    
	  Represents the protein hits with additional parameters.
		
		@ingroup Metadata
  */
  class ProteinIdentification
  {
  public:

    /** @name constructors,destructors,assignment operator <br> */
    //@{
    
    /// default constructor
    ProteinIdentification();
    
    /// destructor
    ~ProteinIdentification();
    
    /// copy constructor
    ProteinIdentification(const ProteinIdentification& source);
    
    /// assignment operator
    ProteinIdentification& operator=(const ProteinIdentification& source);
    //@}
    
		/// Equality operator
		bool operator == (const ProteinIdentification& rhs) const;
		
		/// Inequality operator
		bool operator != (const ProteinIdentification& rhs) const;

		/// Sets the peptide and protein hits
    void setProteinHits(const std::vector<ProteinHit>& protein_hits);

    /// read access to protein hits
    const std::vector<ProteinHit>& getProteinHits() const;

		/// Inserts a protein hit into a container
    void insertProteinHit(const ProteinHit& input);

		/// retrival of the protein significance threshold value
    float getProteinSignificanceThreshold() const;

		/// setting of the protein significance threshold value
		void setProteinSignificanceThreshold(float value);

		/// assigns ranks to the stored protein hits according to their score
    void assignRanks();

		/// retrival of the date of the identification
    DateTime& getDateTime();

		/// retrival of the date of the identification
    const DateTime& getDateTime() const;

		/// retrival of the date of the identification
    void setDateTime(const DateTime& date);

		/// clears all information of this instance
    void clear();

		/// sorts the protein hits according to their score
		void sort();

		/// tests wether there is no information stored
		bool empty() const;

  protected:

    /// predicate for sorting the hits     
    class RankLess
    {
    public:
    	template<typename Arg>
      inline bool operator()(const Arg& a, const Arg& b)
      {
        return a.getRank() < b.getRank();
      }
    };
    
    /// predicate for sorting the hits     
    class ScoreMore
    {
    public:
    	template<typename Arg>
      bool operator()(const Arg& a, const Arg& b)
      {
        return a.getScore() > b.getScore();
      }
    };
   
    DateTime date_;																		///< the date when the search was performed
    std::vector<ProteinHit> protein_hits_;						///< a list containing the protein hits
    float protein_significance_threshold_;						///< the protein significance threshold
  };
} //namespace OpenMS
#endif //OPENMS_METADATA_PROTEINIDENTIFICATION_H
