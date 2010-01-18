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
// $Maintainer: Nico Pfeifer, Chris Bielow $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PROTEINIDENTIFICATION_H
#define OPENMS_METADATA_PROTEINIDENTIFICATION_H

#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <set>

namespace OpenMS
{   	
  /**
    @brief Representation of a peptide/protein ProteinIdentification
    
    This class stores the general information and the protein hits of an ProteinIdentification run.
    
    The actual peptide hits are stored in PeptideIdentification instances that are part of 
    spectra or features. 
    
    In order to be able to connect the ProteinIdentification and the corresponding peptide identifications, both
    classes have a string identifier. We recommend using the search engine name and the date as identifier.
    Setting this identifier is especially important, when there are several ProteinIdentification runs for a map
    i.e. several ProteinIdentification instances.

    @todo Add MetaInfoInterface to modifications => update IdXMLFile and ProteinIdentificationVisualizer (Andreas)
     
		@ingroup Metadata
  */
  class OPENMS_DLLAPI ProteinIdentification
  	: public MetaInfoInterface
  {
	  public:
	 		/// Hit type definition
	 		typedef ProteinHit HitType;
			
			/**
				@brief Connects multiple proteins which are indistinguishable
			*/
			struct ProteinGroup
			{
				bool operator == (const ProteinGroup rhs) const
				{
					return (id == rhs.id &&
									probability == rhs.probability &&
									indices == rhs.indices);
				}
			
				/// id of the group
				String id;
				/// probability of this group
				DoubleReal probability;
				/// Indices to (indistinguishable) proteins belonging to the same group
				std::set<Size> indices;
			};
			
			/// Peak mass type
			enum PeakMassType
			{
				MONOISOTOPIC,
				AVERAGE,
				SIZE_OF_PEAKMASSTYPE
			};
			/// Names corresponding to peak mass types
			static const std::string NamesOfPeakMassType[SIZE_OF_PEAKMASSTYPE];
			
			
			enum DigestionEnzyme
			{
				TRYPSIN,
				PEPSIN_A,
				PROTEASE_K,
				CHYMOTRYPSIN,
				NO_ENZYME,
				UNKNOWN_ENZYME,
				SIZE_OF_DIGESTIONENZYME
			};
			/// Names corresponding to digestion enzymes
			static const std::string NamesOfDigestionEnzyme[SIZE_OF_DIGESTIONENZYME];
			
			/// Search parameters of the DB search
			struct SearchParameters
				: public MetaInfoInterface
			{
				String db; ///< The used database
				String db_version; ///< The database version
				String taxonomy; ///< The taxonomy restriction
				String charges; ///< The allowed charges for the search
				PeakMassType mass_type; ///< Mass type of the peaks
				std::vector<String> fixed_modifications; ///< Used fixed modifications
				std::vector<String> variable_modifications; ///< Allowed variable modifications
				DigestionEnzyme enzyme; ///< The enzyme used for cleavage
				UInt missed_cleavages; ///< The number of allowed missed cleavages
				DoubleReal peak_mass_tolerance; ///< Mass tolerance of fragment ions (Dalton)
				DoubleReal precursor_tolerance; ///< Mass tolerance of precursor ions (Dalton)
				
				SearchParameters()
					: db(),
						db_version(),
						taxonomy(),
						charges(),
						mass_type(MONOISOTOPIC),
						fixed_modifications(),
						variable_modifications(),
						enzyme(UNKNOWN_ENZYME),
						missed_cleavages(0),
						peak_mass_tolerance(0.0),
						precursor_tolerance(0.0)
				{
				};

				bool operator == (const SearchParameters& rhs) const
				{
					return 	db == rhs.db &&
									db_version == rhs.db_version &&
									taxonomy == rhs.taxonomy &&
									charges == rhs.charges &&
									mass_type == rhs.mass_type &&
									fixed_modifications == rhs.fixed_modifications &&
									variable_modifications == rhs.variable_modifications &&
									enzyme == rhs.enzyme &&
									missed_cleavages == rhs.missed_cleavages &&
									peak_mass_tolerance == rhs.peak_mass_tolerance &&
									precursor_tolerance == rhs.precursor_tolerance;
				}

				bool operator != (const SearchParameters& rhs) const
				{
					return !(*this == rhs);
				}
			};
	  	
	  	
	    /** @name constructors,destructors,assignment operator <br> */
	    //@{
	    /// default constructor
	    ProteinIdentification();
	    /// destructor
	    virtual ~ProteinIdentification();
	    /// copy constructor
	    ProteinIdentification(const ProteinIdentification& source);
	    /// assignment operator
	    ProteinIdentification& operator=(const ProteinIdentification& source);
			/// Equality operator
			bool operator == (const ProteinIdentification& rhs) const;		
			/// Inequality operator
			bool operator != (const ProteinIdentification& rhs) const;
	    //@}	 

	   	///@name Protein hit information
	  	//@{	
	    /// returns the protein hits
	    const std::vector<ProteinHit>& getHits() const;
	    /// returns the protein hits (mutable)
	    std::vector<ProteinHit>& getHits();
			/// Appends a protein hit
	    void insertHit(const ProteinHit& input);
			/// Sets the peptide and protein hits
	    void setHits(const std::vector<ProteinHit>& hits); 

			/// returns the protein groups
			const std::vector<ProteinGroup>& getProteinGroups() const;
			/// returns the protein groups (mutable)
			std::vector<ProteinGroup>& getProteinGroups();
			/// appends a new protein group
			void insertGroup (const ProteinGroup& group);

			/// returns the peptide significance threshold value
	    DoubleReal getSignificanceThreshold() const;
			/// Sets the peptide significance threshold value
			void setSignificanceThreshold(DoubleReal value);
	    /// Returns the protein score type
	    const String& getScoreType() const;   
	    /// Sets the protein score type
	    void setScoreType(const String& type);
	    /// Returns true if a higher score represents a better score
	    bool isHigherScoreBetter() const;   
	    /// Sets the orientation of the score (higher is better?)
	    void setHigherScoreBetter(bool higher_is_better);
			/// Sorts the peptide and protein hits according to their score
			void sort();
			/// Sorts the peptide hits by score and assigns ranks (best score has rank of 1)
	    void assignRanks();
			//@}

	   	///@name General information
	  	//@{
			/// returns the date of the ProteinIdentification
	    const DateTime& getDateTime() const;
			/// sets the date of the ProteinIdentification
	    void setDateTime(const DateTime& date);
			/// sets the search engine type
			void setSearchEngine(const String& search_engine);
			/// returns the type of search engine used
			const String& getSearchEngine() const;
			/// sets the search engine version
			void setSearchEngineVersion(const String& search_engine_version);
			/// returns the search engine version
			const String& getSearchEngineVersion() const;
			/// sets the search parameters
			void setSearchParameters(const SearchParameters& search_parameters);
			/// returns the search parameters 
			const SearchParameters& getSearchParameters() const; 
	    /// returns the identifier
	    const String& getIdentifier() const;
	    /// sets the identifier
	    void setIdentifier(const String& id);
			//@}
			
	  protected:
	  
	  	template<class T, class PRED> 
			struct ProteinCmp 
			{
				ProteinCmp(T arr, const PRED& func) : arr_(arr), func_(func) 
				{}

				bool operator()(const size_t a, const size_t b)
				{
					return func_(arr_[a], arr_[b]);
				}
				T arr_;
				PRED func_;
			};

			///@name General information (search engine, parameters and DB)
	  	//@{
			String id_;
			String search_engine_;
			String search_engine_version_;
			SearchParameters search_parameters_;
	    DateTime date_;
	    //@}
	   	
			///@name Protein hit information
	  	//@{
	    String protein_score_type_;   
			bool higher_score_better_;
		  std::vector<ProteinHit> protein_hits_; 
			std::vector<ProteinGroup> protein_groups_;
			DoubleReal protein_significance_threshold_;
	    //@}
  };

} //namespace OpenMS
#endif // OPENMS_METADATA_PROTEINIDENTIFICATION_H
