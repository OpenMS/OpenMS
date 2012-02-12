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
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PROTEINIDENTIFICATION_H
#define OPENMS_METADATA_PROTEINIDENTIFICATION_H

#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <set>

namespace OpenMS
{
  class PeptideIdentification;
  /**
    @brief Representation of a protein identification run
    
    This class stores the general information and the protein hits of a protein identification run.
    
    The actual peptide hits are stored in PeptideIdentification instances that are part of spectra or features. 
    
    In order to be able to connect the ProteinIdentification and the corresponding peptide identifications, both classes have a string identifier. We recommend using the search engine name and the date as identifier.
    Setting this identifier is especially important when there are several protein identification runs for a map, i.e. several ProteinIdentification instances.

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
				@brief Bundles multiple (e.g. indistinguishable) proteins in a group
			*/
			struct ProteinGroup
			{
				/// Probability of this group
				DoubleReal probability;
				/// Accessions of (indistinguishable) proteins that belong to the same group
				StringList accessions;

				ProteinGroup(): probability(0.0), accessions()
				{}

				bool operator==(const ProteinGroup rhs) const
				{
					return (probability == rhs.probability &&
									accessions == rhs.accessions);
				}
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
        }

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
	  	
	  	
	    /** @name Constructors, destructors, assignment operator <br> */
	    //@{
	    /// Default constructor
	    ProteinIdentification();
	    /// Destructor
	    virtual ~ProteinIdentification();
	    /// Copy constructor
	    ProteinIdentification(const ProteinIdentification& source);
	    /// Assignment operator
	    ProteinIdentification& operator=(const ProteinIdentification& source);
			/// Equality operator
			bool operator == (const ProteinIdentification& rhs) const;		
			/// Inequality operator
			bool operator != (const ProteinIdentification& rhs) const;
	    //@}	 

	   	///@name Protein hit information (public members)
	  	//@{
	    /// Returns the protein hits
	    const std::vector<ProteinHit>& getHits() const;
	    /// Returns the protein hits (mutable)
	    std::vector<ProteinHit>& getHits();
			/// Appends a protein hit
	    void insertHit(const ProteinHit& input);
			/// Sets the protein hits
	    void setHits(const std::vector<ProteinHit>& hits); 
			/// Finds a protein hit by accession (returns past-the-end iterator if not found)
			std::vector<ProteinHit>::iterator findHit(const String& accession);

			/// Returns the protein groups
			const std::vector<ProteinGroup>& getProteinGroups() const;
			/// Returns the protein groups (mutable)
			std::vector<ProteinGroup>& getProteinGroups();
			/// Appends a new protein group
			void insertProteinGroup (const ProteinGroup& group);

			/// Returns the indistinguishable proteins
			const std::vector<ProteinGroup>& getIndistinguishableProteins() const;
			/// Returns the indistinguishable proteins (mutable)
			std::vector<ProteinGroup>& getIndistinguishableProteins();
			/// Appends new indistinguishable proteins
			void insertIndistinguishableProteins(const ProteinGroup& group);

			/// Returns the protein significance threshold value
	    DoubleReal getSignificanceThreshold() const;
			/// Sets the protein significance threshold value
			void setSignificanceThreshold(DoubleReal value);
	    /// Returns the protein score type
	    const String& getScoreType() const;   
	    /// Sets the protein score type
	    void setScoreType(const String& type);
	    /// Returns true if a higher score represents a better score
	    bool isHigherScoreBetter() const;   
	    /// Sets the orientation of the score (is higher better?)
	    void setHigherScoreBetter(bool higher_is_better);
			/// Sorts the protein hits according to their score
			void sort();
			/// Sorts the protein hits by score and assigns ranks (best score has rank 1)
	    void assignRanks();
      /**
				 @brief Compute the coverage (in percent) of all ProteinHits given PeptideHits
      
				 @throws Exception::MissingInformation if ProteinsHits do not have sequence information

				 @return The number of Proteins referenced by the @p pep_ids that are not contained in this ProteinIdentification set (should be 0)
			*/
      Size computeCoverage(const std::vector<PeptideIdentification>& pep_ids);
			//@}

	   	///@name General information
	  	//@{
			/// Returns the date of the protein identification run
	    const DateTime& getDateTime() const;
			/// Sets the date of the protein identification run
	    void setDateTime(const DateTime& date);
			/// Sets the search engine type
			void setSearchEngine(const String& search_engine);
			/// Returns the type of search engine used
			const String& getSearchEngine() const;
			/// Sets the search engine version
			void setSearchEngineVersion(const String& search_engine_version);
			/// Returns the search engine version
			const String& getSearchEngineVersion() const;
			/// Sets the search parameters
			void setSearchParameters(const SearchParameters& search_parameters);
			/// Returns the search parameters 
			const SearchParameters& getSearchParameters() const; 
	    /// Returns the identifier
	    const String& getIdentifier() const;
	    /// Sets the identifier
	    void setIdentifier(const String& id);
			//@}
			
	  protected:
			///@name General information (search engine, parameters and database)
	  	//@{
			String id_;
			String search_engine_;
			String search_engine_version_;
			SearchParameters search_parameters_;
	    DateTime date_;
	    //@}
	   	
			///@name Protein hit information (protected members)
	  	//@{
	    String protein_score_type_;   
			bool higher_score_better_;
		  std::vector<ProteinHit> protein_hits_; 
			std::vector<ProteinGroup> protein_groups_;
			/// Indistinguishable proteins: @p accessions[0] is "group leader", @p probability is meaningless
			std::vector<ProteinGroup> indistinguishable_proteins_;
			DoubleReal protein_significance_threshold_;
	    //@}
  };

} //namespace OpenMS
#endif // OPENMS_METADATA_PROTEINIDENTIFICATION_H
