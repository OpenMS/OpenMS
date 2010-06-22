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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_ID_IDFILTER_H
#define OPENMS_FILTERING_ID_IDFILTER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief used to filter identifications by different criteria
    
    The identifications are filtered by significance thresholds and
    by sequences. The filtering by significance thresholds looks for the 
    best ProteinIdentification that fullfills the significance threshold criterium.
    score > significance-threshold * significance_fraction. 
    The filtering by sequences looks for the best ProteinIdentification that
    is contained in one of the protein sequences.
  */
  class OPENMS_DLLAPI IDFilter
  {
    public:

      /// Constructor
      IDFilter();
      
      /// Destructor
      virtual ~IDFilter();

      /// filters a ProteinIdentification or PeptideIdentification by only allowing peptides/proteins which reach a score above @p threshold_fraction * SignificanceThreshold
      template <class IdentificationType>
			void filterIdentificationsByThreshold(const IdentificationType& identification, DoubleReal threshold_fraction, IdentificationType& filtered_identification)
			{
				typedef typename IdentificationType::HitType HitType;
				std::vector<HitType> temp_hits;
				std::vector<HitType> filtered_hits;

				filtered_identification = identification;
				filtered_identification.setHits(std::vector<HitType>());
				
				for(typename std::vector<HitType>::const_iterator it = identification.getHits().begin();
					it != identification.getHits().end(); 
					++it)
				{
					if (it->getScore() >= threshold_fraction * identification.getSignificanceThreshold())
					{	
	 					filtered_hits.push_back(*it);
					}	
				}

				if (filtered_hits.size() > 0)
				{
		  		filtered_identification.setHits(filtered_hits);																						
					filtered_identification.assignRanks();  																								
				}
			}
			
		  /**
		    @brief filters a ProteinIdentification or PeptideIdentification corresponding to the @p threshold_score
		    
		    If the method higherScoreBetter() returns true for the IdentificationType all hits with a score
		    smaller than @p threshold_score are removed. Otherwise all hits with a score bigger than
		    @p threshold_score are removed.
		  */
      template <class IdentificationType>
			void filterIdentificationsByScore(const IdentificationType& identification, DoubleReal threshold_score, IdentificationType& filtered_identification)
			{
				typedef typename IdentificationType::HitType HitType;
				std::vector<HitType> temp_hits;
				std::vector<HitType> filtered_hits;

				filtered_identification = identification;
				filtered_identification.setHits(std::vector<HitType>());
				
				for(typename std::vector<HitType>::const_iterator it = identification.getHits().begin();
					it != identification.getHits().end(); 
					++it)
				{
					if (identification.isHigherScoreBetter())
					{
						if (it->getScore() >= threshold_score)
						{	
	 						filtered_hits.push_back(*it);
						}
					}
					else
					{
						if (it->getScore() <= threshold_score)
						{	
	 						filtered_hits.push_back(*it);
						}
					}	
				}

				if (filtered_hits.size() > 0)
				{
		  		filtered_identification.setHits(filtered_hits);																						
					filtered_identification.assignRanks();  																								
				}
			}
			
		  /**
		    @brief filters a ProteinIdentification or PeptideIdentification corresponding to the score.
		    
		    If the method higherScoreBetter() returns true for the IdentificationType the 
		    n highestscoring hits are kept. Otherwise the n lowest scoring hits are kept. 
		  */
      template <class IdentificationType>
			void filterIdentificationsByBestNHits(const IdentificationType& identification, Size n, IdentificationType& filtered_identification)
			{
				typedef typename IdentificationType::HitType HitType;
				std::vector<HitType> temp_hits;
				std::vector<HitType> filtered_hits;
				Size count = 0;
				
				IdentificationType temp_identification = identification;
				temp_identification.sort(); // .. by score

        filtered_identification = identification;
				filtered_identification.setHits(std::vector<HitType>());
				
				
				typename std::vector<HitType>::const_iterator it = temp_identification.getHits().begin();
				while(it != temp_identification.getHits().end()
							&& count < n)
				{
					filtered_hits.push_back(*it);
					++it;
					++count;	
				}

				if (filtered_hits.size() > 0)
				{
		  		filtered_identification.setHits(filtered_hits);																						
					filtered_identification.assignRanks();  																								
				}
			}
			
      /// filters a PeptideIdentification keeping only the best scoring hits (if strict is set, keeping only the best hit only if it is the only hit with that score)
			void filterIdentificationsByBestHits(const PeptideIdentification& identification, PeptideIdentification& filtered_identification, bool strict = false);

      /// filters a PeptideIdentification corresponding to the given proteins
      /// PeptideHits with no matching @em proteins are removed.
      /// Matching is done either based on accessions or on sequence (if no accessions are given, or @em no_protein_identifiers is set)
			void filterIdentificationsByProteins(const PeptideIdentification& identification, const std::vector< FASTAFile::FASTAEntry >& proteins, PeptideIdentification& filtered_identification, bool no_protein_identifiers = false);

      /// filters a ProteinIdentification corresponding to the given proteins
      /// ProteinHits with no matching @em proteins are removed.
      /// Matching is done based on accessions only
			void filterIdentificationsByProteins(const ProteinIdentification& identification, const std::vector< FASTAFile::FASTAEntry >& proteins, ProteinIdentification& filtered_identification);
																														
			/// removes all peptide hits having a sequence equal to a element in <code>peptides</code>
			void filterIdentificationsByExclusionPeptides(const PeptideIdentification& identification, const std::set<String>& peptides, PeptideIdentification& filtered_identification);
																														
		  /// only peptides having a length equal to or greater than 'length' will be kept
		  void filterIdentificationsByLength(const PeptideIdentification& identification, Size length, PeptideIdentification& filtered_identification);

			/// only protein hits in 'identification' which are referenced by a peptide in 'peptide_identifications' are kept
			void removeUnreferencedProteinHits(const ProteinIdentification& 	identification, const std::vector<PeptideIdentification> peptide_identifications, ProteinIdentification& 	filtered_identification);

			/// if a peptide hit occurs more than once, only one instance is kept
		  void filterIdentificationsUnique(const PeptideIdentification& identification, PeptideIdentification& filtered_identification);
		  /**
				@brief Filters the peptide hits according to their predicted rt p-values
				
				Filters the peptide hits of this ProteinIdentification by the 
				probability (p-value) of a correct ProteinIdentification having a deviation between 
				observed and predicted rt equal or bigger than allowed.
			*/
			void filterIdentificationsByRTPValues(const PeptideIdentification& identification, PeptideIdentification& 				filtered_identification, DoubleReal p_value = 0.05);

		  /**
				@brief Filters the peptide hits according to their predicted rt p-values of the first dimension
				
				Filters the peptide hits of this ProteinIdentification by the 
				probability (p-value) of a correct ProteinIdentification having a deviation between 
				observed and predicted rt equal or bigger than allowed. 
			*/
			void filterIdentificationsByRTFirstDimPValues(const PeptideIdentification& 	identification,
																										PeptideIdentification& 				filtered_identification,
																										DoubleReal 										p_value = 0.05);

      /// filters an MS/MS experiment corresponding to the threshold_fractions
			template <class PeakT>
			void filterIdentificationsByThresholds(MSExperiment< PeakT >& experiment, DoubleReal peptide_threshold_fraction, DoubleReal protein_threshold_fraction)
			{
				//filter protein hits
				ProteinIdentification temp_protein_identification;				
				std::vector<ProteinIdentification> filtered_protein_identifications;
					
				for (Size j = 0; j < experiment.getProteinIdentifications().size(); j++)
				{
					filterIdentificationsByThreshold(experiment.getProteinIdentifications()[j], protein_threshold_fraction, temp_protein_identification);
					if (!temp_protein_identification.getHits().empty())
					{
						filtered_protein_identifications.push_back(temp_protein_identification);
					}
				}
				experiment.setProteinIdentifications(filtered_protein_identifications);
				
				//filter peptide hits
				PeptideIdentification temp_identification;
				std::vector<PeptideIdentification> filtered_identifications;
				
				for (Size i = 0; i < experiment.size(); i++)
				{
					for (Size j = 0; j < experiment[i].getPeptideIdentifications().size(); j++)
					{
						filterIdentificationsByThreshold(experiment[i].getPeptideIdentifications()[j], peptide_threshold_fraction, temp_identification);
						if (!temp_identification.getHits().empty())
						{
							filtered_identifications.push_back(temp_identification);
						}
					}
					experiment[i].setPeptideIdentifications(filtered_identifications);
					filtered_identifications.clear();					
				}				
			}
																															      
      /// filters an MS/MS experiment corresponding to the threshold_fractions
			template <class PeakT>
			void filterIdentificationsByScores(MSExperiment< PeakT >& experiment, DoubleReal peptide_threshold_score, DoubleReal protein_threshold_score)
			{
				//filter protein hits
				ProteinIdentification temp_protein_identification;				
				std::vector<ProteinIdentification> filtered_protein_identifications;
					
				for (Size j = 0; j < experiment.getProteinIdentifications().size(); j++)
				{
					filterIdentificationsByScore(experiment.getProteinIdentifications()[j], protein_threshold_score, temp_protein_identification);
					if (!temp_protein_identification.getHits().empty())
					{
						filtered_protein_identifications.push_back(temp_protein_identification);
					}
				}
				experiment.setProteinIdentifications(filtered_protein_identifications);
				
				//filter peptide hits
				PeptideIdentification temp_identification;
				std::vector<PeptideIdentification> filtered_identifications;
				
				for (Size i = 0; i < experiment.size(); i++)
				{
					for (Size j = 0; j < experiment[i].getPeptideIdentifications().size(); j++)
					{
						filterIdentificationsByScore(experiment[i].getPeptideIdentifications()[j], peptide_threshold_score, temp_identification);
						if (!temp_identification.getHits().empty())
						{
							filtered_identifications.push_back(temp_identification);
						}
					}
					experiment[i].setPeptideIdentifications(filtered_identifications);
					filtered_identifications.clear();					
				}				
			}
																															      
      /// filters an MS/MS experiment corresponding to the best n hits for every spectrum
			template <class PeakT>
			void filterIdentificationsByBestNHits(MSExperiment< PeakT >& experiment, Size n)
			{
				//filter protein hits
				ProteinIdentification temp_protein_identification;				
				std::vector<ProteinIdentification> filtered_protein_identifications;
					
				for (Size j = 0; j < experiment.getProteinIdentifications().size(); j++)
				{
					filterIdentificationsByBestNHits(experiment.getProteinIdentifications()[j], n, temp_protein_identification);
					if (!temp_protein_identification.getHits().empty())
					{
						filtered_protein_identifications.push_back(temp_protein_identification);
					}
				}
				experiment.setProteinIdentifications(filtered_protein_identifications);
				
				//filter peptide hits
				PeptideIdentification temp_identification;
				std::vector<PeptideIdentification> filtered_identifications;
				
				for (Size i = 0; i < experiment.size(); i++)
				{
					for (Size j = 0; j < experiment[i].getPeptideIdentifications().size(); j++)
					{
						filterIdentificationsByBestNHits(experiment[i].getPeptideIdentifications()[j], n, temp_identification);
						if (!temp_identification.getHits().empty())
						{
							filtered_identifications.push_back(temp_identification);
						}
					}
					experiment[i].setPeptideIdentifications(filtered_identifications);
					filtered_identifications.clear();					
				}				
			}
																															      
      /// filters an MS/MS experiment corresponding to the given proteins
			template <class PeakT>
			void filterIdentificationsByProteins(MSExperiment< PeakT >& experiment, 
																					 const std::vector< FASTAFile::FASTAEntry >& proteins)
			{
				std::vector<PeptideIdentification> temp_identifications;
				std::vector<PeptideIdentification> filtered_identifications;
				PeptideIdentification temp_identification;

				for (Size i = 0; i < experiment.size(); i++)
				{
					if (experiment[i].getMSLevel() == 2)
					{
						temp_identifications = experiment[i].getPeptideIdentifications();
						for (Size j = 0; j < temp_identifications.size(); j++)
						{
							filterIdentificationsByProteins(temp_identifications[j], proteins, temp_identification);
							if (!temp_identification.getHits().empty())
							{
								filtered_identifications.push_back(temp_identification);
              }
						}
						experiment[i].setPeptideIdentifications(filtered_identifications);					
						filtered_identifications.clear();					
					}
				}				
			}
  };
 
} // namespace OpenMS

#endif // OPENMS_FILTERING_ID_IDFILTER_H
