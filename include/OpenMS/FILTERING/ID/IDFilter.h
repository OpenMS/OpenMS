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

#ifndef OPENMS_FILTERING_ID_IDFILTER_H
#define OPENMS_FILTERING_ID_IDFILTER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief used to filter identifications by different criteria
    
    The identifications are filtered by significance thresholds and
    by sequences. The filtering by significance thresholds looks for the 
    best identification that fullfills the significance threshold criterium.
    score > significance-threshold * significance_fraction. 
    The filtering by sequences looks for the best identification that
    is contained in one of the protein sequences.
  */
  class IDFilter
  {
    public:

      /// Constructor
      IDFilter();

      /// Copy constructor
      IDFilter(const IDFilter& source);

      /// Destructor
      ~IDFilter();
      
      /// Assignment operator
      IDFilter& operator = (const IDFilter& source);
      
      /// returns a const reference to the peptide threshold fraction
      const double& getPeptideThresholdFraction() const;

      /// sets the peptide threshold fraction
      void setPeptideThresholdFraction(const double& peptide_threshold_fraction);

      /// returns a const reference to the protein threshold fraction
      const double& getProteinThresholdFraction() const;

      /// sets the protein threshold fraction
      void setProteinThresholdFraction(const double& protein_threshold_fraction);

      /// returns a const reference to the protein sequences used for filtering
      const std::vector< std::pair<String, String> >& getProteins() const;

      /// sets the protein sequences used for filtering
      void setProteins(const std::vector< std::pair<String, String> >& proteins);

      /// filters an identification corresponding to the threshold_fractions
			void filterIdentificationsByThresholds(const Identification& identification, const double& peptide_threshold_fraction, const double& protein_threshold_fraction, Identification& filtered_identification);

      /// filters an identification corresponding to the threshold_fractions
			void filterIdentificationsByThresholds(const Identification& identification, Identification& filtered_identification);
			
      /// filters an identification keeping only the best scoring hits (if strict is set, keeping only the best hit only if it is the only hit with that score)
			void filterIdentificationsByBestHits(const Identification& identification, Identification& filtered_identification, bool strict = false);

      /// filters an identification corresponding to the given proteins
			void filterIdentificationsByProteins(const Identification& identification, std::vector< std::pair<String, String> > proteins, Identification& filtered_identification);

      /// filters an identification corresponding to the given proteins
			void filterIdentificationsByProteins(const Identification& identification, Identification& filtered_identification);

      /// filters a protein identification corresponding to the given proteins only proteins with the same accession
			void filterIdentificationsByProteins(const ProteinIdentification&     identification, 
                                           ProteinIdentification&           filtered_identification);
        
      /// filters the peptides of a identification corresponding to the retention times
			void filterIdentificationsByRetentionTimes(const Identification& identification, const std::map<String, double>& predicted_retention_times, double measured_retention_time, double predicted_sigma, double allowed_deviation, double total_gradient_time, Identification& filtered_identification);
																														
			/// removes all peptide hits having a sequence equal to a element in <code>peptides</code>
			void filterIdentificationsByExclusionPeptides(const Identification& identification, std::vector<String> peptides, Identification& filtered_identification);
																														

      /// filters an MS/MS experiment corresponding to the threshold_fractions
			template <class PeakT>
			void filterIdentificationsByThresholds(MSExperiment< PeakT >& experiment,
																						 double peptide_threshold_fraction,
																						 double protein_threshold_fraction)
			{
				std::vector<Identification> temp_identifications;
				std::vector<Identification> filtered_identifications;
				Identification temp_identification;
				
				peptide_threshold_fraction_ = peptide_threshold_fraction;
				protein_threshold_fraction_ = protein_threshold_fraction;		
				
				for(unsigned int i = 0; i < experiment.size(); i++)
				{		
					
					if (experiment[i].getMSLevel() == 2)
					{
						temp_identifications = experiment[i].getIdentifications();
						if (temp_identifications.size() > 0)
						{
							for(UnsignedInt j = 0; j < temp_identifications.size(); j++)
							{
								filterIdentificationsByThresholds(temp_identifications[j], temp_identification);
								if (!temp_identification.empty())
								{
									filtered_identifications.push_back(temp_identification);
								}
							}
							experiment[i].setIdentifications(filtered_identifications);
							filtered_identifications.clear();					
						}
					}
				}				
			}
																															      
      /// filters an MS/MS experiment corresponding to the given proteins
			template <class PeakT>
			void filterIdentificationsByProteins(MSExperiment< PeakT >& experiment, 
																					 std::vector< std::pair<String, String> >proteins)
			{
				std::vector<Identification> temp_identifications;
				std::vector<Identification> filtered_identifications;
				Identification temp_identification;
				
				setProteins(proteins);
		
				for(unsigned int i = 0; i < experiment.size(); i++)
				{		
					
					if (experiment[i].getMSLevel() == 2)
					{
						temp_identifications = experiment[i].getIdentifications();
						if (temp_identifications.size() > 0)
						{
							for(UnsignedInt j = 0; j < temp_identifications.size(); j++)
							{
								filterIdentificationsByProteins(temp_identifications[j], temp_identification);
								if (!temp_identification.empty())
								{
									filtered_identifications.push_back(temp_identification);
                }
							}
							experiment[i].setIdentifications(filtered_identifications);					
							filtered_identifications.clear();					
						}
					}
				}				
			}
			
    protected:
      double peptide_threshold_fraction_;			/// the fraction of the significance threshold a score should have to be kept
      double protein_threshold_fraction_;			/// the fraction of the significance threshold a score should have to be kept
      std::vector< std::pair <String, String> > proteins_; /// the proteins (restriction of result space) 
  };
 
} // namespace OpenMS

#endif // OPENMS_FILTERING_ID_IDFILTER_H
