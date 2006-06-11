// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: MSExperimentAnnotator.h,v 1.8 2006/06/10 06:40:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_MSEXPERIMENTANNOTATOR_H
#define OPENMS_ANALYSIS_ID_MSEXPERIMENTANNOTATOR_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Identification.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Annotates MSExperiments with Identification instances
    
    The identifications stored in a Identification instance can be added to the
    corresponding spectrum. Furthermore the annotations that are present
    can be retrieved.
  */
  class MSExperimentAnnotator
  {
    public:

      /// Constructor
      MSExperimentAnnotator();
      /// Copy constructor
      MSExperimentAnnotator(const MSExperimentAnnotator& source);
      /// Destructor
      ~MSExperimentAnnotator();
      
      /// Assignment operator
      MSExperimentAnnotator& operator = (const MSExperimentAnnotator& source);

			/**
				@brief Annotates the spectra belonging to the experiment
				
				The retention times and mz-values are used to find the 
				spectrum in experiment to which the corresponding 
				identification belongs. The Identification is then added to the spectrum      					
      */
      template <class PeakT>				
      UnsignedInt annotate(MSExperiment< PeakT >& experiment,
 			      					 	   const std::vector<Identification>& identifications,
			      				 			 const std::vector<float>& precursor_retention_times,
			      				 			 const std::vector<float>& precursor_mz_values,
			  				 				   float precision = 0.01f)
  		{
				std::multimap<float, UnsignedInt> experiment_precursors;
				std::multimap<float, UnsignedInt> identifications_precursors;
				std::multimap<float, UnsignedInt>::iterator experiment_iterator;
				std::multimap<float, UnsignedInt>::iterator identifications_iterator;
				float temp_experiment_value;
				float temp_db_search_value;
				DPosition< 1 >::CoordinateType experiment_precursor_position;
				DPosition< 1 >::CoordinateType identifications_precursor_position;
				std::vector<Identification> temp_identifications;
				float temp;
				float actual_retention_time = 0;
				UnsignedInt counter = 0;
					
				for(UnsignedInt i = 0; i < experiment.size(); i++)
				{
					actual_retention_time = experiment[i].getRetentionTime();
					experiment_precursors.insert(std::make_pair(actual_retention_time, i));
				}
				for(UnsignedInt i = 0; i < precursor_retention_times.size(); i++)
				{
					identifications_precursors.insert(std::make_pair(precursor_retention_times[i], i));
				}
				experiment_iterator = experiment_precursors.begin();
				identifications_iterator = identifications_precursors.begin();
				while(experiment_iterator != experiment_precursors.end()
							&& identifications_iterator != identifications_precursors.end())
				{
					temp_experiment_value = experiment_iterator->first;
					temp_db_search_value = identifications_iterator->first;
					temp = temp_experiment_value - temp_db_search_value;
					
					/// testing whether the retention times are within the precision threshold
					if (((temp < precision) && temp >= 0)
						 || ((-1 * temp < precision) && (-1 * temp >= 0)))
					{
						experiment_precursor_position = 
							experiment[experiment_iterator->second].getPrecursorPeak().getPosition()[0];
						identifications_precursor_position = precursor_mz_values[identifications_iterator->second];				
						temp = identifications_precursor_position - experiment_precursor_position;
						if (((temp < precision) && temp >= 0)
						 || ((-1 * temp < precision) && (-1 * temp >= 0)))
						{
							if (!identifications[identifications_iterator->second].empty())
							{
								temp_identifications = experiment[experiment_iterator->second].getIdentification();
								temp_identifications.push_back(identifications[identifications_iterator->second]);
								experiment[experiment_iterator->second].setIdentification(temp_identifications);
								counter++;
							}
						}
						identifications_iterator++;
					}
					else if (temp_experiment_value < temp_db_search_value)
					{
						experiment_iterator++;
					}
					else
					{
						identifications_iterator++;
					}
				}
				return counter;  			
  		}
      				 
			/**
				@brief stores all annotations of the spectra in the specific member
				variables							
      	
				Every non-empty Identification that is associated with a spectrum is
				stored as well as the retention times and the mz-values.				
      */				
      template <class PeakT>				
      void getAnnotations(const MSExperiment< PeakT >& experiment,
      						 			  std::vector<Identification>* identifications,
      				       			std::vector<float>* precursor_retention_times,
      				       			std::vector<float>* precursor_mz_values)
      {
				*identifications = std::vector<Identification>();  					 // clear information
				*precursor_retention_times = std::vector<float>();  // clear information
				*precursor_mz_values = std::vector<float>();  			 // clear information
				
				std::vector<Identification> temp_identifications;
				float temp_precursor_retention_time;
				float temp_precursor_mz_value;
				
				for(UnsignedInt i = 0; i < experiment.size(); i++)
				{
					temp_identifications = experiment[i].getIdentification();
					if (temp_identifications.size() > 0)
					{
						temp_precursor_retention_time = experiment[i].getRetentionTime();				
						temp_precursor_mz_value = experiment[i].getPrecursorPeak().getPosition()[0];
						for(UnsignedInt j = 0; j < temp_identifications.size(); j++)
						{
							if (!temp_identifications[j].empty())
							{
								identifications->push_back(temp_identifications[j]);
								precursor_retention_times->push_back(temp_precursor_retention_time);
								precursor_mz_values->push_back(temp_precursor_mz_value);
							}
						}
					}
				}							      	
      }
      
    
    protected:
      
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_MSEXPERIMENTANNOTATOR_H
