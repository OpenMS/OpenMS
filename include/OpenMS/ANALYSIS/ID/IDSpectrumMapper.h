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

#ifndef OPENMS_ANALYSIS_ID_IDSPECTRUMMAPPER_H
#define OPENMS_ANALYSIS_ID_IDSPECTRUMMAPPER_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Identification.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Annotate a MSExperiment instances with Identification instances
    
    The identifications stored in a Identification instance can be added to the
    corresponding spectrum. Furthermore the annotations that are present
    can be retrieved.
  */
  class IDSpectrumMapper
  {
    public:

      /// Constructor
      IDSpectrumMapper();
      
			/**
				@brief Annotates the spectra belonging to the experiment
				
				The retention times and mz-values are used to find the 
				spectrum in experiment to which the corresponding 
				identification belongs. The Identification is then added to the spectrum      					
      */
      template <class PeakT>				
      UnsignedInt annotate(MSExperiment< PeakT >& experiment, const std::vector<IdentificationData>& identifications, float precision = 0.01f)
  		{
				std::multimap<float, UnsignedInt> experiment_precursors;
				std::multimap<float, UnsignedInt> identifications_precursors;
				std::multimap<float, UnsignedInt>::iterator experiment_iterator;
				std::multimap<float, UnsignedInt>::iterator identifications_iterator;
				float temp_experiment_value;
				float temp_identification_value;
				DPosition< 1 >::CoordinateType experiment_precursor_position;
				DPosition< 1 >::CoordinateType identifications_precursor_position;
				float temp;
				float actual_retention_time = 0;
				UnsignedInt counter = 0;
					
				for(UnsignedInt i = 0; i < experiment.size(); i++)
				{
					actual_retention_time = experiment[i].getRetentionTime();
					experiment_precursors.insert(std::make_pair(actual_retention_time, i));
				}
				for(UnsignedInt i = 0; i < identifications.size(); i++)
				{
					identifications_precursors.insert(std::make_pair(identifications[i].rt, i));
				}
				experiment_iterator = experiment_precursors.begin();
				identifications_iterator = identifications_precursors.begin();	
				
				while(experiment_iterator != experiment_precursors.end()
							&& identifications_iterator != identifications_precursors.end())
				{
					temp_experiment_value = experiment_iterator->first;
					while(identifications_iterator != identifications_precursors.end())
					{
						temp_identification_value = identifications_iterator->first;
						temp = temp_experiment_value - temp_identification_value;
						
						/// testing whether the retention times are within the precision threshold
						if (((temp < precision) && temp >= 0) || (((-1 * temp) < precision) && temp < 0))
						{
							experiment_precursor_position = experiment[experiment_iterator->second].getPrecursorPeak().getPosition()[0];
							identifications_precursor_position = identifications[identifications_iterator->second].mz;				
							temp = identifications_precursor_position - experiment_precursor_position;
							if (((temp < precision) && temp >= 0) || ((-1 * temp < precision) && temp < 0))
							{
								if (!identifications[identifications_iterator->second].id.empty())
								{
									experiment[experiment_iterator->second].getIdentifications().push_back(identifications[identifications_iterator->second].id);
									counter++;
								}
							}
						}
						++identifications_iterator;
					}
					identifications_iterator = identifications_precursors.begin();	
					++experiment_iterator;
				}
				return counter;  			
  		}
      				 
			/**
				@brief stores all annotations of the spectra in the specific member variables							
      	
				Every non-empty Identification that is associated with a spectrum is
				stored as well as the retention times and the mz-values.				
      */				
      template <class PeakT>				
      void getAnnotations(const MSExperiment< PeakT >& experiment, std::vector<IdentificationData>& identifications)
      {
				identifications.clear();
				std::vector<Identification> temp_identifications;
				IdentificationData tmp_id;
				for(UnsignedInt i = 0; i < experiment.size(); i++)
				{
					temp_identifications = experiment[i].getIdentifications();
					if (temp_identifications.size() > 0)
					{
						tmp_id.rt = experiment[i].getRetentionTime();				
						tmp_id.mz = experiment[i].getPrecursorPeak().getPosition()[0];
						for(UnsignedInt j = 0; j < temp_identifications.size(); j++)
						{
							if (!temp_identifications[j].empty())
							{
								tmp_id.id = temp_identifications[j];
								identifications.push_back(tmp_id);
							}
						}
					}
				}							      	
      }
      
    
    protected:
      
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDSPECTRUMMAPPER_H
