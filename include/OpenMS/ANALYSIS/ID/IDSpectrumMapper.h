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
				
		 		The retention time and mass-to-charge ratio of the PeptideIdentification have to 
		 		be given in the MetaInfoInterface ('MZ' and 'RT').   					
      	
      	The exception MissingInformation is thrown if the MetaInfoInterface of @p identifications does not contain 'MZ' and 'RT'.
      */
      template <class PeakT>				
      UInt annotate(MSExperiment< PeakT >& experiment, const std::vector<PeptideIdentification>& identifications, DoubleReal precision = 0.01f) throw (Exception::MissingInformation)
  		{
				std::multimap<DoubleReal, UInt> experiment_precursors;
				std::multimap<DoubleReal, UInt> identifications_precursors;
				std::multimap<DoubleReal, UInt>::iterator experiment_iterator;
				std::multimap<DoubleReal, UInt>::iterator identifications_iterator;
				DoubleReal temp_experiment_value;
				DoubleReal temp_identification_value;
				DoubleReal experiment_precursor_position;
				DoubleReal identifications_precursor_position;
				DoubleReal temp;
				DoubleReal actual_retention_time = 0;
				UInt counter = 0;
					
				for(UInt i = 0; i < experiment.size(); i++)
				{
					actual_retention_time = experiment[i].getRT();
					experiment_precursors.insert(std::make_pair(actual_retention_time, i));
				}
				for(UInt i = 0; i < identifications.size(); i++)
				{
					if (!identifications[i].metaValueExists("RT"))
					{
						throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDSpectrumMapper: MetaValue 'RT' missing!"); 
					}

					if (!identifications[i].metaValueExists("MZ"))
					{
						throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDSpectrumMapper: MetaValue 'MZ' missing!"); 
					}
					identifications_precursors.insert(std::make_pair(identifications[i].getMetaValue("RT"), i));
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
							identifications_precursor_position = identifications[identifications_iterator->second].getMetaValue("MZ");				
							temp = identifications_precursor_position - experiment_precursor_position;
							if (((temp < precision) && temp >= 0) || ((-1 * temp < precision) && temp < 0))
							{
								if (!(identifications[identifications_iterator->second].empty()))
								{
									experiment[experiment_iterator->second].getPeptideIdentifications().push_back(identifications[identifications_iterator->second]);
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
    
    protected:
      
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDSPECTRUMMAPPER_H
