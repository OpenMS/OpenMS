// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>
#include <cmath>

namespace OpenMS 
{
  /**
    @brief Annotate a MSExperiment instances with ProteinIdentification instances
    
    The identifications stored in a ProteinIdentification instance can be added to the
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
  			//store mapping of scan RT to index
				std::multimap<DoubleReal, UInt> experiment_precursors;
				for(UInt i = 0; i < experiment.size(); i++)
				{
					experiment_precursors.insert(std::make_pair(experiment[i].getRT(), i));
				}
				
				//store mapping of identification RT to index
				std::multimap<DoubleReal, UInt> identifications_precursors;
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
				
				//calculate the actual mapping
				std::multimap<DoubleReal, UInt>::iterator experiment_iterator = experiment_precursors.begin();
				std::multimap<DoubleReal, UInt>::iterator identifications_iterator = identifications_precursors.begin();	
				UInt counter = 0;				
				while(experiment_iterator != experiment_precursors.end() && identifications_iterator != identifications_precursors.end())
				{
					while(identifications_iterator != identifications_precursors.end())
					{
						// testing whether the retention times are within the precision threshold
						if (fabs(experiment_iterator->first - identifications_iterator->first) < precision)
						{
							//std::cout << "RT matching (scan/id) " << experiment_iterator->first << " / " << identifications_iterator->first << std::endl;	
							
							// testing wheather the m/z fits
							if (fabs((DoubleReal)(identifications[identifications_iterator->second].getMetaValue("MZ")) -  (DoubleReal)(experiment[experiment_iterator->second].getPrecursorPeak().getPosition()[0])) < precision)
							{
								//std::cout << "MZ matching (scan/id) " << (DoubleReal)(experiment[experiment_iterator->second].getPrecursorPeak().getPosition()[0]) << " / " << (DoubleReal)(identifications[identifications_iterator->second].getMetaValue("MZ")) << std::endl;	
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
