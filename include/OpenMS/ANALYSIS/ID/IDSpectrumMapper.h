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
    @brief Annotates the spectra of an MSExperiment with peptide identifications
    
    The identifications stored in a PeptideIdentification instance can be added to the
    corresponding spectrum. Furthermore the annotations that are present
    can be retrieved.

	  ProteinIdentifications are assigned to the whole map.
	  
 		The retention time and mass-to-charge ratio of the PeptideIdentification have to 
 		be given in the MetaInfoInterface ('MZ' and 'RT').
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
      	
			  @param map MSExperiment to receive the identifications
			  @param ids PeptideIdentification for the ConsensusFeatures
			  @param protein_ids ProteinIdentification for the ConsensusMap
			  @param mz_delta Allowed m/z deviation
			  @param rt_delta Allowed RT deviation
      	
				@throws Exception::MissingInformation if the MetaInfoInterface of @p ids does not contain 'MZ' and 'RT'.
      */
      template <typename PeakType>				
      void annotate(MSExperiment< PeakType >& map,
                    const std::vector<PeptideIdentification>& ids,
		                const std::vector<ProteinIdentification>& protein_ids,
                    DoubleReal rt_delta = 0.01,
                    DoubleReal mz_delta = 0.01)
  		{
				//append protein identifications
				map.getProteinIdentifications().insert(map.getProteinIdentifications().end(),protein_ids.begin(),protein_ids.end());
		
  			//store mapping of scan RT to index
				std::multimap<DoubleReal, UInt> experiment_precursors;
				for(UInt i = 0; i < map.size(); i++)
				{
					experiment_precursors.insert(std::make_pair(map[i].getRT(), i));
				}
				
				//store mapping of identification RT to index
				std::multimap<DoubleReal, UInt> identifications_precursors;
				for(UInt i = 0; i < ids.size(); i++)
				{
					if (!ids[i].metaValueExists("RT"))
					{
						throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDSpectrumMapper: meta data value 'RT' missing for peptide identification!"); 
					}

					if (!ids[i].metaValueExists("MZ"))
					{
						throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDSpectrumMapper: meta data value 'MZ' missing for peptide identification!"); 
					}
					identifications_precursors.insert(std::make_pair(ids[i].getMetaValue("RT"), i));
				}
				
				//calculate the actual mapping
				std::multimap<DoubleReal, UInt>::iterator experiment_iterator = experiment_precursors.begin();
				std::multimap<DoubleReal, UInt>::iterator identifications_iterator = identifications_precursors.begin();	
				while(experiment_iterator != experiment_precursors.end() && identifications_iterator != identifications_precursors.end())
				{
					while(identifications_iterator != identifications_precursors.end())
					{
						// testing whether the retention times are within the precision threshold
						if (fabs(experiment_iterator->first - identifications_iterator->first) < rt_delta)
						{
							//std::cout << "RT matching (scan/id) " << experiment_iterator->first << " / " << identifications_iterator->first << std::endl;	
							
							// testing whether the m/z fits
							if (fabs((DoubleReal)(ids[identifications_iterator->second].getMetaValue("MZ")) -  (DoubleReal)(map[experiment_iterator->second].getPrecursorPeak().getPosition()[0])) < mz_delta)
							{
								//std::cout << "MZ matching (scan/id) " << (DoubleReal)(map[experiment_iterator->second].getPrecursorPeak().getPosition()[0]) << " / " << (DoubleReal)(ids[identifications_iterator->second].getMetaValue("MZ")) << std::endl;	
								if (!(ids[identifications_iterator->second].empty()))
								{
									map[experiment_iterator->second].getPeptideIdentifications().push_back(ids[identifications_iterator->second]);
								}
							}
						}
						++identifications_iterator;
					}
					identifications_iterator = identifications_precursors.begin();	
					++experiment_iterator;
				}
  		}
    
    protected:
      
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDSPECTRUMMAPPER_H
