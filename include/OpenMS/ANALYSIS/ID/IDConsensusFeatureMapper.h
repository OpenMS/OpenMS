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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_IDCONSENSUSFEATUREMAPPER_H
#define OPENMS_ANALYSIS_ID_IDCONSENSUSFEATUREMAPPER_H

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Annotates a ConsensusMap instance with PeptideIdentification instances
 		
		The mapping of PeptideIdentification to ConsensusFeatures is done by comparing the MZ and RT
	  coordinates, allowing for a small deviation.
	  ProteinIdentifications are assigned to the whole ConsensusMap
	 
	  The retention time and mass-to-charge ratio of the PeptideIdentification have to 
 		be given in the MetaInfoInterface ('MZ' and 'RT').
	 
  */
  class IDConsensusFeatureMapper
  {
    public:

			typedef Peak2D::CoordinateType CoordinateType;
			
      /// Default constructor
      IDConsensusFeatureMapper();
      
			/**
				@brief This method does the actual mapping
			
			  @param cm ConsensusMap to receive the identifications
			  @param ids PeptideIdentification for the ConsensusFeatures
			  @param ids ProteinIdentification for the ConsensusMap
			  @param mz_delta Allowed m/z deviation
			  @param rt_delta Allowed RT deviation
			  @param measure_from_subelements Do distance estimate from FeatureHandles instead of Centroid
			 
				@exception Exception::MissingInformation is thrown if the MetaInfoInterface of @p ids does not contain 'MZ' and 'RT'
			*/
		  void annotate(ConsensusMap& cm,
										const std::vector<PeptideIdentification>& ids,
										const std::vector<ProteinIdentification>& protein_ids,
										CoordinateType mz_delta=0.05,
										CoordinateType rt_delta=0.5,
										bool measure_from_subelements=false);      
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDCONSENSUSFEATUREMAPPER_H
