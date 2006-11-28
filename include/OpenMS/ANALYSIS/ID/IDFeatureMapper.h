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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_IDFEATUREMAPPER_H
#define OPENMS_ANALYSIS_ID_IDFEATUREMAPPER_H

#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/METADATA/Identification.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Annotates a DPEatureMap instance with Identification instances
    
    The identifications stored in a Identification instance can be added to the
    corresponding spectrum. Furthermore the annotations that are present
    can be retrieved.
  */
  class IDFeatureMapper
  {
    public:

      /// Constructor
      IDFeatureMapper();
      
			/**
				@brief Annotates the spectra belonging to the experiment
				
				The retention times and mz-values are used to find the 
				spectrum in experiment to which the corresponding 
				identification belongs. The Identification is then added to the spectrum      					
      */
      template <class PeakT>				
      UnsignedInt annotate(DFeatureMap<2> fm,
 			      					 	   const std::vector<Identification>& identifications,
			      				 			 const std::vector<float>& precursor_retention_times,
			      				 			 const std::vector<float>& precursor_mz_values,
			  				 				   float precision = 0.01f)
  		{
				
  		}
      
    protected:
      
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDFEATUREMAPPER_H
