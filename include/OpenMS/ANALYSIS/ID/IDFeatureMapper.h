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
    @brief Annotates a DFeatureMap instance with Identification instances
 		
 		
  */
  class IDFeatureMapper
  {
    public:

      /// Default constructor
      IDFeatureMapper();
      
			///This method does the actual mapping
		  void annotate(DFeatureMap<2>& fm, const std::vector<IdentificationData>& ids, const std::vector<ProteinIdentification>& protein_ids) throw (Exception::Precondition);      
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDFEATUREMAPPER_H
