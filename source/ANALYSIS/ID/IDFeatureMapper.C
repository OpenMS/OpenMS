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

#include <OpenMS/ANALYSIS/ID/IDFeatureMapper.h>

namespace OpenMS 
{
  IDFeatureMapper::IDFeatureMapper() 
  {
    
  }
 
  void IDFeatureMapper::annotate(DFeatureMap<2> fm, const std::vector<Identification>& identifications, const std::vector<float>& precursor_retention_times, const std::vector<float>& precursor_mz_values)
	{
		for(DFeatureMap<2>::Iterator it = fm.begin(); it!=fm.end(); ++it)
		{
			
		}
	}
} // namespace OpenMS
