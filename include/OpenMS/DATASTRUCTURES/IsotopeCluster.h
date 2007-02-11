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

#ifndef OPENMS_DATASTRUCTURES_ISOTOPECLUSTER_H
#define OPENMS_DATASTRUCTURES_ISOTOPECLUSTER_H

#include <OpenMS/CONCEPT/Types.h>

#include <set>

namespace OpenMS
{	
	///Stores information about an isotopic cluster (i.e. potential peptide charge variants)
  struct IsotopeCluster
  {
  	//An index in an MSExperiment
  	typedef std::pair<UnsignedInt,UnsignedInt> IDX;
  	
    IsotopeCluster()
      : charge_(0), 
      	peaks_(), 
      	scans_()
    {
    }
    
    /// predicted charge state of this peptide
    UnsignedInt charge_;
    
    /// peaks in this cluster
    std::set<IDX> peaks_;
    
    /// the scans of this cluster
    std::vector<UnsignedInt> scans_;
  };

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_ISOTOPECLUSTER_H
