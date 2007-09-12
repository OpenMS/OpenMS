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
#include <vector>
#include <set>

namespace OpenMS
{	
	///Stores information about an isotopic cluster (i.e. potential peptide charge variants)
  struct IsotopeCluster
  {
  	/// An index in an MSExperiment
  	typedef std::pair<UInt,UInt> IDX;
		/// Int set
		typedef std::set<IDX> IndexSet;
		/// Charge score
		typedef DoubleReal ProbabilityType;	
		
		/**@brief index set with associated charge estimate
		 *
		 *@todo support a list of charge predictions (FF20)
		 */
		struct ChargedIndexSet 
		: public IndexSet
		{
			ChargedIndexSet() : charge_(0), max_charge_score_() { }
			/// charge estimate (convention: zero means "no charge estimate")
			UInt charge_;			
			/// Score of highest scoring charge state
			ProbabilityType max_charge_score_; 
		};
		  	
    IsotopeCluster()
      : peaks_(), 
      	scans_()
    {
    }
    /// peaks in this cluster
    ChargedIndexSet peaks_;
    
    /// the scans of this cluster
    std::vector<UInt> scans_;
  };

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_ISOTOPECLUSTER_H
