// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_ISOTOPECLUSTER_H
#define OPENMS_DATASTRUCTURES_ISOTOPECLUSTER_H

#include <OpenMS/CONCEPT/Types.h>
#include <vector>
#include <set>

namespace OpenMS
{
	///Stores information about an isotopic cluster (i.e. potential peptide charge variants)
  struct OPENMS_DLLAPI IsotopeCluster
  {
  	/// An index e.g. in an MSExperiment
  	typedef std::pair<Size,Size> IndexPair;
		/// A set of index pairs, usually referring to an MSExperiment.
		typedef std::set<IndexPair> IndexSet;

		///index set with associated charge estimate
		struct ChargedIndexSet
			: public IndexSet
		{
			ChargedIndexSet()
			: charge(0)
			{
			}
			/// charge estimate (convention: zero means "no charge estimate")
			Int charge;
		};

    IsotopeCluster()
      : peaks(),
      	scans()
    {
    }
    /// peaks in this cluster
    ChargedIndexSet peaks;

    /// the scans of this cluster
    std::vector<Size> scans;
  };

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_ISOTOPECLUSTER_H
