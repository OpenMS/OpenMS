// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DECHARGING_MIPWRAPPER_H
#define OPENMS_ANALYSIS_DECHARGING_MIPWRAPPER_H

#include <vector>
#include <OpenMS/DATASTRUCTURES/ChargePair.h>


namespace OpenMS {

	class MassExplainer;
	
	class OPENMS_DLLAPI MIPWrapper
	{

	public:
		typedef std::vector<ChargePair> PairsType;
		typedef PairsType::size_type PairsIndex;
		
		///Constructor
		MIPWrapper();

		///Destructor
		virtual ~MIPWrapper();

		/// compute optimal solution and return value of objective function
		/// @return value of objective function
		///		and @p pairs will have all realized edges set to "active"
		DoubleReal compute(const MassExplainer& me, const FeatureMap<> fm, PairsType& pairs);
		
	private:
		/// slicing the problem into subproblems
		DoubleReal computeSlice_(const MassExplainer& me,
														 const FeatureMap<> fm,
														 PairsType& pairs, 
														 const PairsIndex margin_left, 
														 const PairsIndex margin_right);

		/// calculate a score for the i_th edge
		DoubleReal getLogScore_(const PairsIndex& i, const PairsType& pairs, const FeatureMap<>& fm);

	}; // !class
    
} // !namespace

#endif // OPENMS_ANALYSIS_DECHARGING_MIPWrapper_H

