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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DECHARGING_ILPDCWRAPPER_H
#define OPENMS_ANALYSIS_DECHARGING_ILPDCWRAPPER_H

#include <vector>
#include <OpenMS/DATASTRUCTURES/ChargePair.h>


namespace OpenMS {

	class MassExplainer;
	
	class OPENMS_DLLAPI ILPDCWrapper
	{

	public:
		typedef std::vector<ChargePair> PairsType;
		typedef PairsType::size_type PairsIndex;
		
		///Constructor
		ILPDCWrapper();

		///Destructor
		virtual ~ILPDCWrapper();

		/// Compute optimal solution and return value of objective function
    /// If the input feature map is empty, a warning is issued and -1 is returned.
		/// @return value of objective function
		///		and @p pairs will have all realized edges set to "active"
		DoubleReal compute(const FeatureMap<> fm, PairsType& pairs, Size verbose_level);
		
	private:
    void 
      writeProblem_(const FeatureMap<> fm,
      PairsType& pairs, 
      const PairsIndex margin_left, 
      const PairsIndex margin_right,
      Size verbose_level,
      String filename);

    /// slicing the problem into subproblems
		DoubleReal computeSlice_(const FeatureMap<> fm,
														 PairsType& pairs, 
														 const PairsIndex margin_left, 
														 const PairsIndex margin_right,
                             Size verbose_level);

		/// calculate a score for the i_th edge
		DoubleReal getLogScore_(const PairsType::value_type& pair, const FeatureMap<>& fm);

	}; // !class
    
} // !namespace

#endif // OPENMS_ANALYSIS_DECHARGING_ILPDCWRAPPER_H

