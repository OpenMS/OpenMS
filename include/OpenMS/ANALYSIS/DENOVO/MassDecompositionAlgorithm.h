// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYIS_ID_DECOMPOSITIONALGORITHM_H
#define OPENMS_ANALYIS_ID_DECOMPOSITIONALGORITHM_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include "MassDecomposition.h"

// ims includes
#ifdef OPENMS_COMPILER_MSVC
	#pragma warning( push )
	#pragma warning( disable : 4290 4267)
#endif

#include <ims/decomp/realmassdecomposer.h>
#include <ims/alphabet.h>
#include <ims/weights.h>
#ifdef OPENMS_COMPILER_MSVC
	#pragma warning( pop )
#endif

#include <vector>

namespace OpenMS
{
	/** @brief Mass decomposition algorithm, given a mass it suggests possible compositions

			A mass decomposition algorithm decomposes a mass or a mass difference into
			possible amino acids and frequencies of them, which add up to the given mass.
			This class is a wrapper for the algorithm published in 
			@todo add reference
	*/
	class OPENMS_DLLAPI MassDecompositionAlgorithm : public DefaultParamHandler
	{
		public:
			
			/** @name constructors and desctructor
			*/
			//@{
			/// Default constructor
			MassDecompositionAlgorithm();

			/// Copy constructor
			MassDecompositionAlgorithm(const MassDecompositionAlgorithm& deco);

			/// Destructor
			virtual ~MassDecompositionAlgorithm();
			//@}

			/** @name Operators
			*/
			//@{
			/// assignment operator
			MassDecompositionAlgorithm& operator = (const MassDecompositionAlgorithm& rhs);

			/// returns the possible decompositions given the weight
			void getDecompositions(std::vector<MassDecomposition>& decomps, DoubleReal weight);
			//@}

		protected:

			void updateMembers_();

			ims::Alphabet* alphabet_;

			ims::RealMassDecomposer* decomposer_;
	};

} // namespace OpenMS

#endif
