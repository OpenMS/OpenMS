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
// $Maintainer: Chris Bielow$
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
 
#ifndef OPENMS_CHEMISTRY_WEIGHTWRAPPER_H
#define OPENMS_CHEMISTRY_WEIGHTWRAPPER_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/AASequence.h> 
 
namespace OpenMS
{
 
	/**
	@brief Encapsulated weight queries to simplify mono vs average weight computation
	
	Supports EmpiricalFormula's and AASequence's getMonoWeight() and getAverageWeight()
	
	*/
	class OPENMS_DLLAPI WeightWrapper
	{
 	
		public:
				
			enum WEIGHTMODE {AVERAGE=0, MONO, SIZE_OF_WEIGHTMODE};
			
			/**
			@brief constructor
			*/
			WeightWrapper();

			/**
			@brief constructor
			*/
			WeightWrapper(const WEIGHTMODE weight_mode);
			
			/**
			@brief destructor
			*/
			virtual ~WeightWrapper();
			
			/**
			@brief copy constructor
			*/
			WeightWrapper(const WeightWrapper & source);
		
		
			/**
			@brief Sets the weight mode (MONO or AVERAGE)

			Sets the mode in which getWeight() calls are answered.
			
			*/
			void setWeightMode(const WEIGHTMODE mode);
		

			/**
			@brief Gets the weight mode (MONO or AVERAGE)

			Gets the mode in which getWeight() calls are answered.
			
			*/
			WEIGHTMODE getWeightMode() const;

		
			/**
			@brief returns the weight of either mono or average value

			Which weight is returned depends on the current weight-mode.
			
			@return DoubleReal weight in u
			*/
			DoubleReal getWeight(const AASequence& aa) const;
			
			/**
			@brief returns the weight of either mono or average value

			Which weight is returned depends on the current weight-mode.
			
			@return DoubleReal weight in u
			*/
			DoubleReal getWeight(const EmpiricalFormula& ef) const;


			/**
			@brief returns the weight of either mono or average value

			Which weight is returned depends on the current weight-mode.
			
			@return DoubleReal weight in u
			*/
			DoubleReal getWeight(const Residue& r, Residue::ResidueType res_type = Residue::Full) const;
			

		private:
			
			WEIGHTMODE weight_mode_; ///< one of WeightWrapper::WEIGHTMODE's values 

		
	};
}
#endif // OPENMS_CHEMISTRY_WeightWrapper_H
