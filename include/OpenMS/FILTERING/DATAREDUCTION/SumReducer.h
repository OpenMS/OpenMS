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

#ifndef OPENMS_FILTERING_DATAREDUCTION_SUMREDUCER_H
#define OPENMS_FILTERING_DATAREDUCTION_SUMREDUCER_H

#include<OpenMS/FILTERING/DATAREDUCTION/DataReducer.h>

namespace OpenMS
{
	/**
		@brief Reduces the amount of data in an experiment by summing the
		intensities of neighboring peaks.
		 
		@ref SumReducer_Parameters are explained on a separate page.
	*/
	class SumReducer
	  : public DataReducer
	{
		public:
		 	///constructor
			SumReducer();
	
			///destructor
			~SumReducer();
	
			/// calculates the reduced MSExperiment
	 		virtual void applyReduction(const ExperimentType& in,ExperimentType& out );
			
			/// returns an instance of this class
			static DataReducer* create()
			{	
				return new SumReducer();
			}
	
			/// returns the name of this module
			static const String getProductName()
			{
				return "sum_reducer";
			}
	};
}
#endif
