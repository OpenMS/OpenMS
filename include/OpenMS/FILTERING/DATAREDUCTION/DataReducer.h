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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_DATAREDUCTION_DATAREDUCER_H
#define OPENMS_FILTERING_DATAREDUCTION_DATAREDUCER_H

#include<OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{	
/** 
	@brief Abstract base class for different datareduction-methodes

	@note every derived class has to implement the static functions
	      "DataReducer* create()" and "const String getProductName()" (see FactoryProduct for details)
*/ 
	class DataReducer
	  : public FactoryProduct
	{	
		public:
			typedef MSExperiment<> ExperimentType;
			typedef ExperimentType::SpectrumType SpectrumType;
			typedef ExperimentType::PeakType PeakType;
	
		
			///constructor
			DataReducer();
	
			///destructor
			~DataReducer();
	
			/// assignment operator
			DataReducer& operator = (const DataReducer& source);
	
			///copyconstructor
			DataReducer(const DataReducer& source);
	
	 		/// calculates the reduced MSExperiment
			virtual void applyReduction(const ExperimentType& in, ExperimentType& out )=0;
	
			/// register all derived classes here 
			static void registerChildren();
	
	};
	
}

#endif
